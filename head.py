from unicodedata import normalize
import h5py as h5
from h5py._hl import dataset
import lsqfit as lsf
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
import os
import multiprocessing


from multiprocessing.pool import ThreadPool
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import norm
from numpy.core.fromnumeric import shape
from numpy.lib.stride_tricks import as_strided
from mpmath import *
from scipy import interpolate
from tqdm import tqdm
from scipy.optimize import curve_fit
from scipy import special as sp
from scipy.misc import derivative as drv
from scipy import integrate

fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
fig_size_lc = (fig_width * 0.8, fig_width * 0.8)
gridspec_sub = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}
plt_axes = [0.1,0.12,0.85,0.8]
plt_axes_small = [0.12,0.15,0.8,0.75]
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none"} # circle
errorb = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle
errorl = {"markersize": 5, "mfc": "none", "capsize": 3, "elinewidth": 1} # circle with line
fs_p = {"fontsize": 13} # font size of text, label, ticks
ls_p = {"labelsize": 13}

gev_fm = 0.1973269631 # 1 = 0.197 GeV . fm #
delta_sys = 0 # 0.002, used to add sys error for pdf data after interpolation

a_milc_ls = [0.1213, 0.0882, 0.0574, 0.0425, 0.0318]
a_rbc_ls = [0.11, 0.0828, 0.0626]

lqcd = 0.1 #Lambda_QCD
lms = 0.24451721864451428  #Lambda_MS
k = 3.320
d_pdf = -0.08183 #-0.1252 # for pdf
d_da = 0.19 # 0.1 for a^2 order, 0.19 for a order, from Yushan
m0_da = gv.gvar(-0.094, 0.024) # from Yushan
mu = 2 # GeV, for renormalization
mu_f = 2 # GeV, for factorization

cf=4/3
nf=3
b0=11-2/3*nf
alphas=2 * np.pi/( b0 * np.log(mu/lms) )

mom_to_pz = 0.215 # mom=8 corresponding to pz=1.72, pz=2pi / (0.09*64) * mom=8 * 0.197
alphas_cf_div_2pi = alphas * cf / (2 * np.pi)

z_ls = np.arange(0.06, 1.26, 0.06) # read bare pdf then interpolated into z_ls
z_ls_extend = np.arange(0.06, 1.56, 0.06) # extend f1 and zR
z_ls_da = np.arange(0.06, 1.39, 0.06) # for interpolation of lambda
zs_pz = 2 # for Z_hyb
extend_length_quasi = 200 # how long for extension before FT, 200 about till lambda=120
extend_length_lc = 1000 # how long for extension before FT, 200 about till lambda=120

## lcda endpoints fit
fit_y1 = 0.05 # [0, y1] and [y4, 1] use fit 
fit_y2 = 0.2 # [y1, y2] and [y3, y4] are fit region
fit_y3 = 0.8 # [y2, y3] keeps originally
fit_y4 = 0.95


color_list = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green']
fmt_ls = ['o', 'D', 's', '^']

# labels
lambda_label = r"$\lambda = z P_z$"
x_label = r"$x$"
z_label = r"$z(fm)$"
t_label = r"$t$"
meff_local_label = r"$\ln(C(z=0, t) / C(z=0, t+1))$"
meff_non_local_label = r"$\ln(C(z, t) / C(z, t+1))$"
ratio_label = r"$C(z, t) / C(z=0, t)$"
hyb_ro_re_label = r'$Re[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
hyb_ro_im_label = r'$Im[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
hyb_re_label = r'$Re[H_{\pi}(z)]$'
hyb_im_label = r'$Im[H_{\pi}(z)]$'

def interp_1d(x_in, y_in, x_out, method="linear"): # interpolation
    f=interpolate.interp1d(x_in, y_in, kind=method)
    y_out = f(x_out)

    return y_out

def bootstrap(conf_ls, N_re):
    N_conf = len(conf_ls)
    conf_re = []
    for times in range(N_re):
        idx_ls = np.random.randint(N_conf, size=N_conf)
        temp = []
        for idx in idx_ls:
            temp.append(conf_ls[idx])
        conf_re.append( np.average(temp, axis=0) )

    return np.array(conf_re)

def ZMS_pdf(z):
    #z fm
    ans = 1 + alphas*cf/(2*np.pi)*( 3/2* np.log(mu**2 * (z/gev_fm)**2 * np.exp(2*np.euler_gamma) / 4) + 5/2 ) #!#
    return ans

def ZMS_da(z):
    #z fm
    ans = 1 + alphas*cf/(2*np.pi)*( 3/2* np.log(mu**2 * (z/gev_fm)**2 * np.exp(2*np.euler_gamma) / 4) + 7/2 ) #!#
    return ans

def turn_into_lambda(hyb_conf_ls, z_ls, lambda_ls, a, L, mom):
    pz = 2*np.pi / (a * L) * mom * gev_fm # pz=2pi / (0.09*64) * mom=8 * 0.197
    lambda_in = z_ls * pz / gev_fm

    hyb_conf_ls_lambda = []
    for n_conf in range(len(hyb_conf_ls)):
        hyb_conf_ls_lambda.append( interp_1d(lambda_in, hyb_conf_ls[n_conf], lambda_ls, method='cubic') )

    return hyb_conf_ls_lambda

def add_sdev(hyb_ls, hyb_avg):
    hyb_conf = hyb_ls # shape = len(z_ls_extend)
    hyb_conf_gv = []
    for idx in range(len(hyb_ls)): # each conf add the same sdev
        hyb_conf_gv.append( gv.gvar(hyb_conf[idx], hyb_avg[idx].sdev) )

    return hyb_conf_gv

def f_matching(z2, mu2):
    res = np.log( z2 * mu2 * np.exp(2*np.euler_gamma) / 4 )
    return res

def sum_ft(x_ls, fx_ls, delta_x, output_k): # coordinate to momentum
    ls = []
    for idx in range(len(x_ls)):
        ls.append( delta_x/(2*np.pi) * np.exp(1j * x_ls[idx] * output_k) * fx_ls[idx] )
    val = np.sum(np.array(ls))
    return val

def sum_ft_inv(k_ls, fk_ls, delta_k, output_x): # momentum to coordinate
    ls = []
    for idx in range(len(k_ls)):
        ls.append( delta_k * np.exp(-1j * k_ls[idx] * output_x) * fk_ls[idx] )
    val = np.sum(np.array(ls))
    return val

def rotate(re_conf_z, im_conf_z, lam_ls, back): # 2D list, 0 for configs, 1 for z_ls
    re_ro_conf_z = []
    im_ro_conf_z = []
    for n_conf in range(len(re_conf_z)):
        re_ro_conf_z.append([])
        im_ro_conf_z.append([])
        for idx in range(len(re_conf_z[0])):
            val = complex(re_conf_z[n_conf][idx], im_conf_z[n_conf][idx])
            val = val * np.exp(1j * lam_ls[idx] / 2)
            if back == True: # whether rotate back
                re_ro_conf_z[n_conf].append( (val.real * np.exp(-1j * lam_ls[idx] / 2)).real )
                im_ro_conf_z[n_conf].append( (val.real * np.exp(-1j * lam_ls[idx] / 2)).imag )
            elif back == False:
                re_ro_conf_z[n_conf].append(val.real)
                im_ro_conf_z[n_conf].append(val.imag)

    return re_ro_conf_z, im_ro_conf_z

def gen_an(x, n):
    factor = 2*(2*n+3) / (3*(n+1)*(n+2))
    if n==1:
        return factor*3*(2*x-1)
    if n==2:
        return factor*3/2*(5*(2*x-1)**2-1)
    if n==4:
        return factor*15/8*(1 - 14*(2*x-1)**2 + 21*(2*x-1)**4)

def calc_an(y_ls, lcda, n):
    val = []
    for idx in range(len(y_ls)):
        if y_ls[idx] >= 0 and y_ls[idx] <= 1:
            temp = lcda[idx] * gen_an(y_ls[idx], n) # an
            val.append(temp)

    return np.sum(val)*(y_ls[1]-y_ls[0])

def norm_check(y_ls, lcda):
    val = []
    for idx in range(len(y_ls)):
        if y_ls[idx] >=0 and y_ls[idx] <=1:
            val.append(lcda[idx].mean)

    return np.sum(val)*(y_ls[1]-y_ls[0])

def large_mom_limit(y_ls, mom_n1_lic_da, mom_n2_lic_da, mom_n3_lic_da, mom_ls, meson):
    large_mom_lic_da = []

    for idx in range(len(y_ls)):
        def fcn(x, p):
            return p['psi'] + p['c2']/(x**2)

        priors = gv.BufferDict()
        priors['psi'] = gv.gvar(0.5, 2)
        priors['c2'] = gv.gvar(1, 10)

        pz_ls = np.array(mom_ls) * mom_to_pz
        lcda_ls = [ mom_n1_lic_da[idx], mom_n2_lic_da[idx], mom_n3_lic_da[idx] ]

        fit_result = lsf.nonlinear_fit(data=(pz_ls, lcda_ls), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

        large_mom_lic_da.append(fit_result.p['psi'])

    # mom=inf
    if meson == 'kaon':
        print('mom=inf, a1: ')
        a1 = calc_an(y_ls, large_mom_lic_da, 1)
        print(a1)
        
    print('mom=inf, a2: ')
    a2 = calc_an(y_ls, large_mom_lic_da, 2)
    print(a2)

    print('Light-cone at inf pz integral within [0, 1]: ')
    normal = norm_check(y_ls, large_mom_lic_da) 
    print(normal)

    return large_mom_lic_da

def sum_rule(meson, x, a1, a2, a4):
    def C1(x, a):
        return 2*a*x
    def C2(x, a):
        return 1/2 * ( 2 * x * (2+a-1) * C1(x, a) - (2+2*a-2) )
    def C3(x, a):
        return 1/3 * ( 2 * x * (3+a-1) * C2(x, a) - (3+2*a-2) * C1(x, a) )
    def C4(x, a):
        return 1/4 * ( 2 * x * (4+a-1) * C3(x, a) - (4+2*a-2) * C2(x, a) )

    if meson == 'pion':
        return 6 * x * (1-x) * ( 1 + C2(2*x-1, 3/2)*a2 + C4(2*x-1, 3/2)*a4 )

    if meson == 'kaon':
        return 6 * x * (1-x) * ( 1 + C1(2*x-1, 3/2)*a1 + C2(2*x-1, 3/2)*a2 )

