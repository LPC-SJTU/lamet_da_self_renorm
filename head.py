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
import timeit


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
fig_size_lc = (fig_width * 0.8, fig_width * 0.6)
fig_size_sq = (fig_width * 0.8, fig_width * 0.8)
gridspec_sub = {'height_ratios': [3, 1], 'left': 0.12, 'right': 0.99, 'bottom': 0.15, 'top': 0.98}
plt_axes = [0.1,0.12,0.85,0.8]
plt_axes_small = [0.12,0.15,0.8,0.75]
errorp = {"markersize": 5, "mfc": "none", "linestyle": "none"} # circle
errorb = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle
errorl = {"markersize": 5, "mfc": "none", "capsize": 3, "elinewidth": 1} # circle with line
fs_p = {"fontsize": 13} # font size of text, label, ticks
fs_p_l = {"fontsize": 19} 
ls_p = {"labelsize": 13}
ls_p_l = {"labelsize": 15.5}

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
mu = 2 # GeV, for renormalization #!# for sys err
mu_f = 2 # GeV, for factorization

cf=4/3
nf=3
b0=11-2/3*nf
alphas=2 * np.pi/( b0 * np.log(mu_f/lms) )

mom_to_pz = 0.215 # mom=8 corresponding to pz=1.72, pz=2pi / (0.09*64) * mom=8 * 0.197
alphas_cf_div_2pi = alphas * cf / (2 * np.pi)

z_ls = np.arange(0.06, 1.26, 0.06) # read bare pdf then interpolated into z_ls
z_ls_extend = np.arange(0.06, 1.56, 0.06) # extend f1 and zR
z_ls_da = np.arange(0.06, 1.39, 0.06) # for interpolation of lambda
zs_pz = 2 # for Z_hyb
extend_length_quasi = 200 # how long for extension before FT, 200 about till lambda=120
extend_length_lc = 1000 # how long for extension before FT, 200 about till lambda=120

## lcda endpoints fit
fit_x1 = 0.1 # [0, x1] and [x4, 1] use fit 
fit_x2 = 0.2 # [x1, x2] and [x3, x4] are fit region
fit_x3 = 0.8 # [x2, x3] keeps originally
fit_x4 = 0.9



grey = "#808080" 
red = "#FF6F6F" 
peach = "#FF9E6F" 
orange = "#FFBC6F" 
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C" 
turquoise = "#4AAB89"
blue = "#508EAD" 
grape = "#635BB1"
violet = "#7C5AB8" 
fuschia = "#C3559F"

#color_list = [orange, blue, violet, red, turquoise, lime, peach, sunkist, green, grey, grape]
color_list = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green']
fmt_ls = ['o', 'D', 's', '^', 'x', '.']



# labels
lambda_label = r"$\lambda = z P_z$"
x_label = r"$x$"
z_label = r"$z(fm)$"
t_label = r"$t$"
meff_local_label = r"$\ln(C(z=0, t) / C(z=0, t+1))$"
meff_non_local_label = r"$\ln(C(z, t) / C(z, t+1))$"
ratio_label = r"$C(z, t) / C(z=0, t)$"
hyb_ro_re_label = r'Re$[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
hyb_ro_im_label = r'Im$[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
hyb_re_label = r'$Re[H_{\pi}(z)]$'
hyb_im_label = r'$Im[H_{\pi}(z)]$'

def interp_1d(x_in, y_in, x_out, method="linear"): # interpolation
    f=interpolate.interp1d(x_in, y_in, kind=method)
    y_out = f(x_out)

    return y_out

# def bootstrap(conf_ls, N_re):
#     N_conf = len(conf_ls)
#     conf_re = []
#     for times in range(N_re):
#         idx_ls = np.random.randint(N_conf, size=N_conf)
#         temp = []
#         for idx in idx_ls:
#             temp.append(conf_ls[idx])
#         conf_re.append( np.average(temp, axis=0) )

#     return np.array(conf_re)

def bootstrap(conf_ls, meson, a_str):
    N_conf = len(conf_ls)
    conf_re = []
    bs_ls = gv.load('bs_ls')[meson+a_str]
    for times in range(len(bs_ls)):
        idx_ls = bs_ls[times]
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
    x_ls = np.array(x_ls)
    fx_ls = np.array(fx_ls)
    val = delta_x/(2*np.pi) * np.sum( np.exp(1j * x_ls * output_k) * fx_ls )

    return val

def sum_ft_inv(k_ls, fk_ls, delta_k, output_x): # momentum to coordinate
    k_ls = np.array(k_ls)
    fk_ls = np.array(fk_ls)
    val = delta_k * np.sum( np.exp(-1j * k_ls * output_x) * fk_ls )

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

def DSE(x):
    return 18.2 * x * (1-x) * ( 1 - 2.33*np.sqrt(x*(1-x)) + 1.79*x*(1-x) )

def DSE_kaon():
    pix_x = [22,66,163,235,307,402,429,547,669,829,906,960,1033,982,1112,1185,1266,1306,1373,1417,1454,1486,1509,1525]
    pix_y = [75,204,433,557,658,751,770,840,880,906,907,905,891,904,864,821,745,692,571,469,364,249,161,80]

    x1 = 1537 # coor of x=1
    y1 = 685 # coor of y=1 

    x = []
    y = []

    for val in pix_x:
        x.append( (x1-val)/x1 )
    
    for val in pix_y:
        y.append( val/y1 )

    x = np.array(x)
    y = np.array(y)

    return x, y

def mellin_moment(x_ls, lc_ls, n):
    x_array = np.array(x_ls)
    zeta_array = 2 * x_array - 1

    ### normalization check ###
    val = []
    for idx in range(len(zeta_array)):
        if zeta_array[idx] >=-1 and zeta_array[idx] <=1:
            val.append(lc_ls[idx].mean)
    print('>>> normalization check:')
    print( np.sum(val) * (zeta_array[1]-zeta_array[0] ) * 1/2 ) # 1/2 because zeta has twice bigger axis than x

    ### n order moment ###
    val = []
    for idx in range(len(zeta_array)):
        if zeta_array[idx] >=-1 and zeta_array[idx] <=1:
            val.append( zeta_array[idx]**n * lc_ls[idx] )

    print('>>> mellin moment order '+str(n)+' :')
    print( np.sum(val) * (zeta_array[1]-zeta_array[0]) * 1/2)

    return

def endpoint_ext(x_ls, lc, meson): # dtype of lc should be gvar
    x_1 = []
    x_2 = []
    x_3 = []

    lc_1 = []
    lc_2 = []
    lc_3 = []

    for idx in range(len(x_ls)):
        x = x_ls[idx]
        y = lc[idx]
        if 0 <= x <= fit_x1 or fit_x4 <= x <= 1: # use fit
            x_1.append(x)
        elif fit_x1 <= x <= fit_x2 or fit_x3 <= x <= fit_x4: # fit region
            x_2.append(x)
            lc_2.append(y)
        elif fit_x2 <= x <= fit_x3: # keep same 
            x_3.append(x)
            lc_3.append(y.mean)

    if meson == 'pion':
        def fcn(x, p):
            return p['c'] * ( x**p['n'] ) * ( (1-x)**p['n'] )

    elif meson == 'kaon':
        def fcn(x, p):
            return p['c'] * ( x**p['n1'] ) * ( (1-x)**p['n2'] )

    priors = gv.BufferDict()
    priors['c'] = gv.gvar(0, 10)
    priors['n'] = gv.gvar(0, 10)
    priors['n1'] = gv.gvar(0, 10)
    priors['n2'] = gv.gvar(0, 10)

    fit_result = lsf.nonlinear_fit(data=(np.array(x_2), lc_2), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    x_ls_new = np.hstack( (np.array(x_1), np.array(x_2), np.array(x_3)) )
    lc_1 = [val.mean for val in fcn(np.array(x_1), fit_result.p)]
    lc_2 = [val.mean for val in lc_2]
    lc_new = np.hstack( (lc_1, lc_2, lc_3) )

    z = zip(x_ls_new, lc_new)
    z = sorted(z)
    x_ls_new, lc_new = zip(*z)

    x_ls_new = np.array(x_ls_new)

    return x_ls_new, lc_new