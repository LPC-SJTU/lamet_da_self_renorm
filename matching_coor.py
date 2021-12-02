# %%
import numpy as np
import gvar as gv
from scipy import interpolate
from scipy import integrate
from tqdm import tqdm
from head import *

mu = 2 # GeV, for renormalization

cf=4/3
nf=3
b0=11-2/3*nf
lms = 0.24451721864451428  
alphas=2 * np.pi/( b0 * np.log(mu/lms) )

alphas_cf_div_2pi = alphas * cf / (2 * np.pi)

def f_matching(z2, mu2):
    res = np.log( z2 * mu2 * np.exp(2*np.euler_gamma) / 4 )
    return res

def inv_matching_coor(lambda_ls, quasi_re_ls, quasi_im_ls): # the shape of quasi_re/im_ls is [N_conf, N_lambda]
    pz = 2.15

    lc_re_ls = []
    lc_im_ls = []

    delta_l = 0.000001 # regulater

    print('>>> integrating for matching')
    for n_conf in tqdm(range(len(quasi_re_ls))):
        lc_re_ls.append([])
        lc_im_ls.append([])

        quasi = [1+0j] # add the point at lambda=0
        for idl in range(len(quasi_re_ls[n_conf])):
            quasi.append( quasi_re_ls[n_conf][idl] + 1j * quasi_im_ls[n_conf][idl] )

        lam_ls_0 = np.insert(lambda_ls, 0, 0) # add the point at lambda=0
        h_tilde = interpolate.interp1d(lam_ls_0, quasi, kind='cubic') # analytic function of quasi

        for idl in range(len(lambda_ls)):
            lam = lambda_ls[idl]
            z = lam / pz # here z has GeV^-1

            lc_re = quasi_re_ls[n_conf][idl] #
            lc_im = quasi_im_ls[n_conf][idl]

            part1 = alphas_cf_div_2pi * 1/2 * (f_matching(z**2, mu**2) - 3) * h_tilde(lam)

            def fp_2_re(lamp, lam):
                z = lam / pz # here z has GeV^-1
                res = alphas_cf_div_2pi / lam * (-1-f_matching(z**2, mu**2)) * (lamp/(lam-lamp)) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                    - 2 * h_tilde(lam)
                )
                return res.real

            def fp_2_im(lamp, lam):
                z = lam / pz # here z has GeV^-1
                res = alphas_cf_div_2pi / lam * (-1-f_matching(z**2, mu**2)) * (lamp/(lam-lamp)) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                    - 2 * h_tilde(lam)
                )
                return res.imag

            part2 = integrate.quad(fp_2_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_2_im, 0, lam-delta_l, args=lam)[0] # [0] is result, [1] is error estimation

            def fp_3_re(lamp, lam):
                res = alphas_cf_div_2pi / lam * np.log(1 - lamp/lam) / (1 - lamp/lam) * (-2) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                    - 2 * h_tilde(lam)
                )
                return res.real
            
            def fp_3_im(lamp, lam):
                res = alphas_cf_div_2pi / lam * np.log(1 - lamp/lam) / (1 - lamp/lam) * (-2) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                    - 2 * h_tilde(lam)
                )
                return res.imag

            part3 = integrate.quad(fp_3_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_3_im, 0, lam-delta_l, args=lam)[0]

            def fp_4_re(lamp, lam):
                z = lam / pz # here z has GeV^-1
                res = alphas_cf_div_2pi / (1j * (lam**2)) * ( 1 - np.exp(-1j * (lam - lamp)) ) * (3-f_matching(z**2, mu**2)) * h_tilde(lamp)
                return res.real

            def fp_4_im(lamp, lam):
                z = lam / pz # here z has GeV^-1
                res = alphas_cf_div_2pi / (1j * (lam**2)) * ( 1 - np.exp(-1j * (lam - lamp)) ) * (3-f_matching(z**2, mu**2)) * h_tilde(lamp)
                return res.imag
            part4 = integrate.quad(fp_4_re, 0, lam, args=lam)[0] + 1j * integrate.quad(fp_4_im, 0, lam, args=lam)[0]

            lc_re = lc_re - (part1+part2+part3+part4).real
            lc_re_ls[n_conf].append(lc_re)

            lc_im = lc_im - (part1+part2+part3+part4).imag
            lc_im_ls[n_conf].append(lc_im)

    return lc_re_ls, lc_im_ls
# %%
meson = 'pion'
mom = 10
pz = int(mom) * 0.215
lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_re_ls') # shape = (N_conf, N_z)
hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_im_ls')

hyb_re_ls, hyb_im_ls = rotate(hyb_re_ls, hyb_im_ls, lam_ls, back=True)

# gv.dump(lam_ls, 'lam_ls')
# gv.dump(hyb_re_ls, 'quasi_ro_re_ls')
# gv.dump(hyb_im_ls, 'quasi_ro_im_ls')

quasi_re_avg = gv.dataset.avg_data(hyb_re_ls, bstrap=True)
quasi_im_avg = gv.dataset.avg_data(hyb_im_ls, bstrap=True)

quasi_re = [[val.mean for val in quasi_re_avg]]
quasi_im = [[val.mean for val in quasi_im_avg]]

lc_re, lc_im = inv_matching_coor(lam_ls, quasi_re, quasi_im)

lc_re_ls, lc_im_ls = inv_matching_coor(lam_ls, hyb_re_ls, hyb_im_ls)

lc_re_avg = gv.dataset.avg_data(lc_re_ls, bstrap=True)
lc_im_avg = gv.dataset.avg_data(lc_im_ls, bstrap=True)

f = open('./lc.txt', 'w')
line = []
line.append('z')
line.append('\t')
line.append('lambda')
line.append('\t')
line.append('lc_Re')
line.append('\t')
line.append('lc_Re_err')
line.append('\t')
line.append('lc_Im')
line.append('\t')
line.append('lc_Im_err')
line.append('\n')
f.writelines(line)

for idz in range(len(z_ls_da)):
    line = []
    line.append(str(z_ls_da[idz]))
    line.append('\t')
    line.append(str(lam_ls[idz]))
    line.append('\t')
    line.append(str(lc_re_avg[idz].mean))
    line.append('\t')
    line.append(str(lc_re_avg[idz].sdev))
    line.append('\t')
    line.append(str(lc_im_avg[idz].mean))
    line.append('\t')
    line.append(str(lc_im_avg[idz].sdev))
    line.append('\n')
    f.writelines(line)
f.close

gv.dump(lc_re_avg, 'lc_re')
gv.dump(lc_im_avg, 'lc_im')
gv.dump(lam_ls, 'lam_ls')

# %%
# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)
# ax.scatter(lam_ls, lc_re[0], color='red', label='lc')
# ax.scatter(lam_ls, quasi_re[0], color='blue', label='quasi')
# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-0.5, 1.5])
# ax.tick_params(direction='in', **ls_p)
# plt.show()

# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)
# ax.scatter(lam_ls, lc_im[0], color='red', label='lc')
# ax.scatter(lam_ls, quasi_im[0], color='blue', label='quasi')
# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-1.2, 0.5])
# ax.tick_params(direction='in', **ls_p)
# plt.show()

# f = open('./match_mean.txt', 'w')
# line = []
# line.append('z')
# line.append('\t')
# line.append('lambda')
# line.append('\t')
# line.append('quasi_Re')
# line.append('\t')
# line.append('quasi_Im')
# line.append('\t')
# line.append('lc_Re')
# line.append('\t')
# line.append('lc_Im')
# line.append('\n')
# f.writelines(line)

# for idz in range(len(z_ls_da)):
#     line = []
#     line.append(str(z_ls_da[idz]))
#     line.append('\t')
#     line.append(str(lam_ls[idz]))
#     line.append('\t')
#     line.append(str(quasi_re[0][idz]))
#     line.append('\t')
#     line.append(str(quasi_im[0][idz]))
#     line.append('\t')
#     line.append(str(lc_re[0][idz]))
#     line.append('\t')
#     line.append(str(lc_im[0][idz]))
#     line.append('\n')
#     f.writelines(line)
# f.close



