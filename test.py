

# %%
# fig_width = 6.75 # in inches, 2x as wide as APS column
# gr        = 1.618034333 # golden ratio
# fig_size  = (fig_width, fig_width / gr)
# fig_size_lc = (fig_width * 0.8, fig_width * 0.8)
# plt_axes = [0.1,0.12,0.85,0.8]
# fs_p = {"fontsize": 13} # font size of text, label, ticks
# ls_p = {"labelsize": 13}

# errorb = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle

# # labels
# lambda_label = r"$\lambda = z P_z$"
# x_label = r"$x$"
# z_label = r"$z(fm)$"
# t_label = r"$t$"
# meff_local_label = r"$\ln(C(z=0, t) / C(z=0, t+1))$"
# meff_non_local_label = r"$\ln(C(z, t) / C(z, t+1))$"
# ratio_label = r"$C(z, t) / C(z=0, t)$"
# hyb_ro_re_label = r'$Re[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
# hyb_ro_im_label = r'$Im[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
# hyb_re_label = r'$Re[H_{\pi}(z)]$'
# hyb_im_label = r'$Im[H_{\pi}(z)]$'

# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)
# ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in quasi_re_avg], [(val.mean - val.sdev) for val in quasi_re_avg], color='blue', alpha=0.7, label='quasi')
# ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in lc_re_avg], [(val.mean - val.sdev) for val in lc_re_avg], color='red', alpha=0.7, label='lc')
# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-0.5, 1.25])
# ax.tick_params(direction='in', **ls_p)
# plt.show()

# %%
# from head import *

# meson = 'pion'
# mom = 10

# da_an_ls = gv.load(meson+'/mom='+str(mom)+'/da_an_ls')

# print(type(da_an_ls))
# print(np.shape(da_an_ls))

# da_re_ls = np.zeros_like(da_an_ls, dtype=float)
# da_im_ls = np.zeros_like(da_an_ls, dtype=float)

# for i in range(len(da_an_ls)):
#     for j in range(len(da_an_ls[0])):
#         for k in range(len(da_an_ls[0][0])):
#             da_re_ls[i][j][k] = da_an_ls[i][j][k].real
#             da_im_ls[i][j][k] = da_an_ls[i][j][k].imag

# a_ls = [0.0574, 0.0882, 0.1213]
# L_ls = [96, 64, 48]

# for ida in range(3):
#     lam_ls = z_ls_extend * ( 2*np.pi / (a_ls[ida]*L_ls[ida]) * mom * gev_fm ) / gev_fm
#     da_re_ls[ida], da_im_ls[ida] = rotate(da_re_ls[ida], da_im_ls[ida], lam_ls, back=False)

# avg_re = [[], [], []]
# avg_im = [[], [], []]

# for ida in range(3):
#     avg_re[ida] = gv.dataset.avg_data(da_re_ls[ida], bstrap=True)
#     avg_im[ida] = gv.dataset.avg_data(da_im_ls[ida], bstrap=True)
    
# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)
# ax.errorbar(z_ls_extend, [val.mean for val in avg_re[0]], [val.sdev for val in avg_re[0]], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
# ax.errorbar(z_ls_extend, [val.mean for val in avg_re[1]], [val.sdev for val in avg_re[1]], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
# ax.errorbar(z_ls_extend, [val.mean for val in avg_re[2]], [val.sdev for val in avg_re[2]], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)
# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-0.5, 1.25])
# ax.tick_params(direction='in', **ls_p)
# ax.set_title(hyb_ro_re_label, **fs_p)

# plt.show()


# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)
# ax.errorbar(z_ls_extend, [val.mean for val in avg_im[0]], [val.sdev for val in avg_im[0]], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
# ax.errorbar(z_ls_extend, [val.mean for val in avg_im[1]], [val.sdev for val in avg_im[1]], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
# ax.errorbar(z_ls_extend, [val.mean for val in avg_im[2]], [val.sdev for val in avg_im[2]], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)
# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-1.2, 0.5])
# ax.tick_params(direction='in', **ls_p)
# ax.set_title(hyb_ro_im_label, **fs_p)

# plt.show()

# %%
# from head import *

# x_ls = gv.load('temp/x_ls')
# lc_8 = gv.load('temp/lc_-8')
# lc_10 = gv.load('temp/lc_-10')
# lc_13 = gv.load('temp/lc_-13')
# lc_15 = gv.load('temp/lc_-15')

# qu_8 = gv.load('temp/quasi_-8')
# qu_10 = gv.load('temp/quasi_-10')
# qu_13 = gv.load('temp/quasi_-13')
# qu_15 = gv.load('temp/quasi_-15')

# qu_16 = gv.load('temp/quasi_-16')

# fig = plt.figure(figsize=fig_size_lc)
# ax = plt.axes(plt_axes)
# label_list = [r'$\lambda=6.57$', r'$\lambda=5.75$',  r'$\lambda=4.52$', r'$\lambda=3.69$']
# i=0
# #[lc_8, lc_10, lc_13, lc_15]:#
# for qu in [qu_16]:
#     ax.fill_between(x_ls, [val.mean+val.sdev for val in qu], [val.mean-val.sdev for val in qu], alpha=0.6, color=color_list[i], label=label_list[i])
#     i = i + 1
# ax.plot(x_ls, [val.mean for val in qu], linewidth=2, linestyle='dotted', color='red')
# ax.axvline(0, color='k', linestyle='--')
# ax.axvline(0.5, color='green', linestyle='--')
# ax.axhline(0, color='k', linestyle='--')
# ax.set_xlabel(x_label, **fs_p)
# ax.set_ylim([-0.19, 1.7])
# ax.set_xlim([-1.5, 2.5])
# #ax.legend(loc='upper right')
# ax.tick_params(direction='in', **ls_p)
# plt.show()
# %%
# from head import *

# x_ls = gv.load('temp/x_ls')
# ext_po = gv.load('temp/ext_po')
# ext_ex200 = gv.load('temp/ext_ex200')
# ext_ex150 = gv.load('temp/ext_ex150')
# ext_ex50 = gv.load('temp/ext_ex50')
# ext_ex20 = gv.load('temp/ext_ex20')


# fig = plt.figure(figsize=fig_size_lc)
# ax = plt.axes(plt_axes)
# label_list = ['Polynomial', r'$\lambda_0=200$',  r'$\lambda_0=150$', r'$\lambda_0=50$', r'$\lambda_0=20$']
# i=0
# #[lc_8, lc_10, lc_13, lc_15]:#
# for lc in [ext_po, ext_ex200, ext_ex150, ext_ex50, ext_ex20]:
#     ax.fill_between(x_ls, [val.mean+val.sdev for val in lc], [val.mean-val.sdev for val in lc], alpha=0.4, color=color_list[2*i], label=label_list[i])
#     i = i + 1
#     #ax.plot(x_ls, [val.mean for val in lc], linewidth=2, linestyle='dashed', color=color_list[i])
# ax.axvline(0, color='k', linestyle='--')
# ax.axvline(0.5, color='green', linestyle='--')
# ax.axhline(0, color='k', linestyle='--')
# ax.set_xlabel(x_label, **fs_p)
# ax.set_ylim([-0.19, 1.5])
# ax.set_xlim([-0.5, 1.5])
# ax.legend(loc='upper right')
# ax.tick_params(direction='in', **ls_p)
# plt.show()


# # %%
# from head import *

# def val(px):
#     return 0.2 + (px - 495)/(844 - 495) * 0.1


# fig = plt.figure(figsize=fig_size_lc)
# plt_axes_ = [0.25,0.15,0.7,0.8]
# ax = plt.axes(plt_axes_)

# ax.errorbar([val(790)], [1], xerr=[val(980) - val(790)], fmt='D', **errorb)
# ax.errorbar([val(774)], [2], xerr=[val(844) - val(774)], fmt='D', **errorb)
# ax.errorbar([val(620)], [3], xerr=[val(646) - val(620)], fmt='D', **errorb)
# ax.errorbar([val(648)], [4], xerr=[val(754) - val(648)], fmt='D', **errorb)
# ax.errorbar([val(531)], [5], xerr=[val(576) - val(531)], fmt='D', color='brown', **errorb)
# ax.errorbar([val(531)], [5], xerr=[val(691) - val(531)], fmt='D', color='brown', **errorb)
# ax.errorbar([val(531)], [5.5], xerr=[val(664) - val(531)], color='violet', fmt='D', **errorb)
# ax.errorbar([val(531)], [5.5], xerr=[val(576) - val(531)], fmt='D', color='violet', **errorb)
# ax.errorbar([0.303], [6.5], xerr=[0.029], fmt='D', color='purple', **errorb)


# ax.axvline(0.303-0.029, color='purple', linestyle='--', lw=0.5)
# ax.axvline(0.303+0.029, color='purple', linestyle='--', lw=0.5)
# ax.axvline(0.303-2*0.029, color='purple', linestyle='dashdot', lw=0.5)

# ax.set_xlabel(x_label, **fs_p)
# ax.set_ylim([0.5, 7])
# ax.set_xlim([0.15, 0.35])
# ax.tick_params(direction='in', labelsize=12)
# ax.set_xticks([0.2, 0.3])
# ax.set_xlabel(r'$<\xi^2>(\mu = 2 \rm{GeV})$')
# plt.yticks([1,2,3,4,5,5.5, 6.5], ['Del Debbio et al.\n(2003)', 'Arthur et al.\n(2011)', 'Bali et al.\n(2019)', 'Zhang et al.\n(2020)', 'HOPE Mom', 'HOPE TMR', 'This work'], fontsize=10)
# plt.show()

# %%
from head import *
x_ls = np.linspace(-0.5, 1.5, 500)
y_ls = DSE(x_ls)

a2 = calc_an(x_ls, y_ls, 2)
print(a2)
# %%
