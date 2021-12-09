# # %%
# from matplotlib.pyplot import title
# from numpy import exp
# from numpy.lib.function_base import select
# from head import *

# print(sc.expi(1.8))

# # %%

# def f_kernel(z2, mu2):
#     res = np.log( z2 * mu2 * np.exp(2 * np.euler_gamma) / 4 )
#     return res

# def matching_kernel_coor():
#     pz = 1.72

#     lam_ls = np.arange(0.06, 1.26, 0.06) # for quasi
#     lamp_ls = np.arange(0.06, 1.26, 0.06) # for light-cone

#     ## delta(1am-lamp) ##
#     delta = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         delta[idx1][idx1] = 1 

#     ## delta * const ##
#     delta_c = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         lam = lam_ls[idx1]
#         z = lam / pz * gev_fm
#         delta_c[idx1][idx1] = alphas_cf_div_2pi/2 * (f_kernel(z**2, mu**2) - 3)

#     ## plus func 1 * coefficient 1 ##
#     p_1 = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1): # lower triangle
#             lam = lam_ls[idx1]
#             lamp = lamp_ls[idx2]
#             p_1[idx1][idx2] = lamp / (lam - lamp)
    
#     for idx1 in range(len(lam_ls)): # diagnoal input
#         if p_1[idx1][idx1] != 0:
#             print('p_1 diagnoal error')
#         p_1[idx1][idx1] = - np.sum( [p_1[i][idx1] for i in range(len(lam_ls))] )

#     coef_1 = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1 + 1): # lower triangle
#             lam = lam_ls[idx1]
#             lamp = lamp_ls[idx2]
#             z = lam / pz * gev_fm
#             coef_1[idx1][idx2] = - ( f_kernel(z**2, mu**2)+1 ) * ( 1 + np.exp(-1j * (lam - lamp)) )

#     Ma_1 = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1 + 1): # lower triangle
#             lam = lam_ls[idx1]
#             Ma_1[idx1][idx2] = p_1[idx1][idx2] * coef_1[idx1][idx2] / lam # variable changing gives a (1/lam)
#     Ma_1 = Ma_1 * alphas_cf_div_2pi

#     ## plus func 2 * coefficient 2 ##
#     p_2 = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1): # lower triangle
#             lam = lam_ls[idx1]
#             lamp = lamp_ls[idx2]
#             p_2[idx1][idx2] = np.log(1 - lamp/lam) / (1 - lamp/lam)
    
#     for idx1 in range(len(lam_ls)): # diagnoal input
#         if p_2[idx1][idx1] != 0:
#             print('p_2 diagnoal error')
#         p_2[idx1][idx1] = - np.sum( [p_2[i][idx1] for i in range(len(lam_ls))] )

#     coef_2 = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1 + 1): # lower triangle
#             lam = lam_ls[idx1]
#             lamp = lamp_ls[idx2]
#             coef_2[idx1][idx2] = -2 * ( 1 + np.exp(-1j * (lam - lamp)) )

#     Ma_2 = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1 + 1): # lower triangle
#             lam = lam_ls[idx1]
#             Ma_2[idx1][idx2] = p_2[idx1][idx2] * coef_2[idx1][idx2] / lam # variable changing gives a (1/lam)
#     Ma_2 = Ma_2 * alphas_cf_div_2pi

#     ## extra ##
#     extra = np.zeros([len(lam_ls), len(lamp_ls)], dtype=complex)
#     for idx1 in range(len(lam_ls)):
#         for idx2 in range(idx1 + 1): # lower triangle
#             lam = lam_ls[idx1]
#             lamp = lamp_ls[idx2]
#             z = lam / pz * gev_fm
#             extra[idx1][idx2] = (1 - np.exp(-1j * (lam - lamp))) / (1j * lam) * (3 - f_kernel(z**2, mu**2))
#             extra[idx1][idx2] = extra[idx1][idx2] / lam
#     extra = extra * alphas_cf_div_2pi


#     C_matrix = delta + delta_c + (Ma_1 + Ma_2 + extra) * (lamp_ls[1] - lamp_ls[0])# multiply by d(lamp) to represent integral
#     C_matrix_inverse = np.linalg.inv(C_matrix)

#     return C_matrix_inverse


# m = matching_kernel_coor()
# print(m)

# %%
from head import *

mom=8
pz = int(mom) * 0.215
lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

lc_re_ls = gv.load('./kaon/mom=8/lc_re_ls')
lc_im_ls = gv.load('./kaon/mom=8/lc_im_ls')
quasi_re_ls = gv.load('./kaon/mom=8/a_hyb_re_ls')
quasi_im_ls = gv.load('./kaon/mom=8/a_hyb_im_ls')

lc_re_ls, lc_im_ls = rotate(lc_re_ls, lc_im_ls, lam_ls, back=False)

re_avg = gv.dataset.avg_data(lc_re_ls, bstrap=True)
im_avg = gv.dataset.avg_data(lc_im_ls, bstrap=True)

q_re_avg = gv.dataset.avg_data(quasi_re_ls, bstrap=True)
q_im_avg = gv.dataset.avg_data(quasi_im_ls, bstrap=True)




#lam_ls = np.linspace(0+0.000001, lambda_ls[-1], 50)

fig = plt.figure(figsize=fig_size)
ax = plt.axes(plt_axes)
ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in re_avg], [(val.mean - val.sdev) for val in re_avg], color='blue', alpha=0.6, label='light-cone')
#ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in q_re_avg], [(val.mean - val.sdev) for val in q_re_avg], color='red', alpha=0.6, label='quasi')
ax.axhline(0, color='k', linestyle='--')
ax.legend(loc='upper right', fontsize=15)
ax.set_xlabel(lambda_label, fontsize=16)
ax.set_ylim([-0.5, 1.25])
ax.set_title('Real part', fontsize=16)
ax.tick_params(direction='in', labelsize=16)
plt.show()

fig = plt.figure(figsize=fig_size)
ax = plt.axes(plt_axes)
ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in im_avg], [(val.mean - val.sdev) for val in im_avg], color='blue', alpha=0.6, label='light-cone')
#ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in q_im_avg], [(val.mean - val.sdev) for val in q_im_avg], color='red', alpha=0.6, label='quasi')
ax.axhline(0, color='k', linestyle='--')
ax.legend(loc='upper right', fontsize=15)
ax.set_xlabel(lambda_label, fontsize=16)
ax.set_ylim([-0.2, 0.2])
ax.set_title('Imag part', fontsize=16)
ax.tick_params(direction='in', labelsize=16)
plt.show()


# %%
