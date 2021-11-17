# %%
from pdf_self_renorm import pdf_zR
from meson_da_hyb_class import MESON_DA_HYB
from head import *

def sum_ft_inverse_gvar(k_ls, fk_ls, delta_k, output_x): # mom to coordinate
    ls_re = []
    ls_im = []
    for idk in range(len(k_ls)):
        ls_re.append( delta_k * (np.exp(-1j * k_ls[idk] * output_x)).real * fk_ls[idk] )
        ls_im.append( delta_k * (np.exp(-1j * k_ls[idk] * output_x)).imag * fk_ls[idk] )
    val_re = np.sum(np.array(ls_re))
    val_im = np.sum(np.array(ls_im))
    return [val_re, val_im]

def test_1(x_ls_matching, lic_da, kernel):
    lcda_tail = []
    for idx in range(len(x_ls_matching)):
        if x_ls_matching[idx] < 0 or x_ls_matching[idx] > 1:
            lcda_tail.append(lic_da[idx])
        else:
            lcda_tail.append(gv.gvar(0, 0))

    inv_kernel = np.linalg.inv(kernel) # lcda to quasi

    quasi = np.dot(inv_kernel, lcda_tail)

    coor_x = np.arange(-100, 100, 0.5) # lambda
    dx = x_ls_matching[1] - x_ls_matching[0]

    coordinate_re = []
    coordinate_im = []
    for y in coor_x:
        coordinate_re.append( sum_ft_inverse_gvar(x_ls_matching, quasi, dx, y)[0] )
        coordinate_im.append( sum_ft_inverse_gvar(x_ls_matching, quasi, dx, y)[1] )

    coordinate_ro = []
    for idx in range(len(coor_x)):
        val = coordinate_re[idx] * (np.exp(1j * coor_x[idx] / 2)).real - coordinate_im[idx] * (np.exp(1j * coor_x[idx] / 2)).imag
        coordinate_ro.append(val)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(x_ls_matching, [(val.mean + val.sdev) for val in lcda_tail], [(val.mean - val.sdev) for val in lcda_tail], color=color_list[1], alpha=0.7, label='Light-cone_tail')
    ax.axvline(0.5, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_ylim([-0.3, 1.6])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(x_ls_matching, [(val.mean + val.sdev) for val in quasi], [(val.mean - val.sdev) for val in quasi], color=color_list[1], alpha=0.7, label='Quasi')
    ax.axvline(0.5, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_ylim([-0.3, 1.6])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(coor_x, [(val.mean + val.sdev) for val in coordinate_re], [(val.mean - val.sdev) for val in coordinate_re], color=color_list[1], alpha=0.7, label='coordinate_re')
    ax.axvline(0, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(coor_x, [(val.mean + val.sdev) for val in coordinate_im], [(val.mean - val.sdev) for val in coordinate_im], color=color_list[1], alpha=0.7, label='coordinate_im')
    ax.axvline(0, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    ###
    meson = 'pion'
    mom = 8
    a06_hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a06_hyb_re_ls')
    a06_hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a06_hyb_im_ls')
    a09_hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a09_hyb_re_ls')
    a09_hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a09_hyb_im_ls')
    a12_hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a12_hyb_re_ls')
    a12_hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a12_hyb_im_ls')

    hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_re_ls') # shape = (N_conf, N_z)
    hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_im_ls')

    pz = int(mom) * 0.215
    lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

    a06_re_ro, a06_im_ro = rotate(a06_hyb_re_ls, a06_hyb_im_ls, lam_ls, back=False)
    a09_re_ro, a09_im_ro = rotate(a09_hyb_re_ls, a09_hyb_im_ls, lam_ls, back=False)
    a12_re_ro, a12_im_ro = rotate(a12_hyb_re_ls, a12_hyb_im_ls, lam_ls, back=False)

    a_re_ro, a_im_ro = rotate(hyb_re_ls, hyb_im_ls, lam_ls, back=False)

    a06_ro_avg = gv.dataset.avg_data(a06_re_ro, bstrap=True)
    a09_ro_avg = gv.dataset.avg_data(a09_re_ro, bstrap=True)
    a12_ro_avg = gv.dataset.avg_data(a12_re_ro, bstrap=True)

    a_ro_avg = gv.dataset.avg_data(a_re_ro, bstrap=True)
    ###


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(coor_x, [(val.mean + val.sdev) for val in coordinate_ro], [(val.mean - val.sdev) for val in coordinate_ro], color=color_list[1], alpha=0.7, label='coordinate_rotated')
    ax.errorbar(lam_ls, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a09_ro_avg], [val.sdev for val in a09_ro_avg], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)
    ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_ro_avg], [(val.mean - val.sdev) for val in a_ro_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_ylim([-0.5, 1.25])
    ax.set_xlim([0, 14])
    ax.set_xlabel(lambda_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()
    
    return

def test_2(x_ls_matching, kernel):
    lcda = []

    # ### normal asymptotic form ###
    # for x in x_ls_matching:
    #     if x>=0 and x<=1:
    #         val = 6 * x * (1 - x) * (1 + 0.1*( 15*(2*x-1)**2 - 3 ) )
    #         lcda.append( val )
    #     else:
    #         lcda.append( 0 )

    ### change the shape of asymptotic form ###
    for x in x_ls_matching:
        if x>=0 and x<=1:
            val = 6 * x * (1 - x) * (1 + 0.1*( 15*(2*x-1)**2 - 3 ) ) + 0.2
            lcda.append( val )
        else:
            val = np.exp(- 5*(x-0.5)**2 )
            lcda.append( val )

    integral = np.sum(lcda) * (x_ls_matching[1] - x_ls_matching[0])

    lcda = np.array(lcda)
    lcda = lcda / integral # normalization

    inv_kernel = np.linalg.inv(kernel) # lcda to quasi

    quasi = np.dot(inv_kernel, lcda)

    dx = x_ls_matching[1] - x_ls_matching[0]

    coor_x = np.arange(-300, 300, 0.5)
    coor_delta = 0.5

    coordinate = []
    for y in coor_x:
        coordinate.append( sum_ft_inv(x_ls_matching, quasi, dx, y) )

    re_quasi = []
    for x in x_ls_matching:
        re_quasi.append( sum_ft(coor_x, coordinate, coor_delta, x) )

    re_lcda = np.dot(kernel, quasi)

    
    # quasi vs lcda
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.plot(x_ls_matching, quasi, color=color_list[1], label='quasi')
    ax.plot(x_ls_matching, lcda, color=color_list[2], label='Light-cone')
    ax.axvline(0.5, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_ylim([-0.3, 1.6])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    # coordinate space 
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.plot(coor_x, [val.real for val in coordinate], color=color_list[1], label='coordinate')
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    # quasi vs re_quasi
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.plot(x_ls_matching, quasi, color=color_list[1], label='quasi')
    ax.plot(x_ls_matching, re_quasi, color='orange', label='re_quasi')
    ax.axvline(0.5, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_ylim([-0.3, 1.6])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()
    
    # lc vs re_lc
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.plot(x_ls_matching, lcda, color=color_list[1], label='Light-cone')
    ax.plot(x_ls_matching, re_lcda, color='orange', label='re_Light-cone')
    ax.axvline(0.5, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_ylim([-0.3, 1.6])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.show()

    return 

# %%
zR_dic, m_pdf_dic = pdf_zR()

meson = 'pion'

x_ls = np.arange(-2-0.01, 3.02, 0.01) # x after ft, for quasi before matching

delta = 0.01

ls_1 = np.arange(-2-0.000001, 0.5-delta, delta)
ls_2 = np.array([0.5])
ls_3 = np.arange(0.5+delta+0.000001, 3+delta, delta)
# ls_1 = np.arange(-2+0.0000001, 0.5, delta)
# ls_2 = np.array([0.5])
# ls_3 = np.arange(0.5+delta-0.0000001, 3, delta)

x_ls_matching = np.hstack((ls_1, ls_2, ls_3)) # x for matching, quasi

ls_4 = np.arange(-2+0.0000001, 0.5, delta)
ls_5 = np.array([0.5])
ls_6 = np.arange(0.5+delta-0.0000001, 3, delta)
y_ls_matching = np.hstack((ls_4, ls_5, ls_6)) # x for matching, light-cone

lambda_ls = z_ls_da * 2*np.pi / (0.0574*96) * 8


extend_point = {}
extend_point['mom=6'] = -2
extend_point['mom=8'] = -5
extend_point['mom=10'] = -8
extend_point['mom=12'] = -8


extend_fit_start = {}
extend_fit_start['mom=6'] = -8
extend_fit_start['mom=8'] = -11
extend_fit_start['mom=10'] = -13
extend_fit_start['mom=12'] = -13

t_dic = {}
t_dic['mom=6'] = {}
t_dic['mom=6']['a06'] = [12, 13, 14]
t_dic['mom=6']['a09'] = [8, 9, 10]
t_dic['mom=6']['a12'] = [6, 7, 8]
t_dic['mom=8'] = {}
t_dic['mom=8']['a06'] = [8, 9, 10]
t_dic['mom=8']['a09'] = [6, 7, 8]
t_dic['mom=8']['a12'] = [4, 5, 6]
t_dic['mom=10'] = {}
t_dic['mom=10']['a06'] = [7, 8, 9]
t_dic['mom=10']['a09'] = [5, 6, 7]
t_dic['mom=10']['a12'] = [3, 4, 5]
t_dic['mom=12'] = {}
t_dic['mom=12']['a06'] = [7, 8, 9]
t_dic['mom=12']['a09'] = [5, 6, 7]
t_dic['mom=12']['a12'] = [3, 4, 5]

pion_mom_8 = MESON_DA_HYB(meson, 8, x_ls, x_ls_matching, y_ls_matching, extend_point['mom=8'], extend_fit_start['mom=8'], t_dic['mom=8'], gs_extract=False, fit_1=False, fit_2=False, rotate=True, lcda_fit=False, constant_fit=False)
mom_8_kernel, mom_8_quasi_da, mom_8_y_ls, mom_8_lic_da = pion_mom_8.main(zR_dic)
kernel = pion_mom_8.matching_kernel()

#test_1(x_ls_matching, mom_8_lic_da, kernel)
test_2(x_ls_matching, kernel)  

# %%
