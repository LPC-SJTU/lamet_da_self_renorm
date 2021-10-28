# %%
from head import *

def paper_plot_discrete_effect(mom, meson='pion', if_rotate=True):
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

    if if_rotate == True:
        a06_re_ro, a06_im_ro = rotate(a06_hyb_re_ls, a06_hyb_im_ls, lam_ls, back=False)
        a09_re_ro, a09_im_ro = rotate(a09_hyb_re_ls, a09_hyb_im_ls, lam_ls, back=False)
        a12_re_ro, a12_im_ro = rotate(a12_hyb_re_ls, a12_hyb_im_ls, lam_ls, back=False)

        a_re_ro, a_im_ro = rotate(hyb_re_ls, hyb_im_ls, lam_ls, back=False)

    elif if_rotate == False:
        a06_re_ro, a06_im_ro = a06_hyb_re_ls, a06_hyb_im_ls
        a09_re_ro, a09_im_ro = a09_hyb_re_ls, a09_hyb_im_ls
        a12_re_ro, a12_im_ro = a12_hyb_re_ls, a12_hyb_im_ls

        a_re_ro, a_im_ro = hyb_re_ls, hyb_im_ls

    ## real part
    a06_ro_avg = gv.dataset.avg_data(a06_re_ro, bstrap=True)
    a09_ro_avg = gv.dataset.avg_data(a09_re_ro, bstrap=True)
    a12_ro_avg = gv.dataset.avg_data(a12_re_ro, bstrap=True)

    a_ro_avg = gv.dataset.avg_data(a_re_ro, bstrap=True)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(lam_ls, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a09_ro_avg], [val.sdev for val in a09_ro_avg], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)

    ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_ro_avg], [(val.mean - val.sdev) for val in a_ro_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

    ax.axhline(0, color='k', linestyle='--')
    ax.legend(loc='upper right')
    ax.set_xlabel(lambda_label, **fs_p)
    ax.set_ylim([-0.5, 1.25])
    ax.tick_params(direction='in', **ls_p)
    if if_rotate == True:
        plt.savefig(meson+'/paper/discrete_effect_rotated_Pz='+str(pz)+'GeV.pdf', transparent=True)
    elif if_rotate == False:
        plt.savefig(meson+'/paper/discrete_effect_Pz='+str(pz)+'GeV.pdf', transparent=True)
    plt.show()

    ## imag part
    a06_ro_avg = gv.dataset.avg_data(a06_im_ro, bstrap=True)
    a09_ro_avg = gv.dataset.avg_data(a09_im_ro, bstrap=True)
    a12_ro_avg = gv.dataset.avg_data(a12_im_ro, bstrap=True)

    a_ro_avg = gv.dataset.avg_data(a_im_ro, bstrap=True)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(lam_ls, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a09_ro_avg], [val.sdev for val in a09_ro_avg], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)

    ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_ro_avg], [(val.mean - val.sdev) for val in a_ro_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

    ax.axhline(0, color='k', linestyle='--')
    ax.legend(loc='upper right')
    ax.set_xlabel(lambda_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)
    if if_rotate == True:
        ax.set_ylim([-0.3, 0.3])
        plt.savefig(meson+'/paper/discrete_effect_rotated_imag_part_Pz='+str(pz)+'GeV.pdf', transparent=True)
    elif if_rotate == False:
        ax.set_ylim([-1, 0.5])
        plt.savefig(meson+'/paper/discrete_effect_imag_part_Pz='+str(pz)+'GeV.pdf', transparent=True)
    plt.show()
    return 

def quasi_vs_lc_plot(x_ls, y_ls, quasi_da, lic_da, pz, meson):
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in quasi_da], [(val.mean - val.sdev) for val in quasi_da], color=color_list[0], alpha=0.5, label='Quasi')
    ax.fill_between(y_ls, [(val.mean + val.sdev) for val in lic_da], [(val.mean - val.sdev) for val in lic_da], color=color_list[1], alpha=0.7, label='Light-cone')
    #ax.plot(y_ls, [6*x*(1-x) for x in y_ls], color=color_list[2], label=r'$y=6x(1-x)$')
    ax.axvline(0.5, color='k', linestyle='--')
    #ax.axvline(0, color='k', linestyle='--')
    #ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA hybrid quasi v.s. light-cone '+plot_type+', Pz='+str(pz), **fs_p)
    ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig(meson+'/paper/quasi_v.s._light-cone, Pz='+str(pz)+'.pdf', transparent=True)
    plt.show()

def continuous_limit_pz_mix(meson, mom_ls):
    lam_dic = {}
    cont_dic = gv.BufferDict()

    for mom in mom_ls:
        hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_re_ls') # shape = (N_conf, N_z)
        hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_im_ls')

        pz = int(mom) * 0.215 # pz makes no difference
        lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

        a_re_ro, a_im_ro = rotate(hyb_re_ls, hyb_im_ls, lam_ls, back=False)

        a_ro_avg = gv.dataset.avg_data(a_re_ro, bstrap=True)

        lam_dic[str(mom)] = lam_ls
        cont_dic[str(mom)] = a_ro_avg


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for mom in mom_ls:
        ax.fill_between(lam_dic[str(mom)], [(val.mean + val.sdev) for val in cont_dic[str(mom)]], [(val.mean - val.sdev) for val in cont_dic[str(mom)]], color=color_list[int(mom/2-3)], alpha=0.4, label='mom='+str(mom))
    ax.axhline(0, color='k', linestyle='--')
    ax.legend(loc='upper right')
    ax.set_xlabel(lambda_label, **fs_p)
    ax.set_ylim([-0.5, 1.25])
    ax.tick_params(direction='in', **ls_p)
    plt.show()

def lcda_mix_pz_plot(meson, x_ls):
    mom_n1_lic_da = gv.load(meson+'/mom=6/mom_6_lic_da')
    mom_n2_lic_da = gv.load(meson+'/mom=8/mom_8_lic_da')
    mom_n3_lic_da = gv.load(meson+'/mom=10/mom_10_lic_da')

    if meson == 'kaon':
        mom_n4_lic_da = gv.load(meson+'/mom=12/mom_12_lic_da')

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    if meson == 'kaon':
        ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n4_lic_da], [(val.mean - val.sdev) for val in mom_n4_lic_da], color=color_list[3], alpha=0.5, label=r'$a \to 0, Pz=2.58GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n3_lic_da], [(val.mean - val.sdev) for val in mom_n3_lic_da], color=color_list[2], alpha=0.4, label=r'$a \to 0, Pz=2.15GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n2_lic_da], [(val.mean - val.sdev) for val in mom_n2_lic_da], color=color_list[1], alpha=0.4, label=r'$a \to 0, Pz=1.72GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n1_lic_da], [(val.mean - val.sdev) for val in mom_n1_lic_da], color=color_list[0], alpha=0.4, label=r'$a \to 0, Pz=1.29GeV$')
    ax.axvline(1, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(0.5, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(x_label, **fs_p)
    ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig(meson+'/paper/lcda_Pz_mix.pdf', transparent=True)
    plt.show()
    return 

def lcda_large_pz_plot(meson, x_ls, mom_n_lic_da, large_mom_lic_da):
    delta_ls = []
    for idx in range(len(x_ls)):
        delta = abs(large_mom_lic_da[idx].mean - mom_n_lic_da[idx].mean) # system error
        delta_ls.append(delta)

    y1 = np.array([(val.mean + val.sdev) for val in large_mom_lic_da]) + np.array(delta_ls)
    y2 = np.array([(val.mean - val.sdev) for val in large_mom_lic_da]) - np.array(delta_ls)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    a2 = gv.gvar(0.25, 0.15) # sum rule
    a4 = gv.gvar(-0.015, 0.025)
    ax.fill_between(x_ls, [sum_rule(x, a2, a4).mean + sum_rule(x, a2, a4).sdev for x in x_ls], [sum_rule(x, a2, a4).mean - sum_rule(x, a2, a4).sdev for x in x_ls], color=color_list[1], label='Sum rule', alpha=0.4)

    ax.fill_between(x_ls, y1, y2, color=color_list[0], alpha=0.5, label='This work')
    ax.plot(x_ls, [6*x*(1-x) for x in x_ls], color=color_list[2], label='Asymptotic') # only plot between 0 and 1

    ax.axvline(0.5, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA light-cone Pz to infty', **fs_p)
    ax.set_xlabel(x_label, **fs_p)
    ax.set_ylim([-0.19, 1.55])
    ax.set_xlim([-0.5, 1.5])
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig(meson+'/paper/lcda_Pz_to_infty.pdf', transparent=True)
    plt.show()
    return