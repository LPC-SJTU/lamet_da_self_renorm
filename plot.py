# %%
from head import *

def paper_plot_extra_ro(n=0.17, mom='10'):#(n=1.15):
    meson = 'pion'
    hyb_re_ls = {}
    hyb_im_ls = {}

    hyb_re_ls['mom='+mom] = gv.load(meson+'/mom='+str(mom)+'/quasi_re_ls') # shape = (N_conf, N_z)
    hyb_im_ls['mom='+mom] = gv.load(meson+'/mom='+str(mom)+'/quasi_im_ls')

    lam_dic = {}
    hyb_ro_avg = {}

    pz = int(mom) * 0.215
    lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * int(mom) * gev_fm ) / gev_fm
    # lam_ls = gv.load(meson+'/mom='+str(mom)+'/lam_ls')


    hyb_re_ro, hyb_im_ro = rotate(hyb_re_ls['mom='+mom], hyb_im_ls['mom='+mom], lam_ls, back=False)

    lam_dic['mom='+mom] = lam_ls
    hyb_ro_avg['mom='+mom] = gv.dataset.avg_data(hyb_re_ro, bstrap=True)

    lam_ls = gv.load('pion/mom='+mom+'/lam_ls')
    hyb_complex = gv.load('pion/mom='+mom+'/hyb_complex')

    re_conf_z = []
    im_conf_z = []
    for n_conf in range(len(hyb_complex)):
        re_conf_z.append([])
        im_conf_z.append([])
        for idx in range(len(hyb_complex[0])):
            re_conf_z[n_conf].append(hyb_complex[n_conf][idx].real)
            im_conf_z[n_conf].append(hyb_complex[n_conf][idx].imag)


    print(np.shape(re_conf_z))
    print(np.shape(lam_ls))
    
    re_ro_conf_z, im_ro_conf_z = rotate(re_conf_z, im_conf_z, lam_ls, back=False)

    bft_ro = gv.dataset.avg_data(re_ro_conf_z, bstrap=True)

    x_ls = np.linspace(0, 1, 500)
    a2_lcda = np.array([(x**n)*((1-x)**n) for x in x_ls])
    integral = np.sum(a2_lcda) * (x_ls[1] - x_ls[0])
    a2_lcda = a2_lcda / integral # normalization

    lam_ls_a2 = np.linspace(6, 15, 50)
    a2_bft = []
    for lam in lam_ls_a2:
        val = sum_ft_inv(x_ls, a2_lcda, x_ls[1]-x_ls[0], lam) * np.exp(1j * lam / 2)
        a2_bft.append(val.real)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(lam_dic['mom='+mom], [val.mean for val in hyb_ro_avg['mom='+mom]], [val.sdev for val in hyb_ro_avg['mom='+mom]], color=color_list[int(mom)-10], label='Lattice data, Pz=2.15 GeV', fmt='o', **errorb)
    ini = 216
    fin = ini+10
    ax.errorbar(lam_ls[ini-2:fin], [val.mean for val in bft_ro[ini-2:fin]], yerr=[val.sdev for val in bft_ro[ini-2:fin]], color=color_list[1], label='Polynomial fit', fmt='D',  **errorb)
    ax.plot(lam_ls_a2, a2_bft, color=color_list[2], label = r'$cx^n (1-x)^n, n = $'+str(round(n,2)))
    ax.axvline(lam_ls[ini], ymin=0.1, ymax=0.45, color='g', linestyle='--', label=r'$\lambda_L$')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA hybrid data mix Pz', **fs_p)
    ax.legend(loc='upper right')
    ax.set_xlabel(lambda_label, **fs_p)
    ax.set_ylim([-0.5, 1.25])
    ax.tick_params(direction='in', **ls_p)
    plt.savefig('pion/paper/data_extrapolation_ro_pz=2.15.pdf', transparent=True)
    plt.show()

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
        ax.set_title(hyb_ro_re_label, **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_rotated_Pz='+str(pz)+'GeV.pdf', transparent=True)
    elif if_rotate == False:
        ax.set_title(hyb_re_label, **fs_p)
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
        ax.set_title(hyb_ro_im_label, **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_rotated_imag_part_Pz='+str(pz)+'GeV.pdf', transparent=True)
    elif if_rotate == False:
        ax.set_ylim([-1.1, 0.5])
        ax.set_title(hyb_im_label, **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_imag_part_Pz='+str(pz)+'GeV.pdf', transparent=True)
    plt.show()
    return 

def quasi_vs_lc_plot(x_ls, quasi_da, lic_da, pz, meson):
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in quasi_da], [(val.mean - val.sdev) for val in quasi_da], color=color_list[0], alpha=0.5, label='Quasi')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in lic_da], [(val.mean - val.sdev) for val in lic_da], color=color_list[1], alpha=0.7, label='Light-cone')
    #ax.plot(x_ls, [6*x*(1-x) for x in x_ls], color=color_list[2], label=r'$y=6x(1-x)$')
    ax.axvline(0.5, color='green', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    #ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA hybrid quasi v.s. light-cone '+plot_type+', Pz='+str(pz), **fs_p)
    #ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig(meson+'/paper/quasi_v.s._light-cone, Pz='+str(pz)+'.pdf', transparent=True)
    plt.show()

def quasi_lc_lc_fit_plot(x_ls, y_ls, y_ls_, quasi_da, lic_da, lic_da_, pz, meson):
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in quasi_da], [(val.mean - val.sdev) for val in quasi_da], color=color_list[0], alpha=0.5, label='Quasi')
    gv.dump(x_ls, 'temp/x_ls_mu_3')
    gv.dump(quasi_da, 'temp/quasi_mu_3')
    ax.fill_between(y_ls_, [(val.mean + val.sdev) for val in lic_da_], [(val.mean - val.sdev) for val in lic_da_], color=color_list[1], alpha=0.7, label='Light-cone')
    gv.dump(y_ls_, 'temp/y_ls_mu_3')
    gv.dump(lic_da_, 'temp/lc_mu_3')
    ax.fill_between(y_ls, [(val.mean + val.sdev) for val in lic_da], [(val.mean - val.sdev) for val in lic_da], color='red', alpha=0.5, label='Extrapolated light-cone')
    #ax.plot(y_ls, [6*x*(1-x) for x in y_ls], color=color_list[2], label=r'$y=6x(1-x)$')
    ax.axvline(0.5, color='green', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    #ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA hybrid quasi v.s. light-cone '+plot_type+', Pz='+str(pz), **fs_p)
    ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig(meson+'/paper/quasi_lc_lc_fit, Pz='+str(pz)+'.pdf', transparent=True)
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
    mom_n1_lic_da = gv.load(meson+'/mom=6/lc_mom_ls')
    mom_n2_lic_da = gv.load(meson+'/mom=8/lc_mom_ls')
    mom_n3_lic_da = gv.load(meson+'/mom=10/lc_mom_ls')

    mom_n1_lic_da = gv.dataset.avg_data(mom_n1_lic_da, bstrap=True)
    mom_n2_lic_da = gv.dataset.avg_data(mom_n2_lic_da, bstrap=True)
    mom_n3_lic_da = gv.dataset.avg_data(mom_n3_lic_da, bstrap=True)


    # mom_n4_lic_da = gv.load(meson+'/mom=12/mom_12_lic_da')

    # if meson == 'kaon':
    #     mom_n4_lic_da = gv.load(meson+'/mom=12/mom_12_lic_da')

    fig = plt.figure(figsize=fig_size_lc)
    ax = plt.axes(plt_axes)
    #if meson == 'kaon':
    #ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n4_lic_da], [(val.mean - val.sdev) for val in mom_n4_lic_da], color=color_list[3], alpha=0.5, label=r'$a \to 0, Pz=2.58GeV$')
    
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n1_lic_da], [(val.mean - val.sdev) for val in mom_n1_lic_da], color=color_list[0], alpha=0.4, label=r'$a \to 0, Pz=1.29GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n2_lic_da], [(val.mean - val.sdev) for val in mom_n2_lic_da], color=color_list[1], alpha=0.4, label=r'$a \to 0, Pz=1.72GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n3_lic_da], [(val.mean - val.sdev) for val in mom_n3_lic_da], color=color_list[2], alpha=0.4, label=r'$a \to 0, Pz=2.15GeV$')
    
    ax.fill_between(np.linspace(-0.5, 0.05, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.6)
    ax.fill_between(np.linspace(0.95, 1.5, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.6)
    #ax.axvline(1, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(0.5, color='green', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(x_label, **fs_p)
    ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.legend(loc='lower center')
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


###
    add_err_gv = [ gv.gvar( (y1[i]+y2[i])/2, (y1[i]-y2[i])/2 ) for i in range(len(y1)) ]

    gv.dump(add_err_gv, 'lc_with_sys')

    print('>>> large mom limit a2:')
    a2 = calc_an(x_ls, add_err_gv, 2)
    print(a2)

    print('>>> large mom limit a2:')
    a4 = calc_an(x_ls, add_err_gv, 4)
    print(a4)

    mellin_moment(x_ls, add_err_gv, 2)
    mellin_moment(x_ls, add_err_gv, 4)
###


    fig = plt.figure(figsize=fig_size_lc)
    ax = plt.axes(plt_axes)

    a1 = gv.gvar(-0.06, 0.03) # sum rule
    a2 = gv.gvar(0.25, 0.15) # same for pion and kaon
    a4 = gv.gvar(-0.015, 0.025)


    ax.fill_between(x_ls[200:302], [sum_rule(meson, x, a1, a2, a4).mean + sum_rule(meson, x, a1, a2, a4).sdev for x in x_ls][200:302], [sum_rule(meson, x, a1, a2, a4).mean - sum_rule(meson, x, a1, a2, a4).sdev for x in x_ls][200:302], color=color_list[1], label='Sum rule', alpha=0.4)

    if meson == 'pion':
        a2 = gv.gvar(0.101, 0.024)
        ope = [sum_rule(meson, x, 0, a2, 0) for x in x_ls]
    elif meson == 'kaon':
        a1 = gv.gvar(-0.0533, 0.0034)
        a2 = gv.gvar(0.090, 0.019)
        ope = [sum_rule(meson, x, a1, a2, 0) for x in x_ls]

    ax.fill_between(x_ls[200:302], [val.mean + val.sdev for val in ope][200:302], [val.mean - val.sdev for val in ope][200:302], color=color_list[2], label='OPE', alpha=0.6)

    if meson == 'pion':
        ax.plot(x_ls[200:302], DSE(x_ls)[200:302], color='blue', label='DSE', linestyle='dashed')

    ax.fill_between(x_ls, y1, y2, color=color_list[0], alpha=0.5)

    ax.plot(x_ls, (y1+y2)/2, color=color_list[0], label='This work', linewidth=2, linestyle='dotted')

    ax.plot(x_ls, [6*x*(1-x) for x in x_ls], color='red', linestyle='dashdot', label='Asymptotic') # only plot between 0 and 1

    ax.fill_between(np.linspace(-0.5, 0.05, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.6)
    ax.fill_between(np.linspace(0.95, 1.5, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.6)

    ## grey v band to cover fit region

    ax.axvline(0.5, color='green', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA light-cone Pz to infty', **fs_p)
    ax.set_xlabel(x_label, **fs_p)
    ax.set_ylim([-0.19, 1.7])
    ax.set_xlim([-0.25, 1.25])
    ax.legend(loc='lower center')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig(meson+'/paper/lcda_Pz_to_infty.pdf', transparent=True)
    plt.show()
    return

if __name__ == '__main__':
    paper_plot_extra_ro()
# %%
