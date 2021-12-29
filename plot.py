# %%
from head import *
from pdf_self_renorm import pdf_zR

def quasi_vs_lc_plot(x_ls, y_ls, quasi_da, lic_da, pz, meson):
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in quasi_da], [(val.mean - val.sdev) for val in quasi_da], color=color_list[0], alpha=0.5, label='Quasi')
    ax.fill_between(y_ls, [(val.mean + val.sdev) for val in lic_da], [(val.mean - val.sdev) for val in lic_da], color=color_list[1], alpha=0.7, label='Light-cone')
    #ax.plot(x_ls, [6*x*(1-x) for x in x_ls], color=color_list[2], label=r'$y=6x(1-x)$')
    ax.axvline(0.5, color='green', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    #ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA hybrid quasi v.s. light-cone '+plot_type+', Pz='+str(pz), **fs_p)
    #ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(x_label, **fs_p_l)
    ax.legend(loc='upper right', **fs_p)
    ax.tick_params(direction='in', **ls_p_l)
    ax.grid(linestyle=':')
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
    ax = plt.axes(plt_axes_small)
    #if meson == 'kaon':
    #ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n4_lic_da], [(val.mean - val.sdev) for val in mom_n4_lic_da], color=color_list[3], alpha=0.5, label=r'$a \to 0, Pz=2.58GeV$')
    
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n1_lic_da], [(val.mean - val.sdev) for val in mom_n1_lic_da], color=color_list[0], alpha=0.4, label=r'$a \to 0, Pz=1.29GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n2_lic_da], [(val.mean - val.sdev) for val in mom_n2_lic_da], color=color_list[1], alpha=0.4, label=r'$a \to 0, Pz=1.72GeV$')
    ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n3_lic_da], [(val.mean - val.sdev) for val in mom_n3_lic_da], color=color_list[2], alpha=0.4, label=r'$a \to 0, Pz=2.15GeV$')
    
    ax.fill_between(np.linspace(-0.5, 0.1, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.4)
    ax.fill_between(np.linspace(0.9, 1.5, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.4)
    #ax.axvline(1, color='k', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(0.5, color='green', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(x_label, **fs_p_l)
    ax.set_ylim([-0.19, 1.5])
    ax.set_xlim([-0.5, 1.5])
    ax.legend(loc='lower center')
    ax.tick_params(direction='in', **ls_p_l)
    plt.savefig(meson+'/paper/lcda_Pz_mix.pdf', transparent=True)
    plt.show()
    return 

def lcda_large_pz_plot(meson, x_ls, mom_n_lic_da, large_mom_lic_da):
    ### replace all [:] with [201:302] for plot with tails ###
    delta_ls = []
    for idx in range(len(x_ls)):
        delta = abs(large_mom_lic_da[idx].mean - mom_n_lic_da[idx].mean) # system error
        delta_ls.append(delta)


    y1 = np.array([(val.mean + val.sdev) for val in large_mom_lic_da]) + np.array(delta_ls)
    y2 = np.array([(val.mean - val.sdev) for val in large_mom_lic_da]) - np.array(delta_ls)

    x_ls = np.hstack((x_ls, np.array([1])))
    y1 = np.hstack((y1, np.array([0])))
    y2 = np.hstack((y2, np.array([0])))

    ###
    add_err_gv = [ gv.gvar( (y1[i]+y2[i])/2, (y1[i]-y2[i])/2 ) for i in range(len(y1)) ]

    # print('>>> large mom limit a2:')
    # a2 = calc_an(x_ls, add_err_gv, 2)
    # print(a2)

    # print('>>> large mom limit a4:')
    # a4 = calc_an(x_ls, add_err_gv, 4)
    # print(a4)

    # mellin_moment(x_ls, add_err_gv, 2)
    # mellin_moment(x_ls, add_err_gv, 4)

    ###


    fig = plt.figure(figsize=fig_size_lc)
    ax = plt.axes(plt_axes_small)

    a1 = gv.gvar(-0.06, 0.03) # sum rule
    a2 = gv.gvar(0.25, 0.15) # same for pion and kaon
    a4 = gv.gvar(-0.015, 0.025)


    ax.fill_between(x_ls[201:302], [sum_rule(meson, x, a1, a2, a4).mean + sum_rule(meson, x, a1, a2, a4).sdev for x in x_ls][201:302], [sum_rule(meson, x, a1, a2, a4).mean - sum_rule(meson, x, a1, a2, a4).sdev for x in x_ls][201:302], color=color_list[1], label='Sum rule', alpha=0.4)

    if meson == 'pion':
        a2 = gv.gvar(0.101, 0.024)
        ope = [sum_rule(meson, x, 0, a2, 0) for x in x_ls]
    elif meson == 'kaon':
        a1 = gv.gvar(-0.0533, 0.0034)
        a2 = gv.gvar(0.090, 0.019)
        ope = [sum_rule(meson, x, a1, a2, 0) for x in x_ls]

    ax.fill_between(x_ls[201:302], [val.mean + val.sdev for val in ope][201:302], [val.mean - val.sdev for val in ope][201:302], color=color_list[2], label='OPE', alpha=0.6)

    if meson == 'pion':
        ax.plot(x_ls[201:302], DSE(x_ls)[201:302], color='blue', label='DSE', linestyle='dashed')

    elif meson == 'kaon':
        dse_x, dse_y = DSE_kaon()
        ax.plot(dse_x, dse_y, color='blue', label='DSE', linestyle='dashed')

    

    ax.fill_between(x_ls, y1, y2, color=color_list[0], alpha=0.5)

    # gv.dump(x_ls, 'temp/k_fit_x')
    # gv.dump(y1, 'temp/k_fit_y1')
    # gv.dump(y2, 'temp/k_fit_y2')

    fit_x_ls = gv.load('temp/fit_x')
    fit_y1 = gv.load('temp/fit_y1')
    fit_y2 = gv.load('temp/fit_y2')

    # ax.fill_between(fit_x_ls, fit_y1, fit_y2, color='green', alpha=0.3) # endpoints fit

    ax.plot(x_ls, (y1+y2)/2, color=color_list[0], label='This work', linewidth=2, linestyle='dotted')

    ax.plot(x_ls[201:302], [6*x*(1-x) for x in x_ls][201:302], color='red', linestyle='dashdot', label='Asymptotic') # only plot between 0 and 1

    ax.fill_between(np.linspace(-0.5, 0.1, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.2)
    ax.fill_between(np.linspace(0.9, 1.5, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.2)

    ## grey v band to cover fit region

    ax.axvline(0.5, color='green', linestyle='--')
    ax.axvline(0, color='k', linestyle='--')
    ax.axvline(1, color='k', linestyle='--')
    ax.axhline(0, color='k', linestyle='--')
    #ax.set_title('DA light-cone Pz to infty', **fs_p)
    ax.set_xlabel(x_label, **fs_p_l)
    ax.set_ylim([-0.19, 1.7])
    ax.set_xlim([-0.25, 1.25])
    ax.legend(loc='lower center')
    ax.tick_params(direction='in', **ls_p_l)
    plt.savefig(meson+'/paper/lcda_Pz_to_infty.pdf', transparent=True)
    plt.show()
    return

def extrapolation_check(title, lam_ls, lc_ls, po_1, po_2):
    lc_avg = gv.dataset.avg_data(lc_ls, bstrap=True)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    ax.errorbar(lam_ls, [val.mean for val in lc_avg], [val.sdev for val in lc_avg], color=blue, marker='D', label='lc_coor', **errorb)

    ax.axvline(po_1, ymin=0.1, ymax=0.45, color='g', linestyle='--', label='fit include')
    ax.axvline(po_2, ymin=0.1, ymax=0.45, color='r', linestyle='dashdot', label='extra start')

    ax.axhline(0, color='k', linestyle='--', lw=0.8)
    ax.set_xlabel(lambda_label, **fs_p)
    ax.set_ylim([-0.99, 1.59])
    ax.set_title(title)
    ax.tick_params(direction='in', **ls_p)
    ax.legend(loc='upper right')
    plt.show()

    return

def fig_1():
    meson = 'pion'
    mom = 10

    zR_dic, m_pdf_dic = pdf_zR()

    data_ls = gv.load(meson+'/mom='+str(mom)+'/da_an_ls')
    quasi_re_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_re_ls')
    quasi_im_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_im_ls')

    key_ls = ['a=0.0574', 'a=0.0882', 'a=0.1213']

    renorm_data_re_ls = []
    renorm_data_im_ls = []

    def renorm(da_conf, key):
        da_re_ls = []
        da_im_ls = []
        for n_conf in range(len(da_conf_ls)):
            temp_re_ls = []
            temp_im_ls = []
            for idx in range(len(da_conf_ls[0])):
                temp_val = (da_conf_ls[n_conf][idx] / zR_dic[key][idx].mean) # divide by zR
                temp_re_ls.append(temp_val.real)
                temp_im_ls.append(temp_val.imag)

            da_re_ls.append(temp_re_ls)
            da_im_ls.append(temp_im_ls)

        return da_re_ls, da_im_ls

    for id in range(len(data_ls)):
        da_conf_ls = data_ls[id]
        key = key_ls[id]
        re, im = renorm(da_conf_ls, key)
        renorm_data_re_ls.append(re)
        renorm_data_im_ls.append(im)

    len_dic = {}
    len_dic['a=0.0574'] = 0.0574*96
    len_dic['a=0.0882'] = 0.0882*64
    len_dic['a=0.1213'] = 0.1213*48

    label_ls= [r'$a: 0.06 fm$', r'$a: 0.09 fm$', r'$a: 0.12 fm$']


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for id in range(len(renorm_data_re_ls)): 
        key = key_ls[id]
        lam_ls = z_ls_extend * ( 2*np.pi / (len_dic[key]) * mom * gev_fm ) / gev_fm
        re_ro, im_ro = rotate(renorm_data_re_ls[id], renorm_data_im_ls[id], lam_ls, back=False)
        y_ls = gv.dataset.avg_data(re_ro, bstrap=True)
        z_ls = z_ls_extend
        ax.errorbar(z_ls, [val.mean for val in y_ls], [val.sdev for val in y_ls], color=color_list[id], marker=fmt_ls[id], label=label_ls[id], **errorb)

    lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm
    q_re_ro, q_im_ro = rotate(quasi_re_ls, quasi_im_ls, lam_ls, back=False)

    q_avg = gv.dataset.avg_data(q_re_ro, bstrap=True)

    for id in range(len(q_avg)):
        q_avg[id] = q_avg[id] * ZMS_da(z_ls_da[id])

    ax.fill_between(z_ls_da, [val.mean+val.sdev for val in q_avg], [val.mean-val.sdev for val in q_avg], color=color_list[3], alpha=0.2, label=r'$a \to 0$')

    ax.plot(z_ls_da, ZMS_da(z_ls_da), 'r-', label='ZMS-bar')

    ax.axhline(0, color='k', linestyle='--', lw=0.8)
    ax.set_xlim([0, 1.65])
    ax.set_xlabel(z_label, **fs_p)
    ax.set_ylim([-0.49, 1.99])
    ax.set_title(hyb_ro_re_label, **fs_p)
    ax.tick_params(direction='in', **ls_p)
    ax.legend(loc='upper right')
    plt.savefig(meson+'/paper/fig1.pdf', transparent=True)
    plt.show()

    return 

def fig_1_1():
    meson = 'pion'
    mom_ls = [6, 8, 10]

    zR_dic, m_pdf_dic = pdf_zR()

    label_dic = {}
    label_dic['6'] = r'$a \to 0, Pz=1.29GeV$'
    label_dic['8'] = r'$a \to 0, Pz=1.72GeV$'
    label_dic['10'] = r'$a \to 0, Pz=2.15GeV$'

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)

    for mom in mom_ls:
        quasi_re_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_re_ls')
        quasi_im_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_im_ls')
    
        lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm
        q_re_ro, q_im_ro = rotate(quasi_re_ls, quasi_im_ls, lam_ls, back=False)

        q_avg = gv.dataset.avg_data(q_re_ro, bstrap=True)

        for id in range(len(q_avg)):
            q_avg[id] = q_avg[id] * ZMS_da(z_ls_da[id])

        ax.fill_between(z_ls_da, [val.mean+val.sdev for val in q_avg], [val.mean-val.sdev for val in q_avg], color=color_list[int(mom/2-3)], alpha=0.4, label=label_dic[str(mom)])

    ax.plot(z_ls_da, ZMS_da(z_ls_da), 'r-', label='ZMS-bar')

    ax.axhline(0, color='k', linestyle='--', lw=0.8)
    ax.set_xlim([0, 1.4])
    ax.set_xlabel(z_label, **fs_p_l)
    ax.set_ylim([-0.49, 1.99])
    ax.set_title(hyb_ro_re_label, **fs_p_l)
    ax.tick_params(direction='in', **ls_p_l)
    ax.grid(linestyle=':')
    ax.legend(loc=(0.55, 0.48), **fs_p)
    plt.savefig(meson+'/paper/fig1-1.pdf', transparent=True)
    plt.show()

def fig_1_2(meson, mom, if_rotate=False):
    lam_ls = gv.load(meson+'/mom='+str(mom)+'/lam_ls')
    re_lam_ls = gv.load(meson+'/mom='+str(mom)+'/re_lam_ls')
    im_lam_ls = gv.load(meson+'/mom='+str(mom)+'/im_lam_ls')

    quasi_re_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_re_ls')
    quasi_im_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_im_ls')


    pz = int(mom) * 0.215
    lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

    if if_rotate == True:
        a06_re_ro, a06_im_ro = rotate(re_lam_ls[0], im_lam_ls[0], lam_ls, back=False)
        a09_re_ro, a09_im_ro = rotate(re_lam_ls[1], im_lam_ls[1], lam_ls, back=False)
        a12_re_ro, a12_im_ro = rotate(re_lam_ls[2], im_lam_ls[2], lam_ls, back=False)

        a_re_ro, a_im_ro = rotate(quasi_re_ls, quasi_im_ls, lam_ls, back=False)

    elif if_rotate == False:
        a06_re_ro, a06_im_ro = re_lam_ls[0], im_lam_ls[0]
        a09_re_ro, a09_im_ro = re_lam_ls[1], im_lam_ls[1]
        a12_re_ro, a12_im_ro = re_lam_ls[2], im_lam_ls[2]

        a_re_ro, a_im_ro = quasi_re_ls, quasi_im_ls

    ## real part
    a06_ro_avg = gv.dataset.avg_data(a06_re_ro, bstrap=True)
    a09_ro_avg = gv.dataset.avg_data(a09_re_ro, bstrap=True)
    a12_ro_avg = gv.dataset.avg_data(a12_re_ro, bstrap=True)

    a_ro_avg = gv.dataset.avg_data(a_re_ro, bstrap=True)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)
    ax.errorbar(lam_ls, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a09_ro_avg], [val.sdev for val in a09_ro_avg], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)

    ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_ro_avg], [(val.mean - val.sdev) for val in a_ro_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(lambda_label, **fs_p_l)
    ax.set_ylim([-0.49, 1.249])
    ax.tick_params(direction='in', **ls_p_l)
    ax.grid(linestyle=':')
    if if_rotate == True:
        ax.set_title(hyb_ro_re_label, **fs_p_l)
        ax.legend(loc='upper right', **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_rotated_Pz='+str(pz)+'GeV.pdf', transparent=True)
    elif if_rotate == False:
        ax.set_title(hyb_re_label, **fs_p_l)
        ax.legend(loc='upper right', **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_Pz='+str(pz)+'GeV.pdf', transparent=True)
    plt.show()

    ## imag part
    a06_ro_avg = gv.dataset.avg_data(a06_im_ro, bstrap=True)
    a09_ro_avg = gv.dataset.avg_data(a09_im_ro, bstrap=True)
    a12_ro_avg = gv.dataset.avg_data(a12_im_ro, bstrap=True)

    a_ro_avg = gv.dataset.avg_data(a_im_ro, bstrap=True)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)
    ax.errorbar(lam_ls, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a09_ro_avg], [val.sdev for val in a09_ro_avg], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)

    ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_ro_avg], [(val.mean - val.sdev) for val in a_ro_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(lambda_label, **fs_p_l)
    ax.tick_params(direction='in', **ls_p_l)
    ax.grid(linestyle=':')
    if if_rotate == True:
        ax.set_ylim([-0.29, 0.29])
        ax.set_title(hyb_ro_im_label, **fs_p_l)
        ax.legend(loc='lower right', ncol=2, **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_rotated_imag_part_Pz='+str(pz)+'GeV.pdf', transparent=True)
    elif if_rotate == False:
        ax.set_ylim([-1.09, 0.49])
        ax.set_title(hyb_im_label, **fs_p_l)
        ax.legend(loc='lower right', **fs_p)
        plt.savefig(meson+'/paper/discrete_effect_imag_part_Pz='+str(pz)+'GeV.pdf', transparent=True)
    plt.show()
    return 

def fig_2(n=0.18, mom='10', lambdaL=8.209736899624895, fit_start=5.473157933083264):
    meson = 'pion'

    # before extrapolation
    lam_ls = gv.load(meson+'/mom='+str(mom)+'/lam_ls')
    lc_re_ls = gv.load(meson+'/mom='+str(mom)+'/lc_re_ls')
    lc_im_ls = gv.load(meson+'/mom='+str(mom)+'/lc_im_ls')

    lc_re_ro, lc_im_ro = rotate(lc_re_ls, lc_im_ls, lam_ls, back=False)

    lc_avg = gv.dataset.avg_data(lc_re_ro, bstrap=True)

    # after extrapolation
    lam_ls_ex = gv.load(meson+'/mom='+str(mom)+'/lam_ls_ex')
    lc_ext_ls = gv.load(meson+'/mom='+str(mom)+'/lc_ext_ls')

    st = int(fit_start / lam_ls[0]) - 1 # this lambda point is included in the fit
    leng = 18 # how many fit points are plotted


    ## cx^n * (1-x)^n ##
    x_ls = np.linspace(0, 1, 500)
    a2_lcda = np.array([(x**n)*((1-x)**n) for x in x_ls])
    integral = np.sum(a2_lcda) * (x_ls[1] - x_ls[0])
    a2_lcda = a2_lcda / integral # normalization

    lam_ls_a2 = np.linspace(6, 18, 80)
    a2_bft = []
    for lam in lam_ls_a2:
        val = sum_ft_inv(x_ls, a2_lcda, x_ls[1]-x_ls[0], lam) * np.exp(1j * lam / 2)
        a2_bft.append(val.real)

    ## fit points ##
    lc_re_ex_ls = []
    lc_im_ex_ls = []
    for n_conf in range(len(lc_ext_ls)):
        lc_re_ex_ls.append([])
        lc_im_ex_ls.append([])
        for idl in range(1001, len(lam_ls_ex)):
            lc_re_ex_ls[n_conf].append(lc_ext_ls[n_conf][idl].real)
            lc_im_ex_ls[n_conf].append(lc_ext_ls[n_conf][idl].imag)

    lc_re_ex_ro, lc_im_ex_ro = rotate(lc_re_ex_ls, lc_im_ex_ls, lam_ls_ex[1001:], back=False)

    lc_ex_avg = gv.dataset.avg_data(lc_re_ex_ro, bstrap=True)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)

    ax.plot(lam_ls_a2, a2_bft, color=color_list[2], label = r'$cx^n (1-x)^n, n = $'+str(round(n,2)))

    ax.errorbar(lam_ls, [val.mean for val in lc_avg], [val.sdev for val in lc_avg], color='red', marker='o', label='Lattice data, Pz=2.15 GeV', **errorb)

    ax.errorbar(lam_ls_ex[1001+st:1001+st+leng], [val.mean for val in lc_ex_avg[st:st+leng]], [val.sdev for val in lc_ex_avg[st:st+leng]], color=color_list[1], marker='D', label='Polynomial fit', **errorb)

    ax.axvline(lambdaL, ymin=0.03, ymax=0.38, color='green', linestyle='--', label=r'$\lambda_L$')

    ax.axhline(0, color='k', linestyle='--', lw=0.8)
    ax.set_xlabel(lambda_label, **fs_p_l)
    ax.set_ylim([-0.29, 1.19])
    #ax.set_title(title)
    ax.tick_params(direction='in', **ls_p_l)
    ax.legend(loc='upper right', **fs_p)
    ax.grid(linestyle=':')
    plt.savefig('pion/paper/fig2.pdf', transparent=True)
    plt.show()

    return 


if __name__ == '__main__':
    #fig_1()
    #fig_1_1()
    #fig_1_2(meson='pion', mom=10, if_rotate=True)

    fig_2()

    
# %%
