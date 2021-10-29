# %%
from head import *
from pdf_self_renorm import pdf_zR

def es_calc(mom, an, L):
    pz = 2*np.pi / (an * L) * mom * gev_fm

    m0 = 0.14
    m1 = gv.gvar(1.42, 0.027) # GeV
    m2 = gv.gvar(2.05, 0.035)
    m3 = gv.gvar(2.28, 0.098)

    e0 = np.sqrt(m0**2 + pz**2)
    e1 = np.sqrt(m1**2 + pz**2)
    e2 = np.sqrt(m2**2 + pz**2)
    e3 = np.sqrt(m3**2 + pz**2)

    dE1 = (e1 - e0) * an / gev_fm
    dE2 = (e2 - e0) * an / gev_fm
    dE3 = (e3 - e0) * an / gev_fm
    
    return dE1, dE2, dE3

def gs_fit_joint(t_ls, Cz_0_re, Cz_0_im, C0_re, mom, an, L):
    dE1, dE2, dE3 = es_calc(mom, an, L)

    re_m = np.average([val.mean for val in Cz_0_re])
    re_s = np.std([val.mean for val in Cz_0_re])
    im_m = np.average([val.mean for val in Cz_0_im])
    im_s = np.std([val.mean for val in Cz_0_im])

    priors = gv.BufferDict()
    priors['gs_re'] = gv.gvar(re_m, 3*re_s)
    priors['gs_im'] = gv.gvar(im_m, 3*im_s)
    priors['a1'] = gv.gvar(0, 10)
    priors['a2'] = gv.gvar(0, 10)
    priors['b1_re'] = gv.gvar(1, 5)
    priors['b2_re'] = gv.gvar(1, 5)
    priors['b1_im'] = gv.gvar(1, 5)
    priors['b2_im'] = gv.gvar(1, 5)
    priors['c'] = gv.gvar(0, 10)
    priors['dE1'] = dE1
    priors['dE2'] = dE2
    priors['e0'] = gv.gvar(0.5, 10)

    def fcn(x, p):
        val = {}
        val['z_0_re'] = p['gs_re'] * ( 1 + p['b1_re'] * np.exp(- p['dE1'] * x['z_0_re']) + p['b2_re'] * np.exp(- p['dE2'] * x['z_0_re']) ) / ( 1 + p['a1'] * np.exp(- p['dE1'] * x['z_0_re']) + p['a2'] * np.exp(- p['dE2'] * x['z_0_re']) )

        val['z_0_im'] = p['gs_im'] * ( 1 + p['b1_im'] * np.exp(- p['dE1'] * x['z_0_im']) + p['b2_im'] * np.exp(- p['dE2'] * x['z_0_im']) ) / ( 1 + p['a1'] * np.exp(- p['dE1'] * x['z_0_im']) + p['a2'] * np.exp(- p['dE2'] * x['z_0_im']) )

        val['0'] = p['c'] * np.exp(- p['e0'] * x['0']) * ( 1 + p['a1'] * np.exp(- p['dE1'] * x['0']) + p['a2'] * np.exp(- p['dE2'] * x['0']) )

        return val

    x_dic = {}
    x_dic['0'], x_dic['z_0_re'], x_dic['z_0_im'] = t_ls, t_ls, t_ls

    y_dic = {}
    y_dic['0'], y_dic['z_0_re'], y_dic['z_0_im'] = C0_re, Cz_0_re, Cz_0_im

    fit_result = lsf.nonlinear_fit(data=(x_dic, y_dic), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    if fit_result.chi2 / fit_result.dof > 5:
        print(fit_result.format(maxline=True))
        print(y_dic)

    return fit_result

def gs_fit_chained(t_ls, Cz_re, Cz_im, C0_re):

    return 

def t_plot(a_str='a06', file_path='DA_new.hdf5'):
    myfile = h5.File(file_path,'r')
    meson = 'pion'

    N_re = 600

    if a_str == 'a06':
        an = 0.0574
        ll = 96
        Nz = 49 # num of z
        Nt = 23
    elif a_str == 'a09':
        an = 0.0882
        ll = 64
        Nz = 33
        Nt = 15
    elif a_str == 'a12':
        an = 0.1213
        ll = 48
        Nz = 25
        Nt = 15
    else:
        print('a input error')

    an_pi = myfile[a_str+'m130_'+meson][int(mom/2)] # conf, z, t, real/imag

    an_pi_bs = bootstrap(an_pi, N_re) # resampling

    for n_conf in range(N_re): # normalization
        for idt in range(Nt):
            for idz in range(1, Nz):
                an_pi_bs[n_conf][idz][idt][0] = an_pi_bs[n_conf][idz][idt][0] / an_pi_bs[n_conf][0][idt][0]

                an_pi_bs[n_conf][idz][idt][1] = an_pi_bs[n_conf][idz][idt][1] / an_pi_bs[n_conf][0][idt][0]

    an_pi_bs_avg = gv.dataset.avg_data(an_pi_bs, bstrap=True)

    plt.figure()
    plt.errorbar(np.arange(1, Nt+1), [val.mean for val in an_pi_bs_avg[0,:,0]], yerr=[val.sdev for val in an_pi_bs_avg[0,:,0]], **errorb)
    plt.title('z=0 real')
    plt.show()

    m_eff = [np.log(an_pi_bs_avg[0,t,0] / an_pi_bs_avg[0,t+1,0]) for t in range(Nt - 1)]

    plt.figure()
    plt.errorbar(np.arange(1, Nt)[:17], [val.mean for val in m_eff][:17], yerr=[val.sdev for val in m_eff][:17], **errorb)
    plt.title('z=0 effective mass')
    #plt.ylim([0, 1])
    plt.show()

    ratio_t_z_re = []
    ratio_t_z_im = []
    for n_conf in range(N_re):
        ratio_t_z_re.append([])
        ratio_t_z_im.append([])
        for idt in range(Nt):
            z_in = an*np.arange(Nz)
            y_re_in = []
            y_im_in = []
            for idz in range(Nz):
                y_re_in.append( an_pi_bs[n_conf][idz][idt][0] )
                y_im_in.append( an_pi_bs[n_conf][idz][idt][1] )
            y_re_out = interp_1d(z_in, y_re_in, z_ls_extend, method='cubic')
            y_im_out = interp_1d(z_in, y_im_in, z_ls_extend, method='cubic')
            ratio_t_z_re[n_conf].append(y_re_out)
            ratio_t_z_im[n_conf].append(y_im_out)

    ratio_avg_t_z_re = gv.dataset.avg_data(ratio_t_z_re, bstrap=True)
    ratio_avg_t_z_im = gv.dataset.avg_data(ratio_t_z_im, bstrap=True)

    for idt in range(Nt):
        for idz in range(len(z_ls_extend)):
            temp = ratio_avg_t_z_re[idt][idz] / zR_dic['a=0.0574'][idz].mean / ZMS_da(z_ls_extend[idz])
            ratio_avg_t_z_re[idt][idz] = temp

            temp = ratio_avg_t_z_im[idt][idz] / zR_dic['a=0.0574'][idz].mean / ZMS_da(z_ls_extend[idz])
            ratio_avg_t_z_im[idt][idz] = temp

    for z in [3]:
        t_max = 16

        constrain_re = [-0.188, 0.018]
        no_constrain_re = [0.667, 0.016]
        constrain_im = []
        no_constrain_im = []
        re_789 = [0.4406, 0.0038]
        re_8910 = [0.4330, 0.0090]
        re_101112 = [0.429, 0.022]

        print('z='+str(z))

        real = [ratio_avg_t_z_re[t,z] for t in range(Nt)]

        imag = [ratio_avg_t_z_im[t,z] for t in range(Nt)]

        t_nc = np.arange(1, 18)

        ft_nc = np.array([0.46209287254343306, 0.4567188618226018, 0.45237384776582223, 0.44890260251869213, 0.44615583887765364, 0.44399877145926353, 0.44231486652327934, 0.44100643999882766, 0.43999344080869035, 0.43921136060981086, 0.43860886599824145, 0.4381454932709064, 0.43778957490400916, 0.4375164620449299, 0.4373070482460298, 0.43714656989303285, 0.4370236465088145])
        ft_nc_err = np.array([0.011894921297856766, 0.006044369781270242, 0.0025328663592327086, 0.0008520313034093377, 0.0007509884183954541, 0.0009547550568640708, 0.001258217863248, 0.001829024213467061, 0.002578529162235665, 0.003384504145792482, 0.004171893527437973, 0.00490108753382541, 0.005553836535301537, 0.006124539908123299, 0.006614921914114232, 0.007030708392814789, 0.007379546316894823])

        ft_pdg_e1 = np.array([0.45744261459640023, 0.45431631568057085, 0.4514317836847471, 0.44878698189061167, 0.4463759119712952, 0.4441894303297636, 0.44221603927076314, 0.4404426137759917, 0.43885503939227416, 0.4374387494303074, 0.4361791596030105, 0.43506200529269834, 0.4340735910946181, 0.4332009646101208, 0.4324320271851501, 0.43175559390782603, 0.4311614141135824])
        ft_pdg_e1_err = np.array([0.0043206530208204945, 0.002673971774074046, 0.0015166344813426576, 0.0008093066005710153, 0.0006365082756993438, 0.0008807130775536214, 0.001253076440928197, 0.0016921336534826666, 0.002193462333909812, 0.00274769113468605, 0.0033391450911918986, 0.0039507994600952836, 0.0045672672514321284, 0.005175921507611261, 0.005767044070239807, 0.0063335615341497415, 0.006870642683168408])

        ft_p0_e1 = np.array([0.4520196617340565, 0.4507652599936559, 0.4494315170062697, 0.44803167884183975, 0.4465823177891829, 0.4451026612706968, 0.4436136296825387, 0.4421366778770889, 0.44069257298682096, 0.4393002540158016, 0.43797590231331507, 0.43673231153039843, 0.43557859243184166, 0.43452019568797606, 0.43355919570427837, 0.43269475644011013, 0.4319236960010304])
        ft_p0_e1_err = np.array([0.0008672468132385717, 0.0006017665295722042, 0.00044296595544091185, 0.0004424661603604554, 0.0005735368363927765, 0.00078113980938197, 0.0010544511679248848, 0.0014011885467779949, 0.0018263682120210531, 0.0023269358325586718, 0.002892841958618943, 0.0035095996741722067, 0.004160647214715377, 0.004829185661018739, 0.005499505906584761, 0.006157864192320706, 0.006792968382631826])

        ft_pdg_e2 = np.array([0.44434895251840656, 0.44948275962952494, 0.45048733570234856, 0.4491501071626159, 0.44681256862645774, 0.4442175868251225, 0.44171041027210645, 0.43942576280119816, 0.4373988764063555, 0.4356220153196186, 0.4340712440706468, 0.4327184788952716, 0.431536576798186, 0.43050122202094177, 0.4295913866539997, 0.4287891932189064, 0.4280795588621008])
        ft_pdg_e2_err = np.array([0.0027649680604949216, 0.0007346711405030844, 0.0005753587674186594, 0.0005876747459829563, 0.0006295304971756503, 0.0008509319615954016, 0.0012289639138392352, 0.0016926694746730044, 0.002189470709277645, 0.002686745486028806, 0.0031657526734087506, 0.0036167279932486744, 0.004035443115031543, 0.00442093919696772, 0.004774097022078928, 0.005096761654314233, 0.005391219385522737])

        fig = plt.figure(figsize=fig_size)
        ax = plt.axes(plt_axes)
        ax.errorbar(np.arange(1, Nt+1)[:t_max], [val.mean for val in real][:t_max], yerr=[val.sdev for val in real][:t_max], color=color_list[3], label='Data', **errorb)
        #ax.fill_between(np.arange(1, Nt+1)[:t_max], np.ones(t_max)*(constrain_re[0]+constrain_re[1]), np.ones(t_max)*(constrain_re[0]-constrain_re[1]), color=color_list[4], label='constrain', alpha=0.4)
        ax.fill_between(t_nc, ft_nc+ft_nc_err, ft_nc-ft_nc_err, color=color_list[5], label='Two state fit 1', alpha=0.4)
        ax.fill_between(t_nc, ft_p0_e1+ft_p0_e1_err, ft_p0_e1-ft_p0_e1_err, color=color_list[7], label='Two state fit 2', alpha=0.4)
        # ax.fill_between(t_nc, ft_pdg_e1+ft_pdg_e1_err, ft_pdg_e1-ft_pdg_e1_err, color=color_list[9], label='E1 from pdg', alpha=0.4)
        # ax.fill_between(t_nc, ft_pdg_e2+ft_pdg_e2_err, ft_pdg_e2-ft_pdg_e2_err, color=color_list[2], label='E1E2 from pdg', alpha=0.4)

        ax.fill_between(np.array([7,8,9]), np.ones(3)*(re_789[0]+re_789[1]), np.ones(3)*(re_789[0]-re_789[1]), color=color_list[10], label='One state fit', alpha=1)
        # ax.fill_between(np.array([8,9,10]), np.ones(3)*(re_8910[0]+re_8910[1]), np.ones(3)*(re_8910[0]-re_8910[1]), color='red', label='fix t=8,9,10', alpha=0.4)
        # ax.axvline(4, linestyle='--', color='orange')
        # ax.axvline(6, linestyle='--', color='orange')
        ax.set_title(r'$\phi_2(z=4, t) / \phi_2(z=0, t)$', **fs_p)
        ax.set_ylim([0.3, 0.6])
        ax.set_xlabel('t', **fs_p)
        ax.tick_params(direction='in', **ls_p)
        ax.legend(loc='upper right')
        plt.savefig(meson+'/paper/gs_fit.pdf', transparent=True)
        plt.show()

        #print(np.average([val.mean for val in real[2:5]]))

        # fig = plt.figure(figsize=fig_size)
        # ax = plt.axes(plt_axes)
        # ax.errorbar(np.arange(1, Nt+1)[:t_max], [val.mean for val in imag][:t_max], yerr=[val.sdev for val in imag][:t_max], color=color_list[3], label='data', **errorb)
        # ax.fill_between(np.arange(1, Nt+1)[:t_max], np.ones(t_max)*(constrain_im[0]+constrain_im[1]), np.ones(t_max)*(constrain_im[0]-constrain_im[1]), color=color_list[4], label='constrain', alpha=0.4)
        # ax.fill_between(np.arange(1, Nt+1)[:t_max], np.ones(t_max)*(no_constrain_im[0]+no_constrain_im[1]), np.ones(t_max)*(no_constrain_im[0]-no_constrain_im[1]), color=color_list[5], label='no constrain', alpha=0.4)
        # # ax.axvline(4, linestyle='--', color='orange')
        # # ax.axvline(6, linestyle='--', color='orange')
        # ax.set_title('ratio imag, z='+str(z), **fs_p)
        # #ax.set_ylim([-0.55, -0.4])
        # ax.set_xlabel('t', **fs_p)
        # ax.tick_params(direction='in', **ls_p)
        # ax.legend(loc='upper right')
        # plt.show()

def stability_analysis_single_z(file_path='DA_new.hdf5'):
    myfile = h5.File(file_path,'r')
    meson = 'pion'
    a_str='a12'

    if a_str == 'a06':
        an = 0.0574
        ll = 96
        N_re = 600 # same times for eliminating discretization effect
        Nz = 49 # num of z
        Nt = 23
        t_ls = np.arange(5, 13)#(3, 11)
    elif a_str == 'a09':
        an = 0.0882
        ll = 64
        N_re = 600
        Nz = 33
        Nt = 15
        t_ls = np.arange(5, 12)#(3, 11)
    elif a_str == 'a12':
        an = 0.1213
        ll = 48
        N_re = 600
        Nz = 25
        Nt = 15
        t_ls = np.arange(3, 7)
    else:
        print('a input error')

    pz = 2*np.pi / (an * ll) * mom * gev_fm

    an_pi = myfile[a_str+'m130_'+meson][int(mom/2)] # conf, z, t, real/imag

    an_pi_bs = bootstrap(an_pi, N_re) # resampling

    for n_conf in range(N_re): # normalize
        for idt in range(Nt):
            for idz in range(1, Nz):
                an_pi_bs[n_conf][idz][idt][0] = an_pi_bs[n_conf][idz][idt][0] / an_pi_bs[n_conf][0][idt][0]

                an_pi_bs[n_conf][idz][idt][1] = an_pi_bs[n_conf][idz][idt][1] / an_pi_bs[n_conf][0][idt][0]

    an_pi_bs_avg = gv.dataset.avg_data(an_pi_bs, bstrap=True)
   
    idz = 15

    chi2_ls = []
    gs_re_ls = []
    gs_im_ls = []
    for i in range(-2, 3):
        tt_ls = np.arange(t_ls[0], t_ls[-1]+1+i)

        print(tt_ls)

        Cz_0_re = []
        Cz_0_im = []
        C0_re = []
        for idt in tt_ls-1:
            C0_re.append( an_pi_bs_avg[0][idt][0] )
            Cz_0_re.append(an_pi_bs_avg[idz][idt][0])
            Cz_0_im.append(an_pi_bs_avg[idz][idt][1])

        fit_result = gs_fit_joint(tt_ls, Cz_0_re, Cz_0_im, C0_re, mom, an, L=ll)
        chi2_ls.append(fit_result.chi2 / fit_result.dof)
        gs_re_ls.append(fit_result.p['gs_re'])
        gs_im_ls.append(fit_result.p['gs_im'])

    x_ls = np.arange(-2, 3)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_sub)
    ax1.tick_params(direction='in', **ls_p)
    ax2.tick_params(direction='in', **ls_p)
    ax2.set_ylabel('chi2/dof')
    ax2.set_ylim(0, 6)
    plt.subplots_adjust(wspace=0, hspace=0)

    ax1.errorbar(x_ls, [val.mean for val in gs_re_ls], yerr=[val.sdev for val in gs_re_ls], color=color_list[0], fmt='o', label='g.s. real', **errorb)
    ax2.scatter(x_ls, chi2_ls, color='red', marker='o')
    ax2.plot(x_ls, np.ones(len(x_ls)) * 5, 'r--')
    ax1.legend(loc='upper right')
    plt.show()

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, gridspec_kw=gridspec_sub)
    ax1.tick_params(direction='in', **ls_p)
    ax2.tick_params(direction='in', **ls_p)
    ax2.set_ylabel('chi2/dof')
    ax2.set_ylim(0, 6)
    plt.subplots_adjust(wspace=0, hspace=0)

    ax1.errorbar(x_ls, [val.mean for val in gs_im_ls], yerr=[val.sdev for val in gs_im_ls], color=color_list[0], fmt='o', label='g.s. imag', **errorb)
    ax2.scatter(x_ls, chi2_ls, color='red', marker='o')
    ax2.plot(x_ls, np.ones(len(x_ls)) * 5, 'r--')
    ax1.legend(loc='upper right')
    plt.show()
   
def stability_analysis(file_path='DA_new.hdf5'):
    myfile = h5.File(file_path,'r')
    meson = 'pion'
    a_str='a12'

    if a_str == 'a06':
        an = 0.0574
        ll = 96
        N_re = 600 # same times for eliminating discretization effect
        Nz = 49 # num of z
        Nt = 23
        t_ls = np.arange(5, 13)#(3, 11)
    elif a_str == 'a09':
        an = 0.0882
        ll = 64
        N_re = 600
        Nz = 33
        Nt = 15
        t_ls = np.arange(5, 12)#(3, 11)
    elif a_str == 'a12':
        an = 0.1213
        ll = 48
        N_re = 600
        Nz = 25
        Nt = 15
        t_ls = np.arange(3, 9)#(4, 9)
    else:
        print('a input error')

    pz = 2*np.pi / (an * ll) * mom * gev_fm

    an_pi = myfile[a_str+'m130_'+meson][int(mom/2)] # conf, z, t, real/imag

    an_pi_bs = bootstrap(an_pi, N_re) # resampling

    for n_conf in range(N_re): # normalize
        for idt in range(Nt):
            for idz in range(1, Nz):
                an_pi_bs[n_conf][idz][idt][0] = an_pi_bs[n_conf][idz][idt][0] / an_pi_bs[n_conf][0][idt][0]

                an_pi_bs[n_conf][idz][idt][1] = an_pi_bs[n_conf][idz][idt][1] / an_pi_bs[n_conf][0][idt][0]

    an_pi_bs_avg = gv.dataset.avg_data(an_pi_bs, bstrap=True)

    varying = range(-2, 2)

    chi2_ls = []
    gs_re_ls = []
    gs_im_ls = []
    for i in varying:
        chi2_ls.append([])
        gs_re_ls.append([])
        gs_im_ls.append([])
        for idz in range(1, Nz):
            tt_ls = np.arange(t_ls[0], t_ls[-1]+1+i)

            Cz_0_re = []
            Cz_0_im = []
            C0_re = []
            for idt in tt_ls-1:
                C0_re.append( an_pi_bs_avg[0][idt][0] )
                Cz_0_re.append(an_pi_bs_avg[idz][idt][0])
                Cz_0_im.append(an_pi_bs_avg[idz][idt][1])

            fit_result = gs_fit_joint(tt_ls, Cz_0_re, Cz_0_im, C0_re, mom, an, L=ll)
            chi2_ls[i+2].append(fit_result.chi2 / fit_result.dof)
            gs_re_ls[i+2].append(fit_result.p['gs_re'])
            gs_im_ls[i+2].append(fit_result.p['gs_im'])

    x_ls = np.array(varying)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for x in x_ls:
        ax.errorbar(np.arange(1, Nz)+x*0.1, [val.mean for val in gs_re_ls[x+2]], [val.sdev for val in gs_re_ls[x+2]], color=color_list[x+2], label=str(x), fmt='o', **errorb)
    ax.legend(loc='upper right')
    ax.set_xlabel('z', **fs_p)
    ax.set_ylabel('real', **fs_p)
    ax.tick_params(direction='in', **ls_p)
    plt.show()
    return 

def da_data_new(a_str='a06', file_path='DA_new.hdf5'):
    myfile = h5.File(file_path,'r')
    meson = 'pion'

    if a_str == 'a06':
        an = 0.0574
        ll = 96
        N_re = 600 # same times for eliminating discretization effect
        Nz = 49 # num of z
        Nt = 23
        t_ls = np.arange(2, 12)#(5, 13)
    elif a_str == 'a09':
        an = 0.0882
        ll = 64
        N_re = 600
        Nz = 33
        Nt = 15
        t_ls = np.arange(1, 10)#(4, 11)
    elif a_str == 'a12':
        an = 0.1213
        ll = 48
        N_re = 600
        Nz = 25
        Nt = 15
        t_ls = np.arange(1, 7)#(3, 9)
    else:
        print('a input error')

    pz = 2*np.pi / (an * ll) * mom * gev_fm

    an_pi = myfile[a_str+'m130_'+meson][int(mom/2)] # conf, z, t, real/imag

    an_pi_bs = bootstrap(an_pi, N_re) # resampling

    for n_conf in range(N_re): # normalize
        for idt in range(Nt):
            for idz in range(1, Nz):
                an_pi_bs[n_conf][idz][idt][0] = an_pi_bs[n_conf][idz][idt][0] / an_pi_bs[n_conf][0][idt][0]

                an_pi_bs[n_conf][idz][idt][1] = an_pi_bs[n_conf][idz][idt][1] / an_pi_bs[n_conf][0][idt][0]

    an_pi_bs_avg = gv.dataset.avg_data(an_pi_bs, bstrap=True)

    an_pi_norm = [] # conf, z, (complex)
    print('>>> fit to extract g.s. of '+a_str)
    for n_conf in tqdm(range(400, 410)):
        an_pi_norm.append([complex(1, 0)])
        for idz in range(1, Nz):
            Cz_0_re = []
            Cz_0_im = []
            C0_re = []
            for idt in t_ls-1:
                c0 = an_pi_bs[n_conf][0][idt][0]
                c0_sdev = an_pi_bs_avg[0][idt][0].sdev
                C0_re.append( gv.gvar(c0, c0_sdev) )
                
                real = an_pi_bs[n_conf][idz][idt][0]
                real_sdev = (an_pi_bs_avg[idz][idt][0]).sdev

                imag = an_pi_bs[n_conf][idz][idt][1]
                imag_sdev = (an_pi_bs_avg[idz][idt][1]).sdev
                
                Cz_0_re.append(gv.gvar(real, real_sdev))
                Cz_0_im.append(gv.gvar(imag, imag_sdev))

            fit_result = gs_fit_joint(t_ls, Cz_0_re, Cz_0_im, C0_re, mom, an, L=ll)

            an_pi_norm[n_conf-400].append( complex(fit_result.p['gs_re'].mean, fit_result.p['gs_im'].mean) )

    da_conf_ls = []
    for n_conf in range(400, 410):
        z_in = an*np.arange(Nz)
        y_in = an_pi_norm[n_conf-400]
        da_conf_ls.append( interp_1d(z_in, y_in, z_ls_extend, method='cubic') )

    da_conf_ls = np.array(da_conf_ls) # conf, z

    return da_conf_ls


zR_dic, m_pdf_dic = pdf_zR()

mom = 8
t_plot()
# stability_analysis()
# da_conf_ls = da_data_new() # conf, z

# zR_dic, m_pdf_dic = pdf_zR()

# def renorm(da_conf_ls, key): 
#     da_re_f_ls = []
#     da_im_f_ls = []
#     da_re_ls = []
#     da_im_ls = []
#     for n_conf in range(len(da_conf_ls)):
#         temp_re_ls = []
#         temp_im_ls = []
#         for idx in range(len(da_conf_ls[0])):

#             temp_val = (da_conf_ls[n_conf][idx] / zR_dic[key][idx].mean) / ZMS_da(z_ls_extend[idx]) # no Z_hyb, divide by zR*ZMS

#             temp_re_ls.append(temp_val.real)
#             temp_im_ls.append(temp_val.imag)

#         da_re_ls.append(temp_re_ls)
#         da_im_ls.append(temp_im_ls)

#     return da_re_ls, da_im_ls


# def make_hyb(da_conf_ls, key): # for each a
#     da_re_ls, da_im_ls = renorm(da_conf_ls, key, zR_dic) 

#     hyb_re_ls = da_re_ls
#     hyb_im_ls = da_im_ls

#     return hyb_re_ls, hyb_im_ls # conf, z_ls_extend

# a12_hyb_re_ls, a12_hyb_im_ls = make_hyb(da_conf_ls, key='a=0.0574')

# lambda_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

# a12_hyb_re_ls = turn_into_lambda(a12_hyb_re_ls, z_ls_extend, lambda_ls, 0.1213, 48, mom)
# a12_hyb_im_ls = turn_into_lambda(a12_hyb_im_ls, z_ls_extend, lambda_ls, 0.1213, 48, mom)


# pz = int(mom) * 0.215
# lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

# a12_re_ro, a12_im_ro = rotate(a12_hyb_re_ls, a12_hyb_im_ls, lam_ls, pz, back=False)

# a12_hyb_re_avg = gv.dataset.avg_data(a12_hyb_re_ls, bstrap=True)
# a12_hyb_im_avg = gv.dataset.avg_data(a12_hyb_im_ls, bstrap=True)
# a12_ro_avg = gv.dataset.avg_data(a12_re_ro, bstrap=True)

# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)
# ax.errorbar(lam_ls, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)
# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-0.5, 1.25])
# ax.tick_params(direction='in', **ls_p)
# plt.show()


# %%
# def da_data_old(a_str='a09', file_path='DA_new.hdf5'):
#     myfile = h5.File(file_path,'r')
#     fix_t_ls = [8, 6, 4]
#     meson = 'pion'
#     mom = 8

#     if a_str == 'a06':
#         an = 0.0574
#         N_re = 600 # same times for eliminating discretization effect
#         Nz = 49 # num of z
#         nt = fix_t_ls[0] # fix t
#     elif a_str == 'a09':
#         an = 0.0882
#         N_re = 600
#         Nz = 33
#         nt = fix_t_ls[1]
#     elif a_str == 'a12':
#         an = 0.1213
#         N_re = 600
#         Nz = 25
#         nt = fix_t_ls[2]
#     else:
#         print('a input error')

#     an_pi = myfile[a_str+'m130_'+meson][int(mom/2)] # conf, z, t, real/imag

#     an_pi_bs = bootstrap(an_pi, N_re) # resampling

#     an_pi_norm = [] # conf, z, (complex)
#     for n_conf in range(N_re):
#         an_pi_norm.append([])
#         for nz in range(Nz): # pick three t vals to avg.
#             real = ( an_pi_bs[n_conf][nz][nt-1][0] / an_pi_bs[n_conf][0][nt-1][0] 
#             + an_pi_bs[n_conf][nz][nt][0] / an_pi_bs[n_conf][0][nt][0] 
#             + an_pi_bs[n_conf][nz][nt+1][0] / an_pi_bs[n_conf][0][nt+1][0] ) / 3

#             imag = ( an_pi_bs[n_conf][nz][nt-1][1] / an_pi_bs[n_conf][0][nt-1][0] 
#             + an_pi_bs[n_conf][nz][nt][1] / an_pi_bs[n_conf][0][nt][0] 
#             + an_pi_bs[n_conf][nz][nt+1][1] / an_pi_bs[n_conf][0][nt+1][0] ) / 3

#             an_pi_norm[n_conf].append( complex(real, imag) )

#     da_conf_ls = []
#     for n_conf in range(N_re):
#         z_in = an*np.arange(Nz)
#         y_in = an_pi_norm[n_conf]
#         da_conf_ls.append( interp_1d(z_in, y_in, z_ls_extend, method='cubic') )

#     da_conf_ls = np.array(da_conf_ls) # conf, z

#     return da_conf_ls