# %%
from head import *

def read_data(a_str, mom):
    file_path='DA_new.hdf5'
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

    # for n_conf in range(N_re): # normalization
    #     for idt in range(Nt):
    #         for idz in range(1, Nz):
    #             an_pi_bs[n_conf][idz][idt][0] = an_pi_bs[n_conf][idz][idt][0] / an_pi_bs[n_conf][0][idt][0]

    #             an_pi_bs[n_conf][idz][idt][1] = an_pi_bs[n_conf][idz][idt][1] / an_pi_bs[n_conf][0][idt][0]

    an_pi_bs_avg = gv.dataset.avg_data(an_pi_bs, bstrap=True)

    return an_pi_bs_avg, Nt, Nz, ll

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

def fcn(x, p):
    val = {}
    val['z_0_re'] = p['gs_re'] * ( 
        1 + p['b1_re'] * np.exp(- p['dE1'] * x['z_0_re']) 
        #+ p['b2_re'] * np.exp(- p['dE2'] * x['z_0_re']) 
        ) / ( 
            1 + p['a1'] * np.exp(- p['dE1'] * x['z_0_re']) 
            #+ p['a2'] * np.exp(- p['dE2'] * x['z_0_re']) 
            )

    val['z_0_im'] = p['gs_im'] * ( 
        1 + p['b1_im'] * np.exp(- p['dE1'] * x['z_0_im']) 
        #+ p['b2_im'] * np.exp(- p['dE2'] * x['z_0_im']) 
        ) / ( 
            1 + p['a1'] * np.exp(- p['dE1'] * x['z_0_im']) 
            #+ p['a2'] * np.exp(- p['dE2'] * x['z_0_im']) 
            )

    val['0'] = p['c'] * np.exp(- p['e0'] * x['0']) * ( 
        1 + p['a1'] * np.exp(- p['dE1'] * x['0']) 
        #+ p['a2'] * np.exp(- p['dE2'] * x['0']) 
        )

    return val

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
    priors['b1_re'] = gv.gvar(1, 10)
    priors['b2_re'] = gv.gvar(1, 10)
    priors['b1_im'] = gv.gvar(1, 10)
    priors['b2_im'] = gv.gvar(1, 10)
    priors['c'] = gv.gvar(0, 10)
    priors['dE1'] = gv.gvar(0, 10)#dE1
    priors['dE2'] = gv.gvar(0, 10)#dE2
    priors['e0'] = gv.gvar(0.5, 10)

    x_dic = {}
    x_dic['0'], x_dic['z_0_re'], x_dic['z_0_im'] = t_ls, t_ls, t_ls

    y_dic = {}
    y_dic['0'], y_dic['z_0_re'], y_dic['z_0_im'] = C0_re, Cz_0_re, Cz_0_im

    fit_result = lsf.nonlinear_fit(data=(x_dic, y_dic), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    print(fit_result.format(100))

    return fit_result

def gs_fit_chained(t_ls, Cz_re, Cz_im, C0_re):

    return 

def t_plot(mom, a_str, z_plot):
    an_pi_bs_avg, Nt, Nz, ll = read_data(a_str, mom)

    m_eff = [np.log(an_pi_bs_avg[0,t,0] / an_pi_bs_avg[0,t+1,0]) for t in range(Nt - 1)]

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)
    ax.errorbar(np.arange(1, Nt)[:13], [val.mean for val in m_eff][:13], yerr=[val.sdev for val in m_eff][:13], label='Local', **errorb)
    ax.set_ylim([-0.5, 1.5])
    ax.set_xlabel(t_label, **fs_p)
    ax.set_ylabel(meff_local_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig('fit_fig/local_meff.pdf', transparent=True)
    plt.show()

    m_eff_z = gv.BufferDict()
    for z in range(1, Nz):
        m_eff_z[str(z)] = [np.log(an_pi_bs_avg[z,t,0] / an_pi_bs_avg[z,t+1,0]) for t in range(Nt - 1)]

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)
    ax.errorbar(np.arange(1, Nt)[:15], [val.mean for val in m_eff_z[str(z_plot)]][:15], yerr=[val.sdev for val in m_eff_z[str(z_plot)]][:15], label='Non local, z='+str(z_plot)+', real', **errorb)
    ax.set_ylim([-0.5, 1.5])
    ax.set_xlabel(t_label, **fs_p)
    ax.set_ylabel(meff_non_local_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig('fit_fig/non_local_meff.pdf', transparent=True)
    plt.show()

def one_state_fit(t_ls, ratio):
    def fcn(x, p):
        return p['gs'] + 0 * x
    priors = gv.BufferDict()
    priors['gs'] = gv.gvar(0, 10)

    fit_result = lsf.nonlinear_fit(data=(t_ls, ratio), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    return fit_result.p['gs']

def pt2_fit(mom, a_str, t_ls, z_fit):
    if a_str == 'a06':
        an = 0.0574
    elif a_str == 'a09':
        an = 0.0882
    elif a_str == 'a12':
        an = 0.1213
    else:
        print('a input error')

    an_pi_bs_avg, Nt, Nz, ll = read_data(a_str, mom)

    C0_re = []
    Cz_0_re = []
    Cz_0_im = []

    for t in t_ls:
        idt = t-1 # t start from 0
        C0_re.append( an_pi_bs_avg[0, idt, 0] )
        Cz_0_re.append( an_pi_bs_avg[z_fit, idt, 0] / an_pi_bs_avg[0, idt, 0] ) # normalization
        Cz_0_im.append( an_pi_bs_avg[z_fit, idt, 1] / an_pi_bs_avg[0, idt, 0] )

    fit_result = gs_fit_joint(t_ls, Cz_0_re, Cz_0_im, C0_re, mom, an, L=ll)

    x_dic = {}
    x_dic['0'], x_dic['z_0_re'], x_dic['z_0_im'] = np.arange(1, Nt)[4:14], np.arange(1, Nt)[4:14], np.arange(1, Nt)[4:14]

    t_ls_ = np.arange(6, 9)
    ratio = [ (an_pi_bs_avg[z_fit, t-1, 0] / an_pi_bs_avg[0, t-1, 0]) for t in t_ls_ ]

    one_state_gs = one_state_fit(t_ls_, ratio)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes_small)
    ax.errorbar(np.arange(1, Nt)[:13], [ (an_pi_bs_avg[z_fit, idt, 0] / an_pi_bs_avg[0, idt, 0]).mean for idt in range(0, Nt-1) ][:13], [ (an_pi_bs_avg[z_fit, idt, 0] / an_pi_bs_avg[0, idt, 0]).sdev for idt in range(0, Nt-1) ][:13], label='Data, z='+str(z_fit), **errorb)

    ax.fill_between(np.arange(1, Nt)[4:14], [(val.mean+val.sdev) for val in fcn(x_dic, fit_result.p)['z_0_re']], [(val.mean-val.sdev) for val in fcn(x_dic, fit_result.p)['z_0_re']], alpha=0.5, label='Two states fit', color='green')
    gs_fit = fit_result.p['gs_re']
    ax.fill_between(np.array(t_ls), np.ones(len(t_ls))*(gs_fit.mean+gs_fit.sdev), np.ones(len(t_ls))*(gs_fit.mean-gs_fit.sdev), alpha=0.5, label='Two states fit, g.s.', color='orange' )

    ax.fill_between(t_ls_, np.ones(len(t_ls_))*(one_state_gs.mean + one_state_gs.sdev), np.ones(len(t_ls_))*(one_state_gs.mean - one_state_gs.sdev), color='red', alpha=0.7, label='One state fit' )

    ax.set_ylim([0.6, 1.2])
    ax.set_xlabel(t_label, **fs_p)
    ax.set_ylabel(ratio_label, **fs_p)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    plt.savefig('fit_fig/fit_result.pdf', transparent=True)
    plt.show()

    return fit_result


if __name__ == "__main__":
    mom=6
    a_str='a12'
    z=1
    t_ls=np.arange(4, 10)

    t_plot(mom, a_str, z_plot=z)

    # an in [0.1213, 0.0882, 0.0574]

    fit_result = pt2_fit(mom, a_str, t_ls, z_fit=z)



# %%
