# %%
from gvar import chi2
from numpy.lib.function_base import average
from head import *

# %%
class MESON_DA_HYB():
    def __init__(self, meson, mom, x_ls, x_ls_matching, y_ls_matching, extend_point, extend_fit_start, t_dic, gs_extract, fit_1, matching, fit_2_quasi, fit_2_lc, rotate, lcda_fit, constant_fit):
        self.meson = meson # str: 'pion', 'kaon', 'etas'
        self.mom = mom # int: 6, 8, 10, 12
        self.pz = mom_to_pz * mom # mom=8 corresponding to pz=1.72, used after the continuous limit
        self.lambda_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * self.mom * gev_fm ) / gev_fm # lambda = z * pz, choose a06's pz, because it has smallest range, interpolate into lambda_ls then do continuous limit
        self.x_ls = x_ls # after FT, use x as the variable
        self.x_ls_matching = x_ls_matching # for quasi
        self.y_ls_matching = y_ls_matching # for light-cone
        self.extend_point = extend_point
        self.lambda_max = self.lambda_ls[extend_point] # for lambda bigger than this max, use extend function
        self.extend_fit_start = extend_fit_start # start to fit extend func from this idx
        self.t_dic = t_dic # t range for avg or const fit
        self.gs_extract = gs_extract # whether extract gs
        self.fit_1 = fit_1 # True/False, whether fit to get continuous limit
        self.matching = matching # whether do the integrals for matching in the coordinate space
        self.fit_2_quasi = fit_2_quasi # True/False, whether fit to extend hybrid data in the coordinate space
        self.fit_2_lc = fit_2_lc # True/False, whether fit to extend hybrid data in the coordinate space
        self.rotate = rotate # True/False, whether do the rotate, for kaon should be False
        self.lcda_fit = lcda_fit # whether do the endpoints fit for lcda plot
        self.constant_fit = constant_fit

    def main(self, zR_dic):
        ### read data, prepare for renormalization ###
        ##############################################
        if self.gs_extract == True:
            da_a06_conf_ls = self.da_data('a06', 'DA_new.hdf5') # conf, z in z_ls_extend
            da_a09_conf_ls = self.da_data('a09', 'DA_new.hdf5')
            da_a12_conf_ls = self.da_data('a12', 'DA_new.hdf5')
            if not os.path.exists(self.meson+'/mom='+str(self.mom)+'/'): # for some temp save
                os.makedirs(self.meson+'/mom='+str(self.mom)+'/') 
            gv.dump(da_a06_conf_ls, self.meson+'/mom='+str(self.mom)+'/da_a06_conf_ls')
            gv.dump(da_a09_conf_ls, self.meson+'/mom='+str(self.mom)+'/da_a09_conf_ls')
            gv.dump(da_a12_conf_ls, self.meson+'/mom='+str(self.mom)+'/da_a12_conf_ls')

        da_a06_conf_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/da_a06_conf_ls')
        da_a09_conf_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/da_a09_conf_ls')
        da_a12_conf_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/da_a12_conf_ls')

        ### renorm -> continuous limit ###
        ##########################################
        a_hyb_re_ls, a_hyb_im_ls, a_hyb_ro_re_avg, a_hyb_ro_im_avg = self.da_hyb(da_a06_conf_ls, da_a09_conf_ls, da_a12_conf_ls, zR_dic) 

        ### matching in the coordinate space ###
        ##########################################
        if self.matching == True:
            lc_re_ls, lc_im_ls = self.inv_matching_coor(a_hyb_re_ls, a_hyb_im_ls)
            gv.dump(lc_re_ls, self.meson+'/mom='+str(self.mom)+'/lc_re_ls')
            gv.dump(lc_im_ls, self.meson+'/mom='+str(self.mom)+'/lc_im_ls')

        lc_re_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/lc_re_ls')
        lc_im_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/lc_im_ls')

        lc_re_avg = gv.dataset.avg_data(lc_re_ls, bstrap=True)
        lc_im_avg = gv.dataset.avg_data(lc_im_ls, bstrap=True)

        ### extrapolation -> FT -> quasi / lcda ###
        ##########################################
        quasi_ls = self.extrapolation_FT(a_hyb_re_ls, a_hyb_im_ls, a_hyb_ro_re_avg, a_hyb_ro_im_avg, mode='quasi', extrapolate=self.fit_2_quasi) # hyb_quasi.shape = (N_conf, len(x_ls))

        lcda_ls = self.extrapolation_FT(lc_re_ls, lc_im_ls, lc_re_avg, lc_im_avg, mode='lcda', extrapolate=self.fit_2_lc) # lic_da.shape = (N_conf, len(x_ls))

        ## gaussian filter
        #lcda_ls = gaussian_filter(lcda_ls, sigma=1)
        
        quasi_avg = gv.dataset.avg_data(quasi_ls, bstrap=True)
        lcda_avg = gv.dataset.avg_data(lcda_ls, bstrap=True)

        return quasi_avg, lcda_avg   
    
    ### read data and normalize, interpolate, prepare for renormalization ###
    ###########################################
    def da_data(self, a_str, file_path): 
        myfile = h5.File(file_path,'r')
        if a_str == 'a06':
            an = 0.0574
            N_re = 600 # resampling times in bootstrap
            Nz = 49 # num of z
        elif a_str == 'a09':
            an = 0.0882
            N_re = 600
            Nz = 33
        elif a_str == 'a12':
            an = 0.1213
            N_re = 600
            Nz = 25
        else:
            print('a input error')

        t_ls = self.t_dic[a_str] # choose some t to average

        an_pi = myfile[a_str+'m130_'+self.meson][int(self.mom/2)] # conf, z, t, real/imag

        an_pi_bs = bootstrap(an_pi, N_re) # resample

        an_pi_bs_avg = gv.dataset.avg_data(an_pi_bs, bstrap=True)

        an_pi_norm = [] # conf, z (complex number)
        print('>>> extracting g.s. from t plot of '+a_str+': ')
        for n_conf in tqdm(range(N_re)):
            an_pi_norm.append([1+0j])
            for idz in range(1, Nz): 
                re_ls = []
                im_ls = []
                re_avg_ls = []
                im_avg_ls = []
                for t in t_ls:
                    idt = t - 1
                    re_ls.append( an_pi_bs[n_conf][idz][idt][0] / an_pi_bs[n_conf][0][idt][0] ) # normalization
                    im_ls.append( an_pi_bs[n_conf][idz][idt][1] / an_pi_bs[n_conf][0][idt][0] )

                    re_avg_ls.append( an_pi_bs_avg[idz][idt][0] / an_pi_bs_avg[0][idt][0] )
                    im_avg_ls.append( an_pi_bs_avg[idz][idt][1] / an_pi_bs_avg[0][idt][0] )
                
                real = np.average(re_ls)
                imag = np.average(im_ls)

                ######## constant fit ########
                if self.constant_fit == True:
                    def fcn(x, p):
                        return p['ratio'] + 0*x

                    priors = gv.BufferDict()
                    priors['ratio'] = gv.gvar(0, 10)

                    y_re = add_sdev(re_ls, re_avg_ls)
                    y_im = add_sdev(im_ls, im_avg_ls)

                    fit_result_re = lsf.nonlinear_fit(data=(np.array(t_ls), y_re), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

                    fit_result_im = lsf.nonlinear_fit(data=(np.array(t_ls), y_im), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

                    real = fit_result_re.p['ratio'].mean
                    imag = fit_result_im.p['ratio'].mean

                an_pi_norm[n_conf].append( complex(real, imag) )

        da_conf_ls = [] # interpolate into z_ls_extend, waiting for renormalization
        for n_conf in range(N_re):
            z_in = an*np.arange(Nz)
            y_in = an_pi_norm[n_conf]
            da_conf_ls.append( interp_1d(z_in, y_in, z_ls_extend, method='cubic') )

        da_conf_ls = np.array(da_conf_ls) # conf, z in z_ls_extend

        return da_conf_ls

    def da_hyb(self, da_a06_conf_ls, da_a09_conf_ls, da_a12_conf_ls, zR_dic):
        x_ls = self.x_ls # after FT, for quasi

        def renorm(da_conf_ls, key): 
            da_re_ls = []
            da_im_ls = []
            for n_conf in range(len(da_conf_ls)):
                temp_re_ls = []
                temp_im_ls = []
                for idx in range(len(da_conf_ls[0])):
                    temp_val = (da_conf_ls[n_conf][idx] / zR_dic[key][idx].mean) / ZMS_da(z_ls_extend[idx]) # divide by zR*ZMS
                    temp_re_ls.append(temp_val.real)
                    temp_im_ls.append(temp_val.imag)

                da_re_ls.append(temp_re_ls)
                da_im_ls.append(temp_im_ls)

            return da_re_ls, da_im_ls
        
        a06_hyb_re_ls, a06_hyb_im_ls = renorm(da_a06_conf_ls, key='a=0.0574')
        a09_hyb_re_ls, a09_hyb_im_ls = renorm(da_a09_conf_ls, key='a=0.0882')
        a12_hyb_re_ls, a12_hyb_im_ls = renorm(da_a12_conf_ls, key='a=0.1213')

        ## turn into lambda axis
        ## from z_ls_extend to lambda_ls(z_ls_da)
        a06_hyb_re_ls = turn_into_lambda(a06_hyb_re_ls, z_ls_extend, self.lambda_ls, 0.0574, 96, self.mom) 
        a06_hyb_im_ls = turn_into_lambda(a06_hyb_im_ls, z_ls_extend, self.lambda_ls, 0.0574, 96, self.mom)
        a09_hyb_re_ls = turn_into_lambda(a09_hyb_re_ls, z_ls_extend, self.lambda_ls, 0.0882, 64, self.mom)
        a09_hyb_im_ls = turn_into_lambda(a09_hyb_im_ls, z_ls_extend, self.lambda_ls, 0.0882, 64, self.mom)
        a12_hyb_re_ls = turn_into_lambda(a12_hyb_re_ls, z_ls_extend, self.lambda_ls, 0.1213, 48, self.mom)
        a12_hyb_im_ls = turn_into_lambda(a12_hyb_im_ls, z_ls_extend, self.lambda_ls, 0.1213, 48, self.mom)

        gv.dump(a06_hyb_re_ls, self.meson+'/mom='+str(self.mom)+'/a06_hyb_re_ls')
        gv.dump(a06_hyb_im_ls, self.meson+'/mom='+str(self.mom)+'/a06_hyb_im_ls')
        gv.dump(a09_hyb_re_ls, self.meson+'/mom='+str(self.mom)+'/a09_hyb_re_ls')
        gv.dump(a09_hyb_im_ls, self.meson+'/mom='+str(self.mom)+'/a09_hyb_im_ls')
        gv.dump(a12_hyb_re_ls, self.meson+'/mom='+str(self.mom)+'/a12_hyb_re_ls')
        gv.dump(a12_hyb_im_ls, self.meson+'/mom='+str(self.mom)+'/a12_hyb_im_ls')

        ## avg over all conf
        a06_hyb_re_avg = gv.dataset.avg_data(a06_hyb_re_ls, bstrap=True)
        a06_hyb_im_avg = gv.dataset.avg_data(a06_hyb_im_ls, bstrap=True)
        a09_hyb_re_avg = gv.dataset.avg_data(a09_hyb_re_ls, bstrap=True)
        a09_hyb_im_avg = gv.dataset.avg_data(a09_hyb_im_ls, bstrap=True)
        a12_hyb_re_avg = gv.dataset.avg_data(a12_hyb_re_ls, bstrap=True)
        a12_hyb_im_avg = gv.dataset.avg_data(a12_hyb_im_ls, bstrap=True)

        ### fit to get continuous limit ###
        ###################################
        if self.fit_1 == True:
            ## eliminate discretization effect again
            a_hyb_re_ls = [] # a to 0
            a_hyb_im_ls = []
            print('>>> fitting f0 for configs')

            f1_ls = []#!#
            f2_ls = []
            chi2_ls = []

            for n_conf in tqdm(range(len(a06_hyb_re_ls))): # all a have same num of configs
                f1_ls.append([])
                f2_ls.append([])

                a_hyb_re_ls.append([])
                a_hyb_im_ls.append([])
                
                a06_hyb_re_gv = add_sdev(a06_hyb_re_ls[n_conf], a06_hyb_re_avg) # shape = len(self.lambda_ls)
                a06_hyb_im_gv = add_sdev(a06_hyb_im_ls[n_conf], a06_hyb_im_avg)
                a09_hyb_re_gv = add_sdev(a09_hyb_re_ls[n_conf], a09_hyb_re_avg) 
                a09_hyb_im_gv = add_sdev(a09_hyb_im_ls[n_conf], a09_hyb_im_avg)
                a12_hyb_re_gv = add_sdev(a12_hyb_re_ls[n_conf], a12_hyb_re_avg) 
                a12_hyb_im_gv = add_sdev(a12_hyb_im_ls[n_conf], a12_hyb_im_avg)

                for idx in range(len(self.lambda_ls)): ## for each conf, at each lambda, fit to eliminate discretization effect
                    def fit_f0(a_ls, hyb_a):
                        def fcn(x, p): #!#
                            return p['f0'] + p['f1'] * (x) #+ p['f2'] * (x**2)
                        priors = gv.BufferDict()
                        priors['f0'] = gv.gvar(1, 10)
                        priors['f1'] = gv.gvar(1, 10)
                        priors['f2'] = gv.gvar(1, 10)

                        fit_result = lsf.nonlinear_fit(data=(a_ls, hyb_a), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

                        return fit_result.p['f0'].mean

                    a_ls = np.array([0.0574, 0.0882, 0.1213])
                    hyb_a_re = [a06_hyb_re_gv[idx], a09_hyb_re_gv[idx], a12_hyb_re_gv[idx]]
                    hyb_a_im = [a06_hyb_im_gv[idx], a09_hyb_im_gv[idx], a12_hyb_im_gv[idx]]

                    a_hyb_re_ls[n_conf].append(fit_f0(a_ls, hyb_a_re)) 
                    a_hyb_im_ls[n_conf].append(fit_f0(a_ls, hyb_a_im))

            gv.dump(a_hyb_re_ls, self.meson+'/mom='+str(self.mom)+'/a_hyb_re_ls')
            gv.dump(a_hyb_im_ls, self.meson+'/mom='+str(self.mom)+'/a_hyb_im_ls')

        a_hyb_re_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/a_hyb_re_ls') # shape = (N_conf, N_lambda)
        a_hyb_im_ls = gv.load(self.meson+'/mom='+str(self.mom)+'/a_hyb_im_ls')

        ### rotate, eliminate imag part, rotate back ###
        ################################################
        if self.rotate == True:
            a_hyb_re_ls, a_hyb_im_ls = rotate(a_hyb_re_ls, a_hyb_im_ls, self.lambda_ls, back=True)

        a_hyb_ro_re_avg = gv.dataset.avg_data(a_hyb_re_ls, bstrap=True)
        a_hyb_ro_im_avg = gv.dataset.avg_data(a_hyb_im_ls, bstrap=True)

        return a_hyb_re_ls, a_hyb_im_ls, a_hyb_ro_re_avg, a_hyb_ro_im_avg

    def inv_matching_coor(self, a_hyb_re_ls, a_hyb_im_ls):
        lc_re_ls = []
        lc_im_ls = []

        delta_l = 0.000001 # regulater

        print('>>> integrating for matching')
        for n_conf in tqdm(range(len(a_hyb_re_ls))):
            lc_re_ls.append([])
            lc_im_ls.append([])

            quasi = [1+0j]
            for idl in range(len(a_hyb_re_ls[n_conf])):
                quasi.append( a_hyb_re_ls[n_conf][idl] + 1j * a_hyb_im_ls[n_conf][idl] )

            lam_ls_0 = np.insert(self.lambda_ls, 0, 0)
            h_tilde = interpolate.interp1d(lam_ls_0, quasi, kind='cubic')

            for idl in range(len(self.lambda_ls)):
                lam = self.lambda_ls[idl]
                z = lam / self.pz # here z has GeV^-1

                lc_re = a_hyb_re_ls[n_conf][idl]
                lc_im = a_hyb_im_ls[n_conf][idl]
 
                part1 = alphas_cf_div_2pi * 1/2 * (f_matching(z**2, mu**2) - 3) * h_tilde(lam)

                def fp_2_re(lamp, lam):
                    z = lam / self.pz # here z has GeV^-1
                    res = alphas_cf_div_2pi / lam * (-1-f_matching(z**2, mu**2)) * (lamp/(lam-lamp)) * ( 
                        ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                        - 2 * h_tilde(lam)
                    )
                    return res.real

                def fp_2_im(lamp, lam):
                    z = lam / self.pz # here z has GeV^-1
                    res = alphas_cf_div_2pi / lam * (-1-f_matching(z**2, mu**2)) * (lamp/(lam-lamp)) * ( 
                        ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                        - 2 * h_tilde(lam)
                    )
                    return res.imag

                part2 = integrate.quad(fp_2_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_2_im, 0, lam-delta_l, args=lam)[0]

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
                    z = lam / self.pz # here z has GeV^-1
                    res = alphas_cf_div_2pi / (1j * (lam**2)) * ( 1 - np.exp(-1j * (lam - lamp)) ) * (3-f_matching(z**2, mu**2)) * h_tilde(lamp)
                    return res.real

                def fp_4_im(lamp, lam):
                    z = lam / self.pz # here z has GeV^-1
                    res = alphas_cf_div_2pi / (1j * (lam**2)) * ( 1 - np.exp(-1j * (lam - lamp)) ) * (3-f_matching(z**2, mu**2)) * h_tilde(lamp)
                    return res.imag
                part4 = integrate.quad(fp_4_re, 0, lam, args=lam)[0] + 1j * integrate.quad(fp_4_im, 0, lam, args=lam)[0]

                lc_re = lc_re - (part1+part2+part3+part4).real
                lc_re_ls[n_conf].append(lc_re)

                lc_im = lc_im - (part1+part2+part3+part4).imag
                lc_im_ls[n_conf].append(lc_im)

        return lc_re_ls, lc_im_ls

    def extrapolation_FT(self, a_hyb_re_ls, a_hyb_im_ls, a_hyb_ro_re_avg, a_hyb_ro_im_avg, mode, extrapolate):
        ### fit z to get DA hybrid over the whole z-axis, then do ft to momentum space ###
        ###########################################
        if extrapolate == True:
            if mode == 'quasi':
                extend_length = extend_length_quasi
            elif mode == 'lcda':
                extend_length = extend_length_lc

            hyb_re_ft = []
            hyb_im_ft = []
            hyb_complex = []
            
            lambda_delta = self.lambda_ls[1] - self.lambda_ls[0]
            lambda_max = self.lambda_max # when lambda bigger than this max, use the extend function

            print('>>> fitting to extend DA hybrid to all lambda, then do ft')
            if self.meson == 'pion':    
                pn_conf_ls = []

            for n_conf in tqdm(range(len(a_hyb_re_ls))):
                a_hyb_re_gv = add_sdev(a_hyb_re_ls[n_conf], a_hyb_ro_re_avg)
                a_hyb_im_gv = add_sdev(a_hyb_im_ls[n_conf], a_hyb_ro_im_avg)

                priors = gv.BufferDict()
                priors['c1'] = gv.gvar(1, 10)
                priors['c2'] = gv.gvar(1, 10)
                priors['x0'] = gv.gvar(1, 10)
                if self.meson == 'pion':
                    priors['n'] = gv.gvar(1, 10)
                elif self.meson == 'kaon':
                    priors['n1'] = gv.gvar(1, 10)
                    priors['n2'] = gv.gvar(1, 10)

                def fcn(x, p):
                    if mode == 'quasi':
                        no_e = False
                    elif mode == 'lcda':
                        no_e = True
                    val = {}
                    if self.meson == 'pion':
                        if no_e == False:
                            val['re'] = ( p['c1']/(x['re']**p['n']) * np.cos(np.pi/2 * p['n']) + p['c2']/(x['re']**p['n']) * np.cos(x['re'] - np.pi/2 * p['n']) ) * np.exp(-x['re'] / 150) #e# ###

                            val['im'] = - ( p['c1']/(x['im']**p['n']) * np.sin(np.pi/2 * p['n']) + p['c2']/(x['im']**p['n']) * np.sin(x['im'] - np.pi/2 * p['n']) ) * np.exp(-x['im'] / 150)

                        elif no_e == True:
                            val['re'] = ( p['c1']/(x['re']**p['n']) * np.cos(np.pi/2 * p['n']) + p['c2']/(x['re']**p['n']) * np.cos(x['re'] - np.pi/2 * p['n']) ) #e# ###

                            val['im'] = - ( p['c1']/(x['im']**p['n']) * np.sin(np.pi/2 * p['n']) + p['c2']/(x['im']**p['n']) * np.sin(x['im'] - np.pi/2 * p['n']) ) 

                    elif self.meson == 'kaon':
                        val['re'] = ( p['c1']/(x['re']**p['n1']) * np.cos(np.pi/2 * p['n1']) + p['c2']/(x['re']**p['n2']) * np.cos(x['re'] - np.pi/2 * p['n2']) ) * np.exp(-x['re'] / p['x0'])

                        val['im'] = - ( p['c1']/(x['im']**p['n1']) * np.sin(np.pi/2 * p['n1']) + p['c2']/(x['im']**p['n2']) * np.sin(x['im'] - np.pi/2 * p['n2']) ) * np.exp(-x['im'] / p['x0'])

                    return val

                lambda_dic = {}
                lambda_dic['re'] = np.array(self.lambda_ls[self.extend_fit_start:])
                lambda_dic['im'] = np.array(self.lambda_ls[self.extend_fit_start:])

                da_dic = {}
                da_dic['re'] = a_hyb_re_gv[self.extend_fit_start:]
                da_dic['im'] = a_hyb_im_gv[self.extend_fit_start:]

                fit_result = lsf.nonlinear_fit(data=(lambda_dic, da_dic), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

                if self.meson == 'pion':  
                    pn_conf_ls.append(fit_result.p['n'].mean)

                hyb_conf = [] # after extend, combine real and imag to do ft
                lam_ls_ex = [] # lambda list after extend
            
                ### piece together to form the renorm_da in the whole coordinate space ###
                ###########################################
                for lam in np.arange(-lambda_delta*extend_length, -lambda_max-0.001, lambda_delta): # use extend func, complex conjugate
                    lam_ls_ex.append(lam)
                    lam_dic = {'re':-lam, 'im':-lam}
                    hyb_conf.append( complex(fcn(lam_dic, fit_result.p)['re'].mean, -fcn(lam_dic, fit_result.p)['im'].mean) )

                for idx in range(len(self.lambda_ls[:self.extend_point])+1): # use data list, complex conjugate
                    lam_ls_ex.append(-self.lambda_max + idx*lambda_delta)
                    hyb_conf.append( complex(a_hyb_re_ls[n_conf][self.extend_point-idx], -a_hyb_im_ls[n_conf][self.extend_point-idx]) )

                lam_ls_ex.append(0) # zero point
                hyb_conf.append(1+0j)

                for idx in range(len(self.lambda_ls[:self.extend_point])+1): # use data list
                    lam_ls_ex.append(self.lambda_ls[idx])
                    hyb_conf.append( complex(a_hyb_re_ls[n_conf][idx], a_hyb_im_ls[n_conf][idx]) )

                for lam in np.arange(lambda_max+lambda_delta, lambda_delta*extend_length+0.001, lambda_delta): # use extend func
                    lam_ls_ex.append(lam)
                    lam_dic = {'re':lam, 'im':lam}
                    hyb_conf.append( complex(fcn(lam_dic, fit_result.p)['re'].mean, fcn(lam_dic, fit_result.p)['im'].mean) )

                hyb_conf_re_ft = []
                hyb_conf_im_ft = []


                ### FT ###
                ##########
                for x in self.x_ls:
                    ft_res = sum_ft(lam_ls_ex, hyb_conf, lambda_delta, x)
                    hyb_conf_re_ft.append(ft_res.real)
                    hyb_conf_im_ft.append(ft_res.imag)

                hyb_complex.append(hyb_conf)
                hyb_re_ft.append(hyb_conf_re_ft)
                hyb_im_ft.append(hyb_conf_im_ft)

            if self.meson == 'pion':  
                avg_n = gv.dataset.avg_data(pn_conf_ls, bstrap=True)
                print('Average n over all configs is '+str(avg_n.mean))

            gv.dump(lam_ls_ex, self.meson+'/mom='+str(self.mom)+'/lam_ls_ex_'+mode) 
            gv.dump(hyb_complex, self.meson+'/mom='+str(self.mom)+'/hyb_complex_'+mode) # before ft
            gv.dump(hyb_re_ft, self.meson+'/mom='+str(self.mom)+'/hyb_re_ft_'+mode) # shape = (N_re, len(x_ls))
            gv.dump(hyb_im_ft, self.meson+'/mom='+str(self.mom)+'/hyb_im_ft_'+mode)

        lam_ls_ex = gv.load(self.meson+'/mom='+str(self.mom)+'/lam_ls_ex_'+mode)
        hyb_complex = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_complex_'+mode)
        hyb_re = []
        hyb_im = []
        for idx in range(len(hyb_complex)):
            hyb_re.append( [val.real for val in hyb_complex[idx]] )
            hyb_im.append( [val.imag for val in hyb_complex[idx]] )
        hyb_re_ft = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_re_ft_'+mode)
        hyb_im_ft = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_im_ft_'+mode)

        return hyb_re_ft

    # def extrapolation_FT(self, a_hyb_re_ls, a_hyb_im_ls, a_hyb_ro_re_avg, a_hyb_ro_im_avg, mode, extrapolate):
    #     ### fit z to get DA hybrid over the whole z-axis, then do ft to momentum space ###
    #     ###########################################
    #     if extrapolate == True:
    #         if mode == 'quasi':
    #             extend_length = extend_length_quasi
    #         elif mode == 'lcda':
    #             extend_length = extend_length_lc

    #         hyb_re_ft = []
    #         hyb_im_ft = []
    #         hyb_complex = []
            
    #         lambda_delta = self.lambda_ls[1] - self.lambda_ls[0]
    #         lambda_max = self.lambda_max # when lambda bigger than this max, use the extend function

    #         print('>>> fitting to extend DA hybrid to all lambda, then do ft')
    #         if self.meson == 'pion':    
    #             pn_conf_ls = []

    #         def _process(params):
    #             a_hyb_re_, a_hyb_im_ = params
    #             a_hyb_re_gv = add_sdev(a_hyb_re_, a_hyb_ro_re_avg)
    #             a_hyb_im_gv = add_sdev(a_hyb_im_, a_hyb_ro_im_avg)

    #             priors = gv.BufferDict()
    #             priors['c1'] = gv.gvar(1, 10)
    #             priors['c2'] = gv.gvar(1, 10)
    #             priors['x0'] = gv.gvar(1, 10)
    #             if self.meson == 'pion':
    #                 priors['n'] = gv.gvar(1, 10)
    #             elif self.meson == 'kaon':
    #                 priors['n1'] = gv.gvar(1, 10)
    #                 priors['n2'] = gv.gvar(1, 10)

    #             def fcn(x, p):
    #                 if mode == 'quasi':
    #                     no_e = False
    #                 elif mode == 'lcda':
    #                     no_e = True
    #                 val = {}
    #                 if self.meson == 'pion':
    #                     if no_e == False:
    #                         val['re'] = ( p['c1']/(x['re']**p['n']) * np.cos(np.pi/2 * p['n']) + p['c2']/(x['re']**p['n']) * np.cos(x['re'] - np.pi/2 * p['n']) ) * np.exp(-x['re'] / 150) #e# ###

    #                         val['im'] = - ( p['c1']/(x['im']**p['n']) * np.sin(np.pi/2 * p['n']) + p['c2']/(x['im']**p['n']) * np.sin(x['im'] - np.pi/2 * p['n']) ) * np.exp(-x['im'] / 150)

    #                     elif no_e == True:
    #                         val['re'] = ( p['c1']/(x['re']**p['n']) * np.cos(np.pi/2 * p['n']) + p['c2']/(x['re']**p['n']) * np.cos(x['re'] - np.pi/2 * p['n']) ) #e# ###

    #                         val['im'] = - ( p['c1']/(x['im']**p['n']) * np.sin(np.pi/2 * p['n']) + p['c2']/(x['im']**p['n']) * np.sin(x['im'] - np.pi/2 * p['n']) ) 

    #                 elif self.meson == 'kaon':
    #                     val['re'] = ( p['c1']/(x['re']**p['n1']) * np.cos(np.pi/2 * p['n1']) + p['c2']/(x['re']**p['n2']) * np.cos(x['re'] - np.pi/2 * p['n2']) ) * np.exp(-x['re'] / p['x0'])

    #                     val['im'] = - ( p['c1']/(x['im']**p['n1']) * np.sin(np.pi/2 * p['n1']) + p['c2']/(x['im']**p['n2']) * np.sin(x['im'] - np.pi/2 * p['n2']) ) * np.exp(-x['im'] / p['x0'])

    #                 return val

    #             lambda_dic = {}
    #             lambda_dic['re'] = np.array(self.lambda_ls[self.extend_fit_start:])
    #             lambda_dic['im'] = np.array(self.lambda_ls[self.extend_fit_start:])

    #             da_dic = {}
    #             da_dic['re'] = a_hyb_re_gv[self.extend_fit_start:]
    #             da_dic['im'] = a_hyb_im_gv[self.extend_fit_start:]

    #             fit_result = lsf.nonlinear_fit(data=(lambda_dic, da_dic), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    #             if self.meson == 'pion':  
    #                 pn_conf_ls.append(fit_result.p['n'].mean)

    #             hyb_conf = [] # after extend, combine real and imag to do ft
    #             lam_ls_ex = [] # lambda list after extend
            
    #             ### piece together to form the renorm_da in the whole coordinate space ###
    #             ###########################################
    #             for lam in np.arange(-lambda_delta*extend_length, -lambda_max-0.001, lambda_delta): # use extend func, complex conjugate
    #                 lam_ls_ex.append(lam)
    #                 lam_dic = {'re':-lam, 'im':-lam}
    #                 hyb_conf.append( complex(fcn(lam_dic, fit_result.p)['re'].mean, -fcn(lam_dic, fit_result.p)['im'].mean) )

    #             for idx in range(len(self.lambda_ls[:self.extend_point])+1): # use data list, complex conjugate
    #                 lam_ls_ex.append(-self.lambda_max + idx*lambda_delta)
    #                 hyb_conf.append( complex(a_hyb_re_[self.extend_point-idx], -a_hyb_im_[self.extend_point-idx]) )

    #             lam_ls_ex.append(0) # zero point
    #             hyb_conf.append(1+0j)

    #             for idx in range(len(self.lambda_ls[:self.extend_point])+1): # use data list
    #                 lam_ls_ex.append(self.lambda_ls[idx])
    #                 hyb_conf.append( complex(a_hyb_re_[idx], a_hyb_im_[idx]) )

    #             for lam in np.arange(lambda_max+lambda_delta, lambda_delta*extend_length+0.001, lambda_delta): # use extend func
    #                 lam_ls_ex.append(lam)
    #                 lam_dic = {'re':lam, 'im':lam}
    #                 hyb_conf.append( complex(fcn(lam_dic, fit_result.p)['re'].mean, fcn(lam_dic, fit_result.p)['im'].mean) )

    #             hyb_conf_re_ft = []
    #             hyb_conf_im_ft = []

    #             ### FT ###
    #             ##########
    #             for x in self.x_ls:
    #                 hyb_conf_re_ft.append(sum_ft(lam_ls_ex, hyb_conf, lambda_delta, x).real)
    #                 hyb_conf_im_ft.append(sum_ft(lam_ls_ex, hyb_conf, lambda_delta, x).imag)

    #             return hyb_conf, hyb_re_ft, hyb_im_ft, lam_ls_ex

    #         thread_pool = ThreadPool(multiprocessing.cpu_count())

    #         res= thread_pool.map(_process, zip(a_hyb_re_ls, a_hyb_im_ls))
    #         hyb_complex, hyb_re_ft, hyb_im_ft, lam_ls_ = tuple(zip(*res))

    #         if self.meson == 'pion':  
    #             avg_n = gv.dataset.avg_data(pn_conf_ls, bstrap=True)
    #             print('Average n over all configs is '+str(avg_n.mean))

    #         gv.dump(lam_ls_[0], self.meson+'/mom='+str(self.mom)+'/lam_ls_ex_'+mode) 
    #         gv.dump(hyb_complex, self.meson+'/mom='+str(self.mom)+'/hyb_complex_'+mode) # before ft
    #         gv.dump(hyb_re_ft, self.meson+'/mom='+str(self.mom)+'/hyb_re_ft_'+mode) # shape = (N_re, len(x_ls))
    #         gv.dump(hyb_im_ft, self.meson+'/mom='+str(self.mom)+'/hyb_im_ft_'+mode)

    #     lam_ls_ex = gv.load(self.meson+'/mom='+str(self.mom)+'/lam_ls_ex_'+mode)
    #     hyb_complex = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_complex_'+mode)
    #     hyb_re = []
    #     hyb_im = []
    #     for idx in range(len(hyb_complex)):
    #         hyb_re.append( [val.real for val in hyb_complex[idx]] )
    #         hyb_im.append( [val.imag for val in hyb_complex[idx]] )
    #     hyb_re_ft = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_re_ft_'+mode)
    #     hyb_im_ft = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_im_ft_'+mode)

    #     return hyb_re_ft