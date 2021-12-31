# %%
from head import *

# %%
### read data and normalize, interpolate, prepare for renormalization ###
###########################################
class DATA_PRE():
    def __init__(self, file_path, meson, const_fit):
        self.file_path = file_path
        self.meson = meson
        self.const_fit = const_fit

    def main(self, mom, t_dic):
        self.mom = mom
        self.t_dic = t_dic # t range for avg or const fit

        da_an_dic = {}
        for a_str in ['a06', 'a09', 'a12']:
            an_meson_norm = self.read_data(self.file_path, a_str)
            da_an_dic[a_str] = self.interp(a_str, an_meson_norm, method='cubic')
        return da_an_dic

    def paras(self, a_str):
        N_re = 600 # resampling times in bootstrap
        if a_str == 'a06':
            an = 0.0574
            Nz = 49
        elif a_str == 'a09':
            an = 0.0882
            Nz = 33
        elif a_str == 'a12':
            an = 0.1213
            Nz = 25
        else:
            print('a_str input error')

        t_ls = self.t_dic[a_str]

        return an, Nz, N_re, t_ls

    def constant_fit(self, t_ls, y_re, y_im):
        def fcn(x, p):
            return p['ratio'] + 0*x

        priors = gv.BufferDict()
        priors['ratio'] = gv.gvar(0, 10)

        fit_result_re = lsf.nonlinear_fit(data=(np.array(t_ls), y_re), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

        fit_result_im = lsf.nonlinear_fit(data=(np.array(t_ls), y_im), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

        real = fit_result_re.p['ratio'].mean
        imag = fit_result_im.p['ratio'].mean
        return real, imag

    def read_data(self, file_path, a_str): # normalization and extract g.s.
        myfile = h5.File(file_path,'r')
        an_meson = myfile[a_str+'m130_'+self.meson][int(self.mom/2)] # conf, z, t, real/imag

        an, Nz, N_re, t_ls = self.paras(a_str)

        an_meson_bs = bootstrap(an_meson, self.meson, a_str) # resample
        an_meson_bs_avg = gv.dataset.avg_data(an_meson_bs, bstrap=True)

        an_meson_norm = [] # conf, z (complex number)
        print('>>> extracting g.s. from t plot of '+a_str+' at mom '+str(self.mom))
        for n_conf in tqdm(range(N_re)):
            an_meson_norm.append([1+0j])
            for idz in range(1, Nz): 
                re_ls = []
                im_ls = []
                re_avg_ls = []
                im_avg_ls = []
                for t in t_ls:
                    idt = t - 1
                    re_ls.append( an_meson_bs[n_conf][idz][idt][0] / an_meson_bs[n_conf][0][idt][0] ) # normalization
                    im_ls.append( an_meson_bs[n_conf][idz][idt][1] / an_meson_bs[n_conf][0][idt][0] )

                    re_avg_ls.append( an_meson_bs_avg[idz][idt][0] / an_meson_bs_avg[0][idt][0] )
                    im_avg_ls.append( an_meson_bs_avg[idz][idt][1] / an_meson_bs_avg[0][idt][0] )
                
                real = np.average(re_ls)
                imag = np.average(im_ls)

                ######## constant fit ########
                if self.const_fit == True:
                    y_re = add_sdev(re_ls, re_avg_ls)
                    y_im = add_sdev(im_ls, im_avg_ls)

                    real, imag = self.constant_fit(t_ls, y_re, y_im)

                an_meson_norm[n_conf].append( complex(real, imag) )

        return an_meson_norm

    def interp(self, a_str, an_meson_norm, method):
        an, Nz, N_re, t_ls = self.paras(a_str)
        da_conf_ls = [] # interpolate into z_ls_extend, waiting for renormalization
        for n_conf in range(N_re):
            z_in = an*np.arange(Nz)
            y_in = an_meson_norm[n_conf]
            da_conf_ls.append( interp_1d(z_in, y_in, z_ls_extend, method) )

        da_conf_ls = np.array(da_conf_ls) # conf, z in z_ls_extend

        return da_conf_ls

class RENORMALIZATION():
    def __init__(self, zR_dic):
        self.zR_dic = zR_dic

    def main(self, data_ls): # data_ls should match with the key_ls
        key_ls = ['a=0.0574', 'a=0.0882', 'a=0.1213']
        renorm_data_re_ls = []
        renorm_data_im_ls = []
        for id in range(len(data_ls)):
            da_conf_ls = data_ls[id]
            key = key_ls[id]
            re, im = self.renorm(da_conf_ls, key)
            renorm_data_re_ls.append(re)
            renorm_data_im_ls.append(im)
        return renorm_data_re_ls, renorm_data_im_ls

    def renorm(self, da_conf_ls, key): 
        da_re_ls = []
        da_im_ls = []
        for n_conf in range(len(da_conf_ls)):
            temp_re_ls = []
            temp_im_ls = []
            for idx in range(len(da_conf_ls[0])):
                temp_val = (da_conf_ls[n_conf][idx] / self.zR_dic[key][idx].mean) / ZMS_da(z_ls_extend[idx]) # divide by zR*ZMS
                temp_re_ls.append(temp_val.real)
                temp_im_ls.append(temp_val.imag)

            da_re_ls.append(temp_re_ls)
            da_im_ls.append(temp_im_ls)

        return da_re_ls, da_im_ls

class CONTINUUM_LIMIT():
    def __init__(self, a_str_ls, rotate):
        self.a_str_ls = a_str_ls
        self.rotate = rotate
        self.a_ls = np.array([0.0574, 0.0882, 0.1213])

    def main(self, mom, renorm_data_re_ls, renorm_data_im_ls):
        self.mom = mom
        self.lambda_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * self.mom * gev_fm ) / gev_fm # lambda = z * pz, choose a06's pz, because it has smallest range, interpolate into lambda_ls then do continuous limit

        ## turn into lambda axis
        ## from z_ls_extend to lambda_ls(z_ls_da)
        re_lam_ls = [] # shape = (N(a_str), N_conf, len(lam_ls))
        im_lam_ls = []
        re_gv_lam_ls = []
        im_gv_lam_ls = []
        for id in range(len(self.a_str_ls)): # each an
            a_str = self.a_str_ls[id]
            an, Nt = self.paras(a_str)
            re = turn_into_lambda(renorm_data_re_ls[id], z_ls_extend, self.lambda_ls, an, Nt, self.mom)
            im = turn_into_lambda(renorm_data_im_ls[id], z_ls_extend, self.lambda_ls, an, Nt, self.mom)

            re_avg = gv.dataset.avg_data(re, bstrap=True)
            im_avg = gv.dataset.avg_data(im, bstrap=True)

            re_gv = []
            im_gv = []
            for n_conf in range(len(re)):
                re_gv.append(add_sdev(re[n_conf], re_avg))
                im_gv.append(add_sdev(im[n_conf], im_avg))

            re_lam_ls.append(re)
            im_lam_ls.append(im)
            re_gv_lam_ls.append(re_gv)
            im_gv_lam_ls.append(im_gv)

        ### continuum limit ###
        cont_re_ls, cont_im_ls = self.continuum_fit(re_gv_lam_ls, im_gv_lam_ls)
        
        ### rotate ###
        if self.rotate == True:
            quasi_re_ls, quasi_im_ls = rotate(cont_re_ls, cont_im_ls, self.lambda_ls, back=True)

        else:
            quasi_re_ls, quasi_im_ls = cont_re_ls, cont_im_ls

        return re_lam_ls, im_lam_ls, self.lambda_ls, quasi_re_ls, quasi_im_ls

        
    def paras(self, a_str):
        if a_str == 'a06':
            an = 0.0574
            Nt = 96
        elif a_str == 'a09':
            an = 0.0882
            Nt = 64
        elif a_str == 'a12':
            an = 0.1213
            Nt = 48
        return an, Nt 

    def continuum_fit(self, re_gv_lam_ls, im_gv_lam_ls):
        cont_re_ls = []
        cont_im_ls = []
        print('>>> fitting to get continuous limit for configs of mom '+str(self.mom))
        for n_conf in tqdm(range(len(re_gv_lam_ls[0]))):
            cont_re_ls.append([])
            cont_im_ls.append([])
            for idl in range(len(self.lambda_ls)):
                def fit_f0(a_ls, hyb_a):
                    def fcn(x, p): #!#
                        return p['f0'] + p['f1'] * (x) #+ p['f2'] * (x**2)
                    priors = gv.BufferDict()
                    priors['f0'] = gv.gvar(1, 10)
                    priors['f1'] = gv.gvar(1, 10)
                    priors['f2'] = gv.gvar(1, 10)

                    fit_result = lsf.nonlinear_fit(data=(a_ls, hyb_a), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

                    return fit_result.p['f0'].mean
                    # return hyb_a[0].mean #!# for continuum limit sys error 

                
                y_re = [ls[n_conf][idl] for ls in re_gv_lam_ls]
                y_im = [ls[n_conf][idl] for ls in im_gv_lam_ls]
                cont_re_ls[n_conf].append(fit_f0(self.a_ls, y_re))
                cont_im_ls[n_conf].append(fit_f0(self.a_ls, y_im))

        return cont_re_ls, cont_im_ls
        
class INV_MATCHING_COOR():
    def __init__(self, pz, lam_ls, quasi_re_ls, quasi_im_ls):
        self.pz = pz
        self.lambda_ls = lam_ls
        self.quasi_re_ls = quasi_re_ls
        self.quasi_im_ls = quasi_im_ls

    def main(self):
        start = timeit.default_timer()
        print('>>> integrating for matching of mom '+str(self.pz / mom_to_pz))

        n_conf = tqdm(np.arange(len(self.quasi_re_ls)))
        p = multiprocessing.Pool(4)#4 multiprocess to integral
        lc_ls = p.map(self.integral, n_conf)
        p.close()
        p.join()

        # lc_ls = []
        # for n_conf in tqdm(np.arange(len(self.quasi_re_ls))):
        #     lc_ls.append(self.integral(n_conf))

        end = timeit.default_timer()
        print('multi processing time:', str(end-start), 's')

        lc_re_ls = np.zeros_like(lc_ls, dtype=float)
        lc_im_ls = np.zeros_like(lc_ls, dtype=float)
        for n_conf in range(len(lc_ls)):
            for idl in range(len(lc_ls[0])):
                lc_re_ls[n_conf][idl] = lc_ls[n_conf][idl].real
                lc_im_ls[n_conf][idl] = lc_ls[n_conf][idl].imag

        return lc_re_ls, lc_im_ls

    def f_matching(self, z2, mu2):
        res = np.log( z2 * mu2 * np.exp(2*np.euler_gamma) / 4 )
        return res

    def integral(self, n_conf):
        delta_l = 0.000001 # regulater
        quasi = [1+0j] # add zero point
        for idl in range(len(self.quasi_re_ls[n_conf])):
            quasi.append( self.quasi_re_ls[n_conf][idl] + 1j * self.quasi_im_ls[n_conf][idl] )

        lam_ls_0 = np.insert(self.lambda_ls, 0, 0)
        h_tilde = interpolate.interp1d(lam_ls_0, quasi, kind='cubic')

        lc_ls = []

        def fp_2_re(lamp, lam):
            z = lam / self.pz # here z has GeV^-1
            res = alphas_cf_div_2pi / lam * (-1-self.f_matching(z**2, mu**2)) * (lamp/(lam-lamp)) * ( 
                ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                - 2 * h_tilde(lam)
            )
            return res.real

        def fp_2_im(lamp, lam):
            z = lam / self.pz # here z has GeV^-1
            res = alphas_cf_div_2pi / lam * (-1-self.f_matching(z**2, mu**2)) * (lamp/(lam-lamp)) * ( 
                ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp)
                - 2 * h_tilde(lam)
            )
            return res.imag

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

        def fp_4_re(lamp, lam):
            z = lam / self.pz # here z has GeV^-1
            res = alphas_cf_div_2pi / (1j * (lam**2)) * ( 1 - np.exp(-1j * (lam - lamp)) ) * (3-self.f_matching(z**2, mu**2)) * h_tilde(lamp)
            return res.real

        def fp_4_im(lamp, lam):
            z = lam / self.pz # here z has GeV^-1
            res = alphas_cf_div_2pi / (1j * (lam**2)) * ( 1 - np.exp(-1j * (lam - lamp)) ) * (3-self.f_matching(z**2, mu**2)) * h_tilde(lamp)
            return res.imag

        for idl in range(len(self.lambda_ls)):
            lam = self.lambda_ls[idl]
            z = lam / self.pz # here z has GeV^-1

            lc_re = self.quasi_re_ls[n_conf][idl]
            lc_im = self.quasi_im_ls[n_conf][idl]

            part1 = alphas_cf_div_2pi * 1/2 * (self.f_matching(z**2, mu**2) - 3) * h_tilde(lam)
            part2 = integrate.quad(fp_2_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_2_im, 0, lam-delta_l, args=lam)[0]
            part3 = integrate.quad(fp_3_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_3_im, 0, lam-delta_l, args=lam)[0]
            part4 = integrate.quad(fp_4_re, 0, lam, args=lam)[0] + 1j * integrate.quad(fp_4_im, 0, lam, args=lam)[0]

            lc_re = lc_re - (part1+part2+part3+part4).real
            lc_im = lc_im - (part1+part2+part3+part4).imag

            lc_ls.append(lc_re + 1j * lc_im)

        return lc_ls

class EXTRAPOLATION_FT():
    def __init__(self, meson, x_ls):
        self.meson = meson
        self.x_ls = x_ls # x_ls in the mom space
        
    def main(self, lam_ls, extend_fit_start, extend_point, da_re_ls, da_im_ls, mode):
        self.lam_ls = lam_ls # different for different pz
        self.extend_fit_start = extend_fit_start # include
        self.extend_point = extend_point
        self.lam_max = self.lam_ls[extend_point] # when lambda bigger (not include) than this max, use the extend function
        self.da_re_ls = da_re_ls
        self.da_im_ls = da_im_ls
        self.mode = mode # 'quasi' or 'lcda'

        re_avg = gv.dataset.avg_data(self.da_re_ls, bstrap=True)
        im_avg = gv.dataset.avg_data(self.da_im_ls, bstrap=True)
        
        da_ext_ls = []
        pn_conf_ls = []
        print('>>> fitting to extend DA hybrid to lambda max of '+self.mode)
        for n_conf in tqdm(range(len(self.da_re_ls))):
            re = self.da_re_ls[n_conf]
            im = self.da_im_ls[n_conf]
            re_gv = add_sdev(re, re_avg)
            im_gv = add_sdev(im, im_avg)
            lam_ls_ex, da_ext, n = self.extrapolation(re, im, re_gv, im_gv)

            da_ext_ls.append(da_ext)
            pn_conf_ls.append(n)

        if self.meson == 'pion':
            avg_n = gv.dataset.avg_data(pn_conf_ls, bstrap=True)
            print('Average n over all configs is '+str(avg_n.mean))

        ### FT ###
        ##########
        if ((lam_ls_ex[1] - lam_ls_ex[0]) - (self.lam_ls[1] - self.lam_ls[0])) > 0.0001:
            print('lambda delta error')
            
        lam_delta = lam_ls_ex[1] - lam_ls_ex[0]
        da_mom_ls = []
        print('>>> FT to momentum space of '+self.mode)
        for n_conf in tqdm(range(len(da_ext_ls))):
            da_mom_ls.append([])
            for x in self.x_ls:
                ft_res = sum_ft(lam_ls_ex, da_ext_ls[n_conf], lam_delta, x)
                da_mom_ls[n_conf].append(ft_res.real)

        return lam_ls_ex, da_ext_ls, da_mom_ls

    def extrapolation(self, re, im, re_gv, im_gv):
        if self.mode == 'quasi':
            extend_length = extend_length_quasi

        elif self.mode == 'lcda':
            extend_length = extend_length_lc

        lam_delta = self.lam_ls[1] - self.lam_ls[0]
        lam_max = self.lam_max # when lambda bigger than this max, use the extend function

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
            val = {}
            if self.meson == 'pion':
                val['re'] = ( p['c1']/(x['re']**p['n']) * np.cos(np.pi/2 * p['n']) + p['c2']/(x['re']**p['n']) * np.cos(x['re'] - np.pi/2 * p['n']) )  #e# ###

                val['im'] = - ( p['c1']/(x['im']**p['n']) * np.sin(np.pi/2 * p['n']) + p['c2']/(x['im']**p['n']) * np.sin(x['im'] - np.pi/2 * p['n']) )

                if self.mode == 'quasi':
                    val['re'] = val['re'] * np.exp(-x['re'] / 150)

                    val['im'] = val['im'] * np.exp(-x['im'] / 150) 

            elif self.meson == 'kaon':
                val['re'] = ( p['c1']/(x['re']**p['n1']) * np.cos(np.pi/2 * p['n1']) + p['c2']/(x['re']**p['n2']) * np.cos(x['re'] - np.pi/2 * p['n2']) )

                val['im'] = - ( p['c1']/(x['im']**p['n1']) * np.sin(np.pi/2 * p['n1']) + p['c2']/(x['im']**p['n2']) * np.sin(x['im'] - np.pi/2 * p['n2']) )

                if self.mode == 'quasi':
                    val['re'] = val['re'] * np.exp(-x['re'] / 150)

                    val['im'] = val['im'] * np.exp(-x['im'] / 150) 

            return val

        lam_dic = {}
        lam_dic['re'] = np.array(self.lam_ls[self.extend_fit_start:])
        lam_dic['im'] = np.array(self.lam_ls[self.extend_fit_start:])

        da_dic = {}
        da_dic['re'] = re_gv[self.extend_fit_start:]
        da_dic['im'] = im_gv[self.extend_fit_start:]

        fit_result = lsf.nonlinear_fit(data=(lam_dic, da_dic), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')
    
        da_ext = [] # after extend, combine real and imag to do ft
        lam_ls_ex = [] # lambda list after extend
    
        ### piece together to form the renorm_da in the whole coordinate space ###
        ###########################################
        for lam in np.arange(-lam_delta*extend_length, -lam_max-0.001, lam_delta): # use extend func, complex conjugate
            lam_ls_ex.append(lam)
            lam_dic = {'re':-lam, 'im':-lam}
            da_ext.append( complex(fcn(lam_dic, fit_result.p)['re'].mean, -fcn(lam_dic, fit_result.p)['im'].mean) )

        for idx in range(len(self.lam_ls[:self.extend_point])+1): # use data list, complex conjugate
            lam_ls_ex.append(-self.lam_max + idx*lam_delta)
            da_ext.append( complex(re[self.extend_point-idx], -im[self.extend_point-idx]) )

        lam_ls_ex.append(0) # zero point
        da_ext.append(1+0j)

        for idx in range(len(self.lam_ls[:self.extend_point])+1): # use data list
            lam_ls_ex.append(self.lam_ls[idx])
            da_ext.append( complex(re[idx], im[idx]) )

        for lam in np.arange(lam_max+lam_delta, lam_delta*extend_length+0.001, lam_delta): # use extend func
            lam_ls_ex.append(lam)
            lam_dic = {'re':lam, 'im':lam}
            da_ext.append( complex(fcn(lam_dic, fit_result.p)['re'].mean, fcn(lam_dic, fit_result.p)['im'].mean) )

        if self.meson == 'pion':
            return lam_ls_ex, da_ext, fit_result.p['n'].mean

        elif self.meson == 'kaon':
            return lam_ls_ex, da_ext, 0

class INV_MATCHING_MOM():
    def __init__(self, x_ls, x_ls_mat, y_ls_mat):
        self.x_ls = x_ls # for original quasi
        self.x_ls_mat = x_ls_mat # for quasi in matching
        self.y_ls_mat = y_ls_mat # for light-cone in matching

    def main(self, pz, quasi_mom_ls):
        self.pz = pz
        kernel = self.matching_kernel()

        lc_mom_ls = []
        print('>>> matching in the momentum space of mom '+str(self.pz / mom_to_pz))
        for n_conf in tqdm(range(len(quasi_mom_ls))):
            quasi = np.array(quasi_mom_ls[n_conf])
            quasi = interp_1d(self.x_ls, quasi, self.x_ls_mat, method='cubic') # interpolation to match a big matrix
            dot = np.dot( kernel, quasi ) ###
            lc = interp_1d(self.y_ls_mat, dot, self.x_ls, method='cubic') # interpolation back to previous x_ls
            lc_mom_ls.append( lc ) 

        return lc_mom_ls

    def matching_kernel(self):
        x_ls = self.x_ls_mat # for quasi
        y_ls = self.y_ls_mat # for light-cone

        def H1(x, y):
            return (1+x-y)/(y-x)*(1-x)/(1-y)*np.log((y-x)/(1-x)) + (1+y-x)/(y-x)*x/y*np.log((y-x)/(-x))
        
        def H2(x, y):
            return (1+y-x)/(y-x)*x/y*np.log(4*x*(y-x)*self.pz**2/mu_f**2) + (1+x-y)/(y-x)*((1-x)/(1-y)*np.log((y-x)/(1-x))-x/y)

        ### CB_matrix ###
        #################
        CB_matrix = np.zeros([len(x_ls), len(y_ls)])
        for idx1 in range(len(x_ls)):
            for idx2 in range(len(y_ls)):
                x = x_ls[idx1]
                y = y_ls[idx2] #!#
                if abs(x-y) > 0.0001:
                    if x < 0 and y > 0 and y < 1:
                        CB_matrix[idx1][idx2] = H1(x, y)
                    elif x > 0 and y > x and y < 1:
                        CB_matrix[idx1][idx2] = H2(x, y)
                    elif y > 0 and y < x and x < 1:
                        CB_matrix[idx1][idx2] = H2(1-x, 1-y)
                    elif y > 0 and y < 1 and x > 1:
                        CB_matrix[idx1][idx2] = H1(1-x, 1-y)

        CB_matrix = CB_matrix * alphas_cf_div_2pi

        for idx in range(len(x_ls)): # diagnoal input
            if CB_matrix[idx][idx] != 0:
                print('CB matrix diagnoal error')
            CB_matrix[idx][idx] = -np.sum([CB_matrix[i][idx] for i in range(len(x_ls))])

        ### extra term related to the modified hybrid method ###
        #################################################
        extra_term = np.zeros([len(x_ls), len(y_ls)])
        for idx1 in range(len(x_ls)):
            for idx2 in range(len(y_ls)):
                x = x_ls[idx1]
                y = y_ls[idx2]
                if y > 0 and y < 1 and abs(x-y) > 0.0001:
                    extra_term[idx1][idx2] = 3/2 * (1 / abs(x-y))

        for idx in range(len(x_ls)): # diagnoal input
            if extra_term[idx][idx] != 0:
                print('extra term matrix diagnoal error')
            extra_term[idx][idx] = -np.sum([extra_term[i][idx] for i in range(len(x_ls))])

        extra_term = extra_term * alphas_cf_div_2pi

        ### delta(x-y) ###
        ##################
        identity = np.zeros([len(x_ls), len(y_ls)])

        for idx in range(len(x_ls)):
            identity[idx][idx] = 1

        C_matrix = (CB_matrix + extra_term) * (x_ls[1]-x_ls[0]) + identity # multiply by da to represent integral

        C_matrix_inverse = np.linalg.inv(C_matrix)

        return C_matrix_inverse

def large_mom_limit(x_ls, mom_da_ls, mom_ls):
    large_mom_da = []
    for idx in range(len(x_ls)):
        def fcn(x, p):
            return p['psi'] + p['c2']/(x**2)

        priors = gv.BufferDict()
        priors['psi'] = gv.gvar(1, 10)
        priors['c2'] = gv.gvar(1, 10)

        pz_ls = np.array(mom_ls) * mom_to_pz
        lcda_ls = [ mom_da[idx] for mom_da in mom_da_ls]

        fit_result = lsf.nonlinear_fit(data=(pz_ls, lcda_ls), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

        large_mom_da.append(fit_result.p['psi'])

    return large_mom_da






# %%
# check with yushan
# meson = 'pion'

# lam_ls = gv.load('lam_ls')
# lc_re_avg = gv.load('lc_re')
# lc_im_avg = gv.load('lc_im')

# lc_re = [val.mean for val in lc_re_avg]
# lc_im = [val.mean for val in lc_im_avg]

# extra = EXTRAPOLATION_FT(meson, lam_ls, -16, -10, lc_re, lc_im, mode='lcda')
# lam_ls_ex, lc_ext, n = extra.extrapolation(lc_re, lc_im, lc_re_avg, lc_im_avg)

# print(np.shape(lc_ext))
# print(np.shape(lam_ls_ex))

# lc_ext_re = []
# lc_ext_im = []

# for id in range(len(lc_ext)):
#     lc_ext_re.append(lc_ext[id].real)
#     lc_ext_im.append(lc_ext[id].imag)
    
# x_ls = np.arange(-2-0.01, 3.02, 0.01) # x after ft, for quasi before matching
# lc_re_ft = []
# lc_im_ft = []
# for x in x_ls:
#     ft_res = sum_ft(lam_ls_ex, lc_ext, lam_ls_ex[1]-lam_ls_ex[0], x)
#     lc_re_ft.append(ft_res.real)
#     lc_im_ft.append(ft_res.imag)

# plt.plot(lam_ls_ex, lc_ext_re)
# plt.show()

# plt.plot(lam_ls_ex, lc_ext_im)
# plt.show()

# plt.plot(x_ls, lc_re_ft)
# plt.xlim([-0.5, 1.5])
# plt.show()

# %%
# f = open('./lc_ext.txt', 'w')
# line = []
# line.append('lambda')
# line.append('\t')
# line.append('lc_Re')
# line.append('\t')
# line.append('lc_Im')
# line.append('\n')
# f.writelines(line)

# for idz in range(len(lam_ls_ex)):
#     line = []
#     line.append(str(lam_ls_ex[idz]))
#     line.append('\t')
#     line.append(str(lc_ext_re[idz]))
#     line.append('\t')
#     line.append(str(lc_ext_im[idz]))
#     line.append('\n')
#     f.writelines(line)
# f.close


# f = open('./lc_ft.txt', 'w')
# line = []
# line.append('x')
# line.append('\t')
# line.append('lc_ft')
# line.append('\n')
# f.writelines(line)

# for idz in range(len(x_ls)):
#     line = []
#     line.append(str(x_ls[idz]))
#     line.append('\t')
#     line.append(str(lc_re_ft[idz]))
#     line.append('\n')
#     f.writelines(line)
# f.close
# %%
