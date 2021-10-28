# %%
from gvar import chi2
from numpy.lib.function_base import average
from head import *

# %%
class MESON_DA_HYB():
    def __init__(self, meson, mom, x_ls, x_ls_matching, y_ls_matching, extend_point, extend_fit_start, t_dic, gs_extract, fit_1, fit_2, rotate, lcda_fit, constant_fit):
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
        self.fit_2 = fit_2 # True/False, whether fit to extend hybrid data in the coordinate space
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

        ### renorm -> continuous limit -> extrapolation -> FT -> quasi ###
        ##########################################
        a_hyb_ro_re_avg, a_hyb_ro_im_avg, hyb_quasi = self.da_hyb(da_a06_conf_ls, da_a09_conf_ls, da_a12_conf_ls, zR_dic) # hyb_quasi.shape = (N_conf, len(x_ls))

        kernel = self.matching_kernel()

        ### quasi -> inv_matching -> light-cone ###
        ###########################################
        quasi_da, y_ls, lic_da = self.inverse_matching(hyb_quasi, kernel)

        return kernel, quasi_da, y_ls, lic_da    
    
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

        ### fit z to get DA hybrid over the whole z-axis, then do ft to momentum space ###
        ###########################################
        if self.fit_2 == True:
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

                def fcn(x, p, no_e=False):
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
                for x in x_ls:
                    hyb_conf_re_ft.append(sum_ft(lam_ls_ex, hyb_conf, lambda_delta, x).real)
                    hyb_conf_im_ft.append(sum_ft(lam_ls_ex, hyb_conf, lambda_delta, x).imag)

                hyb_complex.append(hyb_conf)
                hyb_re_ft.append(hyb_conf_re_ft)
                hyb_im_ft.append(hyb_conf_im_ft)

            if self.meson == 'pion':  
                avg_n = gv.dataset.avg_data(pn_conf_ls, bstrap=True)
                print('Average n over all configs is '+str(avg_n.mean))

            gv.dump(lam_ls_ex, self.meson+'/mom='+str(self.mom)+'/lam_ls_ex') 
            gv.dump(hyb_complex, self.meson+'/mom='+str(self.mom)+'/hyb_complex') # before ft
            gv.dump(hyb_re_ft, self.meson+'/mom='+str(self.mom)+'/hyb_re_ft') # shape = (N_re, len(x_ls))
            gv.dump(hyb_im_ft, self.meson+'/mom='+str(self.mom)+'/hyb_im_ft')

        lam_ls_ex = gv.load(self.meson+'/mom='+str(self.mom)+'/lam_ls_ex')
        hyb_complex = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_complex')
        hyb_re = []
        hyb_im = []
        for idx in range(len(hyb_complex)):
            hyb_re.append( [val.real for val in hyb_complex[idx]] )
            hyb_im.append( [val.imag for val in hyb_complex[idx]] )
        hyb_re_ft = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_re_ft')
        hyb_im_ft = gv.load(self.meson+'/mom='+str(self.mom)+'/hyb_im_ft')

        return a_hyb_ro_re_avg, a_hyb_ro_im_avg, hyb_re_ft

    def matching_kernel(self):
        x_ls = self.x_ls_matching # for quasi
        y_ls = self.y_ls_matching # for light-cone

        def H1(x, y):
            return (1+x-y)/(y-x)*(1-x)/(1-y)*np.log((y-x)/(1-x)) + (1+y-x)/(y-x)*x/y*np.log((y-x)/(-x))
        
        def H2(x, y):
            return (1+y-x)/(y-x)*x/y*np.log(4*x*(y-x)*self.pz**2/mu**2) + (1+x-y)/(y-x)*((1-x)/(1-y)*np.log((y-x)/(1-x))-x/y)

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

        for idx in range(len(x_ls)): # diagonal input
            if CB_matrix[idx][idx] != 0:
                print('CB matrix diagonal error')
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

        for idx in range(len(x_ls)): # diagonal input
            if extra_term[idx][idx] != 0:
                print('extra term matrix diagonal error')
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

    def inverse_matching(self, hyb_quasi, kernel):
        x_ls = self.x_ls
        x_ls_matching = self.x_ls_matching
        y_ls_matching = self.y_ls_matching

        hyb_lic = []

        for ls in hyb_quasi: # each config
            quasi = np.array(ls)
            quasi = interp_1d(x_ls, quasi, x_ls_matching, method='cubic') # interpolation to match a big matrix
            dot = np.dot( kernel, quasi ) ###
            hyb_lic.append( dot ) 

        ## gaussian filter to eliminate some discontinuity at x=0 and x=1
        #hyb_lic = gaussian_filter(hyb_lic, sigma=1)

        hyb_quasi_avg = gv.dataset.avg_data(hyb_quasi, bstrap=True) # before inv matching, so corresponding to x_ls

        hyb_lic_avg = gv.dataset.avg_data(hyb_lic, bstrap=True)

        ### fit the region near the end points ###
        ##########################################
        if self.lcda_fit == True:
            lic_fit = []
            for ls in hyb_lic:
                y_1 = [] # [0, 0.05] and [0.95, 1] use fit results
                y_2 = [] # [0.05, 0.2] and [0.8, 0.95] fit region
                y_3 = [] # [0.2, 0.8] keep original 

                lcda_1 = []
                lcda_2 = []
                lcda_3 = []
                
                for idy in range(len(y_ls_matching)):
                    y = y_ls_matching[idy]
                    if 0 <= y <= fit_y1 or fit_y4 <= y <= 1: # use fit
                        y_1.append(y)
                    elif fit_y1 <= y <= fit_y2 or fit_y3 <= y <= fit_y4: # fit region
                        y_2.append(y)
                        lcda_2.append(gv.gvar(ls[idy], hyb_lic_avg[idy].sdev))
                    elif fit_y2 <= y <= fit_y3: # keep originally
                        y_3.append(y)
                        lcda_3.append(ls[idy])

                if self.meson == 'pion':
                    def fcn(x, p):
                        return p['c'] * ( x**p['n'] ) * ( (1-x)**p['n'] )

                elif self.meson == 'kaon':
                    def fcn(x, p):
                        return p['c'] * ( x**p['n1'] ) * ( (1-x)**p['n2'] )

                priors = gv.BufferDict()
                priors['c'] = gv.gvar(0, 10)
                priors['n'] = gv.gvar(0, 10)
                priors['n1'] = gv.gvar(0, 10)
                priors['n2'] = gv.gvar(0, 10)

                fit_result = lsf.nonlinear_fit(data=(np.array(y_2), lcda_2), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

                y_ls = np.hstack( (np.array(y_1), np.array(y_2), np.array(y_3)) )
                lcda_1 = [val.mean for val in fcn(np.array(y_1), fit_result.p)]
                lcda_2 = [val.mean for val in lcda_2]
                lic_da = np.hstack( (lcda_1, lcda_2, lcda_3) )

                z = zip(y_ls, lic_da)
                z = sorted(z)
                y_ls, lic_da = zip(*z)

                lic_fit.append(lic_da)
                ##

            lic_fit_avg = gv.dataset.avg_data(lic_fit, bstrap=True)

        if self.lcda_fit == False:
            ### normalization check ###
            ###########################
            print('Light-cone integral within [0, 1]: ')
            normal = norm_check(y_ls_matching, hyb_lic_avg) 
            print(normal)

            ### a2 check ###
            ################
            print('Pz='+str(self.pz)+', a2: ')
            a2 = calc_an(y_ls_matching, hyb_lic_avg, 2)
            print(a2)

            return hyb_quasi_avg, y_ls_matching, hyb_lic_avg

        elif self.lcda_fit == True:
            ### normalization check ###
            ###########################
            print('Light-cone integral within [0, 1]: ')
            normal = norm_check(y_ls, lic_fit_avg) 
            print(normal)

            ### a2 check ###
            ################
            print('Pz='+str(self.pz)+', a2: ')
            a2 = calc_an(y_ls, lic_fit_avg, 2)
            print(a2)

            return hyb_quasi_avg, y_ls, lic_fit_avg
