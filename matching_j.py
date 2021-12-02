#coding=utf-8
# %%
import numpy as np
from numpy.core.defchararray import count
import math as m
import gvar as gv
import lsqfit
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
from scipy import integrate
from tqdm import tqdm
import multiprocessing
import timeit


def main():
    #MT = MATCHING_PROCESS()
    #quasi_re_ls, quasi_im_ls, lambda_ls = MT.read(plot=True)
    #MT.matching(quasi_re_ls, quasi_im_ls, lambda_ls, plot=True)
    LFT = LIGHTCONE_FT()
    #LFT.plot_read(plot=True)
    #LFT.extrapolation(st_point=7, plrg=5, max_point=500, plot=True)
    LFT.phi_FT(max_point=500, plot=True)


class MATCHING_PROCESS(object):
    def __init__(self):
        self.pz = 2.15
        self.mu = 2 # GeV, for renormalization
        self.cf=4/3
        self.nf=3
        self.b0=11-2/3*self.nf
        self.lms = 0.24451721864451428  
        self.alphas=2 * np.pi/( self.b0 * np.log(self.mu/self.lms) )
        self.alphas_cf2pi = self.alphas * self.cf / (2 * np.pi)
        self.quasi_re_ls=[]
        self.quasi_im_ls=[]
        self.lambda_ls=[]

    def read(self, plot=False):
        quasi_re = np.array(gv.load('quasi_ro_re_ls'))
        quasi_im = np.array(gv.load('quasi_ro_im_ls'))
        lam_ls = np.array(gv.load('lam_ls'))

        q_re_avg = gv.dataset.avg_data(quasi_re, bstrap=True)
        q_im_avg = gv.dataset.avg_data(quasi_im, bstrap=True)
        if plot==True:  plt.fill_between(lam_ls, [(val.mean + val.sdev) for val in q_re_avg], [(val.mean - val.sdev) for val in q_re_avg], color='pink', alpha=0.8)
        #plt.show()
        def plot_phase(lam_ls, quasi_re, quasi_im):
            n_conf = quasi_re.shape[0]
            combine_re = [];    combine_im = []
            for iz in np.arange(0, quasi_re.shape[1]):
                phase = (iz+1)*lam_ls[0]/2
                combine_re = np.append(combine_re, quasi_re[:,iz] * np.cos(phase) - quasi_im[:,iz] * np.sin(phase))
                combine_im = np.append(combine_im, quasi_im[:,iz] * np.cos(phase) + quasi_re[:,iz] * np.sin(phase))
            combine_re = np.swapaxes(combine_re.reshape(quasi_re.shape[1], n_conf), 1, 0)
            combine_im = np.swapaxes(combine_im.reshape(quasi_re.shape[1], n_conf), 1, 0)
            cv_re = np.mean(combine_re, axis=0);    err_re  = np.std(combine_re, axis=0)
            cv_im = np.mean(combine_im, axis=0);    err_im  = np.std(combine_im, axis=0)
            if plot==True:    plt.errorbar(lam_ls, cv_re, err_re, capsize=4.0,fmt='D', ms=3.0, color=col[0], label = 're')
            if plot==True:    plt.errorbar(lam_ls, cv_im, err_im, capsize=4.0,fmt='D', ms=3.0, color=col[1], label = 'im')
            plt.show()
        #plot_phase(lam_ls, quasi_re, quasi_im)

        return quasi_re, quasi_im, lam_ls

    def matching(self, quasi_re_ls, quasi_im_ls, lambda_ls, plot=False):
        self.quasi_re_ls = quasi_re_ls
        self.quasi_im_ls = quasi_im_ls
        self.lambda_ls = lambda_ls

        ### multi process ###
        start = timeit.default_timer()
        #items_cfg = [x for x in range(len(quasi_re_ls))]
        items_cfg = tqdm(np.arange(len(quasi_re_ls)))
        p = multiprocessing.Pool(6)#4 multiprocess to integral
        lc = p.map(self.integral_fun, items_cfg)
        p.close()
        p.join()
        end = timeit.default_timer()
        print('multi processing time:', str(end-start), 's')

        lc = np.array(lc)
        lc_re = np.empty([lc.shape[0], lc.shape[1]]);   lc_im = np.empty([lc.shape[0], lc.shape[1]])
        for icg in tqdm(range(lc.shape[0])):
            for iz in range(lc.shape[1]):
                lc_re[icg, iz] = lc[icg, iz].real
                lc_im[icg, iz] = lc[icg, iz].imag
        gv.dump(lc_re, 'lc_re')
        gv.dump(lc_im, 'lc_im')
        lc_re_avg = gv.dataset.avg_data(lc_re, bstrap=True)
        lc_im_avg = gv.dataset.avg_data(lc_im, bstrap=True)
        if plot==True:  plt.fill_between(self.lambda_ls, [(val.mean + val.sdev) for val in lc_re_avg], [(val.mean - val.sdev) for val in lc_re_avg], color='green', alpha=0.5)
        print(lc_re_avg)
        plt.show()
        
    def f_matching(self, z2, mu2):
        res = np.log( z2 * mu2 * np.exp(2*np.euler_gamma) / 4 )
        return res

    def integral_fun(self, n_conf):
        delta_l = 0.000001 # regulater
        quasi = [1+0j] # add the point at lambda=0
        for idl in range(len(self.quasi_re_ls[n_conf])):
            quasi.append( self.quasi_re_ls[n_conf][idl] + 1j * self.quasi_im_ls[n_conf][idl] )
        lam_ls_0 = np.insert(self.lambda_ls, 0, 0) # add the point at lambda=0
        h_tilde = interpolate.interp1d(lam_ls_0, quasi, kind='cubic') # analytic function of quasi

        lc_re = [];     lc_im = []
        for idl in range(len(self.lambda_ls)):
            lam = self.lambda_ls[idl]
            z = lam / self.pz # here z has GeV^-1
            
            lc_re_idl = self.quasi_re_ls[n_conf][idl] #
            lc_im_idl = self.quasi_im_ls[n_conf][idl]

            part1 = self.alphas_cf2pi * 1/2 * (self.f_matching(z**2, self.mu**2) - 3) * h_tilde(lam)

            def fp_2_re(lamp, lam):
                z = lam / self.pz # here z has GeV^-1
                res = self.alphas_cf2pi / lam * (-1-self.f_matching(z**2, self.mu**2)) * (lamp/(lam-lamp)) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp) - 2 * h_tilde(lam) )
                return res.real
            def fp_2_im(lamp, lam):
                z = lam / self.pz # here z has GeV^-1
                res = self.alphas_cf2pi / lam * (-1-self.f_matching(z**2, self.mu**2)) * (lamp/(lam-lamp)) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp) - 2 * h_tilde(lam) )
                return res.imag
            part2 = integrate.quad(fp_2_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_2_im, 0, lam-delta_l, args=lam)[0] # [0] is result, [1] is error estimation

            def fp_3_re(lamp, lam):
                res = self.alphas_cf2pi / lam * np.log(1 - lamp/lam) / (1 - lamp/lam) * (-2) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp) - 2 * h_tilde(lam) )
                return res.real
            def fp_3_im(lamp, lam):
                res = self.alphas_cf2pi / lam * np.log(1 - lamp/lam) / (1 - lamp/lam) * (-2) * ( 
                    ( 1 + np.exp(-1j * (lam - lamp)) ) * h_tilde(lamp) - 2 * h_tilde(lam) )
                return res.imag
            part3 = integrate.quad(fp_3_re, 0, lam-delta_l, args=lam)[0] + 1j * integrate.quad(fp_3_im, 0, lam-delta_l, args=lam)[0]

            def fp_4_re(lamp, lam):
                z = lam / self.pz # here z has GeV^-1
                res = self.alphas_cf2pi * ( 1 - np.exp(-1j * (lam - lamp)) ) / (1j * (lam**2)) * (3-self.f_matching(z**2, self.mu**2)) * h_tilde(lamp)
                return res.real
            def fp_4_im(lamp, lam):
                z = lam / self.pz # here z has GeV^-1
                res = self.alphas_cf2pi * ( 1 - np.exp(-1j * (lam - lamp)) ) / (1j * (lam**2))* (3-self.f_matching(z**2, self.mu**2)) * h_tilde(lamp)
                return res.imag
            part4 = integrate.quad(fp_4_re, 0, lam, args=lam)[0] + 1j * integrate.quad(fp_4_im, 0, lam, args=lam)[0]

            lc_re_idl = lc_re_idl - (part1+part2+part3+part4).real
            lc_im_idl = lc_im_idl - (part1+part2+part3+part4).imag

            lc_re = np.append(lc_re, lc_re_idl)
            lc_im = np.append(lc_im, lc_im_idl)
        lc = lc_re + 1j*lc_im
        return lc

class LIGHTCONE_FT(object):
    def plot_read(self, plot=False):
        quasi_re = np.array(gv.load('quasi_ro_re_ls'))
        quasi_im = np.array(gv.load('quasi_ro_im_ls'))
        lam_ls = np.array(gv.load('lam_ls'))

        lc_re = np.array(gv.load('lc_re'))
        lc_im = np.array(gv.load('lc_im'))
        print(lam_ls)

        q_re_avg = gv.dataset.avg_data(quasi_re, bstrap=True)
        q_im_avg = gv.dataset.avg_data(quasi_im, bstrap=True)
        
        l_re_avg = gv.dataset.avg_data(lc_re, bstrap=True)
        l_im_avg = gv.dataset.avg_data(lc_im, bstrap=True)

        if plot==True:
            plt.fill_between(lam_ls, [(val.mean + val.sdev) for val in q_re_avg], [(val.mean - val.sdev) for val in q_re_avg], color='pink', alpha=0.8, label='quasi')
            plt.fill_between(lam_ls, [(val.mean + val.sdev) for val in l_re_avg], [(val.mean - val.sdev) for val in l_re_avg], color='green', alpha=0.5, label='LCDA')
            plt.axhline(0,linestyle='--',color='black') 
            plt.xlabel('lambda')
            plt.legend();   plt.show()

        return lam_ls, quasi_re, quasi_im, lc_re, lc_im

    def extrapolation(self, st_point=7, plrg=5, max_point=200, plot=False):
        lam_ls, quasi_re, quasi_im, lc_re, lc_im = self.plot_read()

        def extra_function(x, paras):#extrapolation function
            ans = {}# ratio 
            ans['re'] = (paras['c1'] /(x['lam'])**paras['d1'] * np.cos(np.pi/2 * paras['d1']) + paras['c2'] /(x['lam'])**paras['d1'] * np.cos(x['lam'] - np.pi/2 * paras['d1']))
            ans['im'] = (- paras['c1'] /(x['lam'])**paras['d1'] * np.sin(np.pi/2 * paras['d1']) - paras['c2'] /(x['lam'])**paras['d1'] * np.sin(x['lam'] - np.pi/2 * paras['d1']))
            return ans

        ############# extraplolation for LCDA ###############
        n_conf = lc_re.shape[0]
        ext_lc_re = np.empty([n_conf, max_point]);    ext_lc_re[:, 0:(st_point+plrg)] = lc_re[:, 0:(st_point+plrg)]#data use orginal data
        ext_lc_im = np.empty([n_conf, max_point]);    ext_lc_im[:, 0:(st_point+plrg)] = lc_im[:, 0:(st_point+plrg)]
        err_lc_r = [val.sdev for val in gv.dataset.avg_data(lc_re[:, st_point:], bstrap=True)]#error of fit range
        err_lc_i = [val.sdev for val in gv.dataset.avg_data(lc_im[:, st_point:], bstrap=True)]

        x_fit = {}; y_fit = {}
        x_fit['lam'] = lam_ls
        x_fit['lam'] = x_fit['lam'][st_point:]
        print('Extraplolate:')
        ############# extraplolation for LCDA ###############
        for i in tqdm(range(0, n_conf)):
            y_fit['re'] = gv.gvar(lc_re[i, st_point:], err_lc_r)
            y_fit['im'] = gv.gvar(lc_im[i, st_point:], err_lc_i)
            prior = gv.gvar(dict(c1='1.0(100)', c2='1.0(100)', d1='1.0(100)', l0='1.0(50)'))
            fit = lsqfit.nonlinear_fit(data=(x_fit, y_fit), svdcut=1e-8, prior=prior, fcn=extra_function)

            ###########  constract data after extraplolation at > st_point+plrg #############
            re_x = {};    re_x['lam'] = np.arange(st_point+plrg+1, max_point+1)*lam_ls[0]
            re_data = extra_function(re_x, fit.p)
            re_real = [];   re_imag =[]
            for iz in np.arange(0, re_x['lam'].shape[0]):
                re_real = np.append(re_real, re_data['re'][iz].mean)
                re_imag = np.append(re_imag, re_data['im'][iz].mean)
            ext_lc_re[i, st_point+plrg:] = re_real
            ext_lc_im[i, st_point+plrg:] = re_imag
        
        ############# extraplolation for quasi ###############
        ext_qs_re = np.empty([n_conf, max_point]);    ext_qs_re[:, 0:(st_point+plrg)] = quasi_re[:, 0:(st_point+plrg)]#data use orginal data
        ext_qs_im = np.empty([n_conf, max_point]);    ext_qs_im[:, 0:(st_point+plrg)] = quasi_im[:, 0:(st_point+plrg)]
        err_qs_r = [val.sdev for val in gv.dataset.avg_data(quasi_re[:, st_point:], bstrap=True)]#error of fit range
        err_qs_i = [val.sdev for val in gv.dataset.avg_data(quasi_im[:, st_point:], bstrap=True)]
        for i in tqdm(range(0, n_conf)):
            y_fit['re'] = gv.gvar(quasi_re[i, st_point:], err_qs_r)
            y_fit['im'] = gv.gvar(quasi_im[i, st_point:], err_qs_i)
            prior = gv.gvar(dict(c1='1.0(100)', c2='1.0(100)', d1='1.0(100)', l0='1.0(50)'))
            fit = lsqfit.nonlinear_fit(data=(x_fit, y_fit), svdcut=1e-8, prior=prior, fcn=extra_function)
            ###########  constract data after extraplolation at > st_point+plrg #############
            re_x = {};    re_x['lam'] = np.arange(st_point+plrg+1, max_point+1)*lam_ls[0]
            re_data = extra_function(re_x, fit.p)
            re_real = [];   re_imag =[]
            for iz in np.arange(0, re_x['lam'].shape[0]):
                re_real = np.append(re_real, re_data['re'][iz].mean)
                re_imag = np.append(re_imag, re_data['im'][iz].mean)
            ext_qs_re[i, st_point+plrg:] = re_real
            ext_qs_im[i, st_point+plrg:] = re_imag


        ##########plot data after extrapolation################
        x = np.arange(1, max_point+1)**lam_ls[0]
        cv_r = np.mean(ext_lc_re, axis=0);    err_r = np.std(ext_lc_re, axis=0)
        cv_i = np.mean(ext_lc_im, axis=0);    err_i = np.std(ext_lc_im, axis=0)
        if plot==True:   plt.fill_between(x, cv_r-err_r, cv_r+err_r, color=col[3], label='ext_lc_re')
        if plot==True:   plt.fill_between(x, cv_i-err_i, cv_i+err_i, color=col[4], label='ext_lc_im')
        if plot==True:   
            plt.axhline(0,linestyle='--',color='black') 
            plt.xlabel('lambda')
            plt.legend();   plt.show()
        
        return ext_qs_re, ext_qs_im, ext_lc_re, ext_lc_im

    def phi_FT(self, max_point=500, plot=True):
        lam_ls = np.array(gv.load('lam_ls'))
        quasi_re, quasi_im, lc_re, lc_im = self.extrapolation(st_point=9, plrg=5, max_point=max_point, plot=False)
        n_conf = quasi_re.shape[0]

        dz = lam_ls[0]
        z = np.arange(-(max_point), max_point+1) * dz
        x = np.arange(-2, 3, 0.01)
        print('FT:')
        ############# FT for LCDA #############
        lcda = np.empty([n_conf, x.shape[0]])######data shape after FT
        for i in tqdm(range(0, n_conf)):
            lc_i_re = lc_re[i, :];    lc_i_im = lc_im[i, :]
            all_re = np.concatenate((np.append(lc_i_re[::-1],1), lc_i_re[0:]), axis=0)
            all_im = -np.concatenate((np.append(-lc_i_im[::-1],0), lc_i_im[0:]), axis=0)
            combine = all_re + all_im*1j
            lcda_i = []
            for ix in x:####### FT function ######
                wave_x = dz * np.sum(1/(2*np.pi) * np.exp(-z*(ix)*1j) * combine)
                lcda_i = np.append(lcda_i, np.real(wave_x))
            lcda[i, :] = lcda_i
        ############# FT for Quasi #############
        qda = np.empty([n_conf, x.shape[0]])######data shape after FT
        for i in tqdm(range(0, n_conf)):
            quasi_i_re = quasi_re[i, :];    quasi_i_im = quasi_im[i, :]
            all_re = np.concatenate((np.append(quasi_i_re[::-1],1), quasi_i_re[0:]), axis=0)
            all_im = -np.concatenate((np.append(-quasi_i_im[::-1],0), quasi_i_im[0:]), axis=0)
            combine = all_re + all_im*1j
            qda_i = []
            for ix in x:####### FT function ######
                wave_x = dz * np.sum(1/(2*np.pi) * np.exp(-z*(ix)*1j) * combine)
                qda_i = np.append(qda_i, np.real(wave_x))
            qda[i, :] = qda_i

        lcda_cv = np.array([val.mean for val in gv.dataset.avg_data(lcda, bstrap=True)])
        lcda_err = np.array([val.sdev for val in gv.dataset.avg_data(lcda, bstrap=True)])
        qda_cv = np.array([val.mean for val in gv.dataset.avg_data(qda, bstrap=True)])
        qda_err = np.array([val.sdev for val in gv.dataset.avg_data(qda, bstrap=True)])

        plt.fill_between(x, qda_cv-qda_err, qda_cv+qda_err, color=col[0], alpha=0.5, label='Quasi')
        plt.fill_between(x, lcda_cv-lcda_err, lcda_cv+lcda_err, color=col[1], alpha=0.5, label='LCDA')
        plt.axhline(0,linestyle='--',color='black') 
        plt.axvline(0,linestyle='--',color='green') 
        plt.axvline(1,linestyle='--',color='green') 
        plt.xlabel('x')
        plt.xlim(-0.5,1.5)
        plt.legend();   plt.show()

        




col = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue',\
'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue',\
'red','green']
fmts = ['D', '^', 'o', 's'] 
colors = ['blue', 'red', 'forestgreen', 'brown', 'darkviolet']

if __name__ == "__main__":
    main()
# %%
