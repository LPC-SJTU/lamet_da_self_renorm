#coding=utf-8

import numpy as np
from numpy.core.defchararray import count
import math as m
import gvar as gv
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
from scipy import integrate
from tqdm import tqdm
import multiprocessing
import timeit


def main():
    MT = MATCHING_PROCESS()
    quasi_re_ls, quasi_im_ls, lambda_ls = MT.read(plot=True)
    MT.matching(quasi_re_ls, quasi_im_ls, lambda_ls, plot=True)

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
        quasi_re = np.array(gv.load('quasi_re'))
        quasi_im = np.array(gv.load('quasi_im'))
        lam_ls = np.array(gv.load('lambda_ls'))

        q_re_avg = gv.dataset.avg_data(quasi_re, bstrap=True)
        q_im_avg = gv.dataset.avg_data(quasi_im, bstrap=True)
        if plot==True:  plt.fill_between(lam_ls, [(val.mean + val.sdev) for val in q_re_avg], [(val.mean - val.sdev) for val in q_re_avg], color='pink', alpha=0.8)
        #plt.show()
        return quasi_re, quasi_im, lam_ls

    def matching(self, quasi_re_ls, quasi_im_ls, lambda_ls, plot=False):
        self.quasi_re_ls = quasi_re_ls
        self.quasi_im_ls = quasi_im_ls
        self.lambda_ls = lambda_ls

        ### multi process ###
        start = timeit.default_timer()
        #items_cfg = [x for x in range(len(quasi_re_ls))]
        items_cfg = tqdm(np.arange(0,60))
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
        lc_re_avg = gv.dataset.avg_data(lc_re, bstrap=True)
        lc_im_avg = gv.dataset.avg_data(lc_im, bstrap=True)
        if plot==True:  plt.fill_between(self.lambda_ls, [(val.mean + val.sdev) for val in lc_re_avg], [(val.mean - val.sdev) for val in lc_re_avg], color='green', alpha=0.5)
        print(lc_re_avg)
        plt.show()
        #l_re_avg = gv.dataset.avg_data(quasi_re, bstrap=True)
        #l_im_avg = gv.dataset.avg_data(quasi_im, bstrap=True)

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



if __name__ == "__main__":
    main()