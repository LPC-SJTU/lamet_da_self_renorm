import multiprocessing
import timeit
import numpy as np
import gvar as gv
from tqdm import tqdm
from scipy import integrate
from scipy import interpolate

mu = 2
lms = 0.24451721864451428  #Lambda_MS
cf=4/3
nf=3
b0=11-2/3*nf
alphas=2 * np.pi/( b0 * np.log(mu/lms) )

mom_to_pz = 0.215 # mom=8 corresponding to pz=1.72, pz=2pi / (0.09*64) * mom=8 * 0.197
alphas_cf_div_2pi = alphas * cf / (2 * np.pi)


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
        p = multiprocessing.Pool(2)#4 multiprocess to integral
        lc_ls = p.map(self.integral, n_conf)
        p.close()
        p.join()

        # lc_ls = []
        # for n_conf in tqdm(np.arange(len(self.quasi_re_ls))):
        #     lc_ls.append(self.integral(n_conf))

        end = timeit.default_timer()
        print('multi processing time:', str(end-start), 's')

        lc_re_ls = np.zeros_like(lc_ls)
        lc_im_ls = np.zeros_like(lc_ls)
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

if __name__ == '__main__':
    lam_ls = gv.load('lam_ls')
    quasi_re_ls = gv.load('quasi_ro_re_ls')
    quasi_im_ls = gv.load('quasi_ro_im_ls')

    inv_matching_coor = INV_MATCHING_COOR(2.15, lam_ls, quasi_re_ls, quasi_im_ls)
    lc_re_ls, lc_im_ls = inv_matching_coor.main()

