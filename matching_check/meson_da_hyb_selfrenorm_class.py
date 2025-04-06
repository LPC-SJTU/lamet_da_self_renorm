# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import timeit
from tqdm import tqdm
import gvar as gv


mom_to_pz = 0.215 # mom=8 corresponding to pz=1.72, pz=2pi / (0.09*64) * mom=8 * 0.197

lqcd = 0.1 #Lambda_QCD
lms = 0.24451721864451428  #Lambda_MS
k = 3.320
d_pdf = -0.08183 #-0.1252 # for pdf
d_da = 0.19 # 0.1 for a^2 order, 0.19 for a order, from Yushan
m0_da = gv.gvar(-0.094, 0.024) # from Yushan
mu = 2 # GeV, for renormalization #!# for sys err
mu_f = 2 # GeV, for factorization
cf=4/3
nf=3
b0=11-2/3*nf
alphas=2 * np.pi/( b0 * np.log(mu_f/lms) )

mom_to_pz = 0.215 # mom=8 corresponding to pz=1.72, pz=2pi / (0.09*64) * mom=8 * 0.197
alphas_cf_div_2pi = alphas * cf / (2 * np.pi)


def interp_1d(x_in, y_in, x_out, method="linear"): # interpolation
    f=interpolate.interp1d(x_in, y_in, kind=method)
    y_out = f(x_out)

    return y_out

# %%  
class INV_MATCHING_COOR():
    def __init__(self, pz, lam_ls, quasi_re_ls, quasi_im_ls):
        self.pz = pz
        self.lambda_ls = lam_ls
        self.quasi_re_ls = quasi_re_ls
        self.quasi_im_ls = quasi_im_ls

    def main(self):
        start = timeit.default_timer()
        print('>>> integrating for matching of mom '+str(self.pz / mom_to_pz))

        from joblib import Parallel, delayed
        lc_ls = Parallel(n_jobs=8)(delayed(self.integral)(n_conf) for n_conf in tqdm(range(len(self.quasi_re_ls)), desc='Processing configurations'))

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
            part2 = integrate.quad(fp_2_re, 0, lam-delta_l, args=lam, full_output=1)[0] + 1j * integrate.quad(fp_2_im, 0, lam-delta_l, args=lam, full_output=1)[0]
            part3 = integrate.quad(fp_3_re, 0, lam-delta_l, args=lam, full_output=1)[0] + 1j * integrate.quad(fp_3_im, 0, lam-delta_l, args=lam, full_output=1)[0]
            part4 = integrate.quad(fp_4_re, 0, lam, args=lam, full_output=1)[0] + 1j * integrate.quad(fp_4_im, 0, lam, args=lam, full_output=1)[0]

            lc_re = lc_re - (part1+part2+part3+part4).real
            lc_im = lc_im - (part1+part2+part3+part4).imag

            lc_ls.append(lc_re + 1j * lc_im)

        return lc_ls

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
        
        #!: add two-loop correction
        Z2m = np.loadtxt("./C2loop_crrctn.csv", delimiter=",")
        C_matrix = C_matrix + Z2m

        C_matrix_inverse = np.linalg.inv(C_matrix)
        
        print("shape of C_matrix_inverse: ", C_matrix_inverse.shape)

        return C_matrix_inverse



# %%
