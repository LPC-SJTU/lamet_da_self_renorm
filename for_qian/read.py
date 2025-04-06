# %%
import gvar as gv
import numpy as np

#! coordinate space quasi, after renormalization, before continuum limit, all lattice spacing are interpolated to the same lambda list
#* momentum = 6, 8, 10
#* lattice spacing = 0.06, 0.09, 0.12
#* N bootstrap samples = 600
#* lambda_ls: the lambda list, shape = (23)
#* pion_mom_re: the real part of the quasi, shape = ( 3 lattice spacing, N of bs samples = 600, len(lambda_ls) = 23 )
#* pion_mom_im: the imaginary part of the quasi, shape = ( 3 lattice spacing, N of bs samples = 600, len(lambda_ls) = 23 )

lambda_ls = gv.load("dump/pion_lam_ls")
pion_mom6_re = gv.load("dump/pion_mom6_re_lam_ls")
pion_mom6_im = gv.load("dump/pion_mom6_im_lam_ls")
print(np.shape(lambda_ls))
print(np.shape(pion_mom6_re))

import matplotlib.pyplot as plt
plt.plot( lambda_ls, np.mean(pion_mom6_re[0], axis=0), label='a=0.06' )
plt.plot( lambda_ls, np.mean(pion_mom6_re[1], axis=0), label='a=0.09' )
plt.plot( lambda_ls, np.mean(pion_mom6_re[2], axis=0), label='a=0.12' )
plt.title('Real part of the quasi at mom=6')
plt.xlabel('lambda')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
plt.show()


plt.plot( lambda_ls, np.mean(pion_mom6_im[0], axis=0), label='a=0.06' )
plt.plot( lambda_ls, np.mean(pion_mom6_im[1], axis=0), label='a=0.09' )
plt.plot( lambda_ls, np.mean(pion_mom6_im[2], axis=0), label='a=0.12' )
plt.title('Imag part of the quasi at mom=6')
plt.xlabel('lambda')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
plt.show()


# %%
#! coordinate space quasi, after continuum limit, rotate and back, waiting for extrapolation
#* momentum = 6, 8, 10
#* lambda_ls: the lambda list, shape = (23)
#* pion_mom6_quasi_re: the real part of the quasi, shape = ( N of bs samples = 600, len(lambda_ls) = 23 )
#* pion_mom6_quasi_im: the imaginary part of the quasi, shape = ( N of bs samples = 600, len(lambda_ls) = 23 )

lambda_ls = gv.load("dump/pion_lam_ls")
pion_mom6_quasi_re = gv.load("dump/pion_mom6_quasi_re_ls")
pion_mom6_quasi_im = gv.load("dump/pion_mom6_quasi_im_ls")
print(np.shape(lambda_ls))
print(np.shape(pion_mom6_quasi_re))

import matplotlib.pyplot as plt
plt.plot( lambda_ls, np.mean(pion_mom6_quasi_re, axis=0), label='mom=6' )
plt.title('Real part of the quasi after continuum limit with mom=6')
plt.xlabel('lambda')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
plt.show()

plt.plot( lambda_ls, np.mean(pion_mom6_quasi_im, axis=0), label='mom=6' )
plt.title('Imag part of the quasi after continuum limit with mom=6')
plt.xlabel('lambda')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
plt.show()



# %%
#! coordinate space quasi, after continuum limit, after extrapolation, waiting for FT
#* momentum = 6, 8, 10
#* lambda_ls_ex: the lambda list after extrapolation, shape = (2001,)
#* pion_mom6_quasi_ext_ls: the real part of the quasi after extrapolation, shape = ( N of bs samples = 600, len(lambda_ls_ex) = 2001 )

pion_mom6_lam_ls_ex = gv.load("dump/pion_mom6_lam_ls_ex")
pion_mom8_lam_ls_ex = gv.load("dump/pion_mom8_lam_ls_ex")
pion_mom10_lam_ls_ex = gv.load("dump/pion_mom10_lam_ls_ex")

pion_mom6_quasi_ext_ls = gv.load("dump/pion_mom6_quasi_ext_ls")
pion_mom8_quasi_ext_ls = gv.load("dump/pion_mom8_quasi_ext_ls")
pion_mom10_quasi_ext_ls = gv.load("dump/pion_mom10_quasi_ext_ls")

import matplotlib.pyplot as plt
plt.plot( pion_mom6_lam_ls_ex, np.mean(np.real(pion_mom6_quasi_ext_ls), axis=0), label='mom=6' )
plt.title('real quasi in the coor space after extrapolation with mom=6')
plt.xlabel('lambda')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
# plt.xlim(-20, 20)
plt.show()

plt.plot( pion_mom6_lam_ls_ex, np.mean(np.imag(pion_mom6_quasi_ext_ls), axis=0), label='mom=6' )
plt.title('imag quasi in the coor space after extrapolation with mom=6')
plt.xlabel('lambda')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
# plt.xlim(-20, 20)
plt.show()


# %%
#! momentum space quasi, after continuum limit, after extrapolation and FT, waiting for inverse matching
#* momentum = 6, 8, 10
#* x_ls: the x list, shape = (503)
#* pion_mom6_mom_ls: the quasi in the momentum space, shape = ( N of bs samples = 600, len(x_ls) = 503 )

x_ls = np.arange(-2-0.01, 3.02, 0.01)

pion_mom6_mom_ls = gv.load("dump/pion_mom6_quasi_mom_ls")
print(np.shape(x_ls))
print(np.shape(pion_mom6_mom_ls))

import matplotlib.pyplot as plt
# plt.plot( x_ls, np.mean(pion_mom6_mom_ls, axis=0), label='mom=6' )
plt.errorbar( x_ls, gv.mean(pion_mom6_mom_ls), yerr=gv.sdev(pion_mom6_mom_ls), label='mom=6' )
plt.title('quasi in the momentum space with mom=6')
plt.xlabel('x')
plt.legend()
plt.tick_params(direction="in", top="on", right="on")
plt.grid(linestyle=":")
plt.show()



# %%
