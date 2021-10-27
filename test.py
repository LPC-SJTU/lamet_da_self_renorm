# %%
import gvar as gv

lam_ls_ex = gv.load('pion'+'/mom='+str(6)+'/lam_ls_ex') 
hyb_complex = gv.load('pion'+'/mom='+str(6)+'/hyb_complex') # before ft
hyb_re_ft = gv.load('pion'+'/mom='+str(6)+'/hyb_re_ft') # shape = (N_re, len(x_ls))
hyb_im_ft = gv.load('pion'+'/mom='+str(6)+'/hyb_im_ft')

print(lam_ls_ex)
# %%
