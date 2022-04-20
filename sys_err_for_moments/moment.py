# %%
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt

from head import mellin_moment

################
### for pion ###
################
x_ls = np.array(gv.load('p_no_sys_err_x'))
y_no_sys = np.array(gv.load('p_no_sys_err_y'))
yerr_no_sys = np.array(gv.load('p_no_sys_err_yerr')) # statistic error 
# x_ls = np.append(x_ls, 1)
# y_no_sys = np.append(y_no_sys, 0)
# yerr_no_sys = np.append(yerr_no_sys, 0)
print(yerr_no_sys.shape)

##################################
## sys error of large mom limit ##
##################################
mom_y = np.array(gv.load('p_mom_y'))
mom_y = np.append(mom_y, 0)
mom_sys_ls = []
for idx in range(len(x_ls)):
    mom_sys = abs(y_no_sys[idx] - mom_y[idx])
    mom_sys_ls.append(mom_sys)
print(mom_y.shape)

################################
## sys error of extrapolation ##
################################
ext_y = np.array(gv.load('p_dif_ext_y'))
ext_sys_ls = []
for idx in range(len(x_ls)):
    ext_sys = abs(y_no_sys[idx] - ext_y[idx])
    ext_sys_ls.append(ext_sys)
print(ext_y.shape)

##################################
## sys error of continuum limit ##
##################################
con_y = np.array(gv.load('p_a06_y'))
con_sys_ls = []
for idx in range(len(x_ls)):
    con_sys = abs(y_no_sys[idx] - con_y[idx])
    con_sys_ls.append(con_sys)
print(con_y.shape)

#####################
## sys error of mu ##
#####################
mu_y = np.array(gv.load('p_mu_y'))
mu_sys_ls = []
for idx in range(len(x_ls)):
    mu_sys = abs(y_no_sys[idx] - mu_y[idx])
    mu_sys_ls.append(mu_sys)
print(mu_y.shape)

def gen_an(x, n):
    factor = 2*(2*n+3) / (3*(n+1)*(n+2))
    if n==1:
        return factor*3*(2*x-1)
    if n==2:
        return factor*3/2*(5*(2*x-1)**2-1)
    if n==3:
        return factor*5/2*(2*x-1)*(7*(2*x-1)**2-3)
    if n==4:
        return factor*15/8*(1 - 14*(2*x-1)**2 + 21*(2*x-1)**4)
    if n==5:
        return factor*21/8*(2*x-1)*(5 -30*(2*x-1)**2 + 33*(2*x-1)**4)
    if n==6:
        return factor*7/16*(-5 + 135*(2*x-1)**2 - 495*(2*x-1)**4 + 429*(2*x-1)**6)
    if n==8:
        return factor*45/128*(7 + 11*(2*x-1)**2*(-28 + 182*(2*x-1)**2 - 364*(2*x-1)**4 + 221*(2*x-1)**6))

##########################################
################# check ##################
print('Statistical Results')
sta_y = gv.gvar(y_no_sys, yerr_no_sys)
x = x_ls


print('>>> Jinchen:')
mellin_moment(x, sta_y, 2) #!#

print('>>> Jun:')
chi2 = (2*x-1)**2 * sta_y ###Xi2
print('Xi_2:' + str((np.sum(chi2)*0.01)))
an={}
for n in [1,2,3,4,6]:
    lim_all =  sta_y * gen_an(x, n )
    sum_lcda = np.sum(lim_all, axis=0)*0.01
    an['a'+str(n)+'_sys'] = sum_lcda
    print('a'+str(n)+':' + str(sum_lcda))


# %%
print('Systematic Error:')
sys_keys = ['mom_y', 'dif_ext_y', 'a06_y', 'mu_y']
for key in sys_keys:
    print('Kind:' + key)
    if key == 'mom_y':          wave = mom_y   
    elif key == 'dif_ext_y':    wave = ext_y
    elif key == 'a06_y':        wave = con_y
    elif key == 'mu_y':         wave = mu_y
    x = x_ls
    chi2 = (2*x-1)**2 * wave
    print('Xi_2:' + str('%.3f'%(np.sum(chi2)*0.01)))
    an={}
    for n in [1,2,3,4,6]:
        lim_all =  wave * gen_an(x, n )
        sum_lcda = np.sum(lim_all, axis=0)*0.01
        an['a'+str(n)+'_sys'] = sum_lcda
        print('a'+str(n)+':' + str('%.3f' %sum_lcda))



############################
## add sys error together ##
############################
#y1 = np.array([( y_no_sys[id] + np.sqrt(yerr_no_sys[id]**2 ) ) for id in range(len(x_ls)) ]) 
#y2 = np.array([( y_no_sys[id] - np.sqrt(yerr_no_sys[id]**2 ) ) for id in range(len(x_ls)) ])
y1 = np.array([( y_no_sys[id] + np.sqrt(yerr_no_sys[id]**2 + mom_sys_ls[id]**2 + ext_sys_ls[id]**2 + con_sys_ls[id]**2 + mu_sys_ls[id]**2) ) for id in range(len(x_ls)) ]) 
y2 = np.array([( y_no_sys[id] - np.sqrt(yerr_no_sys[id]**2 + mom_sys_ls[id]**2 + ext_sys_ls[id]**2 + con_sys_ls[id]**2 + mu_sys_ls[id]**2) ) for id in range(len(x_ls)) ])

#x_ls = np.hstack((x_ls, np.array([1]))) 
#y1 = np.hstack((y1, np.array([0])))
#y2 = np.hstack((y2, np.array([0])))

plt.fill_between(x_ls, y1, y2, color='red', alpha=0.5)
plt.show()


# %%
