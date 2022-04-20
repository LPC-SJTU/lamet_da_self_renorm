# %%
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt

################
### for pion ###
################
x_ls = gv.load('p_no_sys_err_x')
y_no_sys = gv.load('p_no_sys_err_y')
yerr_no_sys = gv.load('p_no_sys_err_yerr') # statistic error 

print(np.shape(x_ls))

##################################
## sys error of large mom limit ##
##################################
mom_y = gv.load('p_mom_y')
mom_sys_ls = []
for idx in range(len(x_ls)):
    mom_sys = abs(y_no_sys[idx] - mom_y[idx])
    mom_sys_ls.append(mom_sys)

print(np.shape(mom_y))

################################
## sys error of extrapolation ##
################################
ext_y = gv.load('p_dif_ext_y')

ext_sys_ls = []
for idx in range(len(x_ls)):
    ext_sys = abs(y_no_sys[idx] - ext_y[idx])
    ext_sys_ls.append(ext_sys)

print(np.shape(ext_y))

##################################
## sys error of continuum limit ##
##################################
con_y = gv.load('p_a06_y')

con_sys_ls = []
for idx in range(len(x_ls)):
    con_sys = abs(y_no_sys[idx] - con_y[idx])
    con_sys_ls.append(con_sys)

print(np.shape(con_y))

#####################
## sys error of mu ##
#####################
mu_y = gv.load('p_mu_y')

mu_sys_ls = []
for idx in range(len(x_ls)):
    mu_sys = abs(y_no_sys[idx] - mu_y[idx])
    mu_sys_ls.append(mu_sys)

print(np.shape(mu_y))

############################
## add sys error together ##
############################
y1 = np.array([( y_no_sys[id] + np.sqrt(yerr_no_sys[id]**2 + mom_sys_ls[id]**2 + ext_sys_ls[id]**2 + con_sys_ls[id]**2 + mu_sys_ls[id]**2) ) for id in range(len(x_ls)) ]) 

y2 = np.array([( y_no_sys[id] - np.sqrt(yerr_no_sys[id]**2 + mom_sys_ls[id]**2 + ext_sys_ls[id]**2 + con_sys_ls[id]**2 + mu_sys_ls[id]**2) ) for id in range(len(x_ls)) ])

x_ls = np.hstack((x_ls, np.array([1]))) 
y1 = np.hstack((y1, np.array([0])))
y2 = np.hstack((y2, np.array([0])))

plt.fill_between(x_ls, y1, y2, color='red', alpha=0.5)
plt.show()


# %%
################
### for kaon ###
################
x_ls = gv.load('k_no_sys_err_x')
y_no_sys = gv.load('k_no_sys_err_y')
yerr_no_sys = gv.load('k_no_sys_err_yerr') # statistic error 

print(np.shape(x_ls))

##################################
## sys error of large mom limit ##
##################################
mom_y = gv.load('k_mom_y')
mom_sys_ls = []
for idx in range(len(x_ls)):
    mom_sys = abs(y_no_sys[idx] - mom_y[idx])
    mom_sys_ls.append(mom_sys)

print(np.shape(mom_y))

################################
## sys error of extrapolation ##
################################
ext_y = gv.load('k_dif_ext_y')

ext_sys_ls = []
for idx in range(len(x_ls)):
    ext_sys = abs(y_no_sys[idx] - ext_y[idx])
    ext_sys_ls.append(ext_sys)

print(np.shape(ext_y))

##################################
## sys error of continuum limit ##
##################################
con_y = gv.load('k_a06_y')

con_sys_ls = []
for idx in range(len(x_ls)):
    con_sys = abs(y_no_sys[idx] - con_y[idx])
    con_sys_ls.append(con_sys)

print(np.shape(con_y))

#####################
## sys error of mu ##
#####################
mu_y = gv.load('k_mu_y')

mu_sys_ls = []
for idx in range(len(x_ls)):
    mu_sys = abs(y_no_sys[idx] - mu_y[idx])
    mu_sys_ls.append(mu_sys)

print(np.shape(mu_y))

############################
## add sys error together ##
############################
y1 = np.array([( y_no_sys[id] + np.sqrt(yerr_no_sys[id]**2 + mom_sys_ls[id]**2 + ext_sys_ls[id]**2 + con_sys_ls[id]**2 + mu_sys_ls[id]**2) ) for id in range(len(x_ls)) ]) 

y2 = np.array([( y_no_sys[id] - np.sqrt(yerr_no_sys[id]**2 + mom_sys_ls[id]**2 + ext_sys_ls[id]**2 + con_sys_ls[id]**2 + mu_sys_ls[id]**2) ) for id in range(len(x_ls)) ])

x_ls = np.hstack((x_ls, np.array([1]))) 
y1 = np.hstack((y1, np.array([0])))
y2 = np.hstack((y2, np.array([0])))

plt.fill_between(x_ls, y1, y2, color='red', alpha=0.5)
plt.show()


# %%
