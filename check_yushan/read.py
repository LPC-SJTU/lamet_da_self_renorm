# %%
import numpy as np
import gvar as gv
import matplotlib.pyplot as plt


x_ls = np.arange(-2-0.01, 3.02, 0.01)

# read data
mom_n1_lic_da = gv.load('lc_mom6_ls')
mom_n2_lic_da = gv.load('lc_mom8_ls')
mom_n3_lic_da = gv.load('lc_mom10_ls')

### example ###
### var = [1,2,3] ### 随便一个变量
### gv.dump(var, './filename') ### 存到 './filename'
### re_var = gv.load('./filename') ### 读回来


# sample average
mom_n1_lic_da = gv.dataset.avg_data(mom_n1_lic_da, bstrap=True)
mom_n2_lic_da = gv.dataset.avg_data(mom_n2_lic_da, bstrap=True)
mom_n3_lic_da = gv.dataset.avg_data(mom_n3_lic_da, bstrap=True)

fig = plt.figure()
ax = plt.axes()
    
ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n1_lic_da], [(val.mean - val.sdev) for val in mom_n1_lic_da], alpha=0.4, label=r'$a \to 0, Pz=1.29GeV$')
ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n2_lic_da], [(val.mean - val.sdev) for val in mom_n2_lic_da], alpha=0.4, label=r'$a \to 0, Pz=1.72GeV$')
ax.fill_between(x_ls, [(val.mean + val.sdev) for val in mom_n3_lic_da], [(val.mean - val.sdev) for val in mom_n3_lic_da], alpha=0.4, label=r'$a \to 0, Pz=2.15GeV$')

ax.fill_between(np.linspace(-0.5, 0.05, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.4)
ax.fill_between(np.linspace(0.95, 1.5, 500), np.ones(500)*-1, np.ones(500)*2, color='grey', alpha=0.4)
#ax.axvline(1, color='k', linestyle='--')
ax.axvline(0, color='k', linestyle='--')
ax.axvline(0.5, color='green', linestyle='--')
ax.axhline(0, color='k', linestyle='--')
ax.set_ylim([-0.19, 1.5])
ax.set_xlim([-0.5, 1.5])
ax.legend(loc='lower center')
ax.tick_params(direction='in')
plt.show()

# %%
