# %%
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt
from head import *


x_ls = gv.load('x_ls_fit')
mom_10_lcda = gv.load('mom_10_lcda_fit')
lmom_lcda = gv.load('large_mom_lcda_fit')
print(np.array(lmom_lcda).shape)

delta_ls = []
for idx in range(len(x_ls)):
        delta = abs(lmom_lcda[idx].mean - mom_10_lcda[idx].mean) # system error
        delta_ls.append(delta)

y1 = np.array([(val.mean + val.sdev) for val in lmom_lcda]) + np.array(delta_ls)
y2 = np.array([(val.mean - val.sdev) for val in lmom_lcda]) - np.array(delta_ls)
ymean = np.array([val.mean for val in lmom_lcda])


x_ls_2 = gv.load('x_ls') # 横坐标 x
mom_10_lcda_2 = gv.load('mom_10_lcda') # mom=10, 用于添加系统误差，gvar 变量
lmom_lcda_2 = gv.load('large_mom_lcda') # 大动量极限结果，gvar 变量

delta_ls_2 = []
for idx in range(len(x_ls_2)):
        delta = abs(lmom_lcda_2[idx].mean - mom_10_lcda_2[idx].mean) # system error
        delta_ls_2.append(delta)

y1_2 = np.array([(val.mean + val.sdev) for val in lmom_lcda_2]) + np.array(delta_ls_2)
y2_2 = np.array([(val.mean - val.sdev) for val in lmom_lcda_2]) - np.array(delta_ls_2)
ymean_2 = np.array([val.mean for val in lmom_lcda_2])

###range -0.2~0.05
s1_1 = np.zeros(26)
s1_1[0:20] = np.abs(y1_2[180:200])
s1_1[20:26] = y1[0:6] + np.abs(ymean_2[200:206] - ymean[0:6]) #error_max + error_fit
s1_2 = np.zeros(26)
for i in range(0,20):
        if y2_2[180+i] > 0:     s1_2[i] = 0#if error_min >0, cover to 0
        else:   s1_2[i] = y2_2[180+i]#if error_min, take error_min
s1_2[20:26] = y2[0:6] - np.abs(ymean_2[200:206] - ymean[0:6]) #error_max - error_fit

###range 0.95~1.2
s2_1 = np.zeros(26)
s2_1[0:6] = y1[95:101] + np.abs(ymean_2[295:301] - ymean[95:101]) 
s2_1[6:26] = np.abs(y1_2[301:321])
s2_2 = np.zeros(26)
s2_2[0:6] = y2[95:101] - np.abs(ymean_2[295:301] - ymean[95:101]) 
for i in range(0,20):
        if y2_2[301+i] > 0:     s2_2[6+i] = 0
        else:   s2_2[6+i] = y2_2[301+i]


x1 = np.arange(-0.2, 0.06, 0.01)
print(x1.shape)
x2 = np.arange(0.95, 1.21, 0.01)

meson = 'pion'


#################
fig = plt.figure(figsize=fig_size_lc_err)
ax = plt.axes(plt_axes)

a1 = gv.gvar(-0.06, 0.03) # sum rule
a2 = gv.gvar(0.25, 0.15) # same for pion and kaon
a4 = gv.gvar(-0.015, 0.025)


ax.fill_between(x_ls, [sum_rule(meson, x, a1, a2, a4).mean + sum_rule(meson, x, a1, a2, a4).sdev for x in x_ls], [sum_rule(meson, x, a1, a2, a4).mean - sum_rule(meson, x, a1, a2, a4).sdev for x in x_ls], color=color_list[1], label='Sum rule', alpha=0.4)

ax.plot(x_ls, [6*x*(1-x) for x in x_ls], color='red', linestyle='dashdot', label='Asymptotic') # only plot between 0 and 1

ax.fill_between(x1, s1_1, s1_2, color='darkseagreen', alpha=0.6, label='Systematic error') # sys err
ax.fill_between(x2, s2_1, s2_2, color='darkseagreen', alpha=0.6) # sys err
#ax.fill_between(x_ls, y1, y2, color='orange', alpha=0.5)


ax.fill_between(x_ls, y1, y2, color=color_list[0], alpha=0.5)
ax.plot(x_ls, (y1+y2)/2, color=color_list[0], label='This work', linewidth=2, linestyle='dotted')

ax.fill_between(x1, np.ones(len(x1))*-1, np.ones(len(x1))*2, color='grey', alpha=0.4)
ax.fill_between(x2, np.ones(len(x2))*-1, np.ones(len(x2))*2, color='grey', alpha=0.4)

## grey v band to cover fit region

ax.axvline(0.5, color='green', linestyle='--')
ax.axvline(0, color='k', linestyle='--')
#ax.axvline(1, color='k', linestyle='--')
ax.axhline(0, color='k', linestyle='--')
#ax.set_title('DA light-cone Pz to infty', **fs_p)
ax.set_xlabel(x_label, **fs_p)
ax.set_ylim([-0.19, 1.7])
ax.set_xlim([-0.2, 1.2])
ax.legend(loc='lower center')
ax.tick_params(direction='in', **ls_p)
plt.savefig('lcda_Pz_to_infty_err.pdf', transparent=True)
plt.show()
########################




# %%
