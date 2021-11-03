# %%
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt

x_ls = gv.load('x_ls') # 横坐标 x
mom_10_lcda = gv.load('mom_10_lcda') # mom=10, 用于添加系统误差，gvar 变量
lmom_lcda = gv.load('large_mom_lcda') # 大动量极限结果，gvar 变量

delta_ls = []
for idx in range(len(x_ls)):
        delta = abs(lmom_lcda[idx].mean - mom_10_lcda[idx].mean) # system error
        delta_ls.append(delta)

y1 = np.array([(val.mean + val.sdev) for val in lmom_lcda]) + np.array(delta_ls)
y2 = np.array([(val.mean - val.sdev) for val in lmom_lcda]) - np.array(delta_ls)

plt.figure()
plt.fill_between(x_ls, y1, y2, color='orange', alpha=0.5)

# %%
x_ls = gv.load('x_ls_fit')
mom_10_lcda = gv.load('mom_10_lcda_fit')
lmom_lcda = gv.load('large_mom_lcda_fit')

delta_ls = []
for idx in range(len(x_ls)):
        delta = abs(lmom_lcda[idx].mean - mom_10_lcda[idx].mean) # system error
        delta_ls.append(delta)

y1 = np.array([(val.mean + val.sdev) for val in lmom_lcda]) + np.array(delta_ls)
y2 = np.array([(val.mean - val.sdev) for val in lmom_lcda]) - np.array(delta_ls)

plt.figure()
plt.fill_between(x_ls, y1, y2, color='orange', alpha=0.5)
# %%
