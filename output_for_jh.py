# %%
from head import *

meson = 'pion'
mom = 10

hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_re_ls') # shape = (N_conf, N_z)
hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a_hyb_im_ls')

a_re_avg = gv.dataset.avg_data(hyb_re_ls, bstrap=True)
a_im_avg = gv.dataset.avg_data(hyb_im_ls, bstrap=True)


pz = int(mom) * 0.215
lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm


fig = plt.figure(figsize=fig_size)
ax = plt.axes(plt_axes)

ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_re_avg], [(val.mean - val.sdev) for val in a_re_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

ax.axhline(0, color='k', linestyle='--')
ax.legend(loc='upper right')
ax.set_xlabel(lambda_label, **fs_p)
ax.set_ylim([-0.5, 1.25])
ax.tick_params(direction='in', **ls_p)
ax.set_title(hyb_re_label, **fs_p)
plt.show()




fig = plt.figure(figsize=fig_size)
ax = plt.axes(plt_axes)

ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_im_avg], [(val.mean - val.sdev) for val in a_im_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

ax.axhline(0, color='k', linestyle='--')
ax.legend(loc='upper right')
ax.set_xlabel(lambda_label, **fs_p)
ax.set_ylim([-1.1, 0.5])
ax.tick_params(direction='in', **ls_p)
ax.set_title(hyb_im_label, **fs_p)
plt.show()


# %%
f = open('./continuum_limit.txt', 'w')
line = []
line.append('z')
line.append('\t')
line.append('lambda')
line.append('\t')
line.append('Re')
line.append('\t')
line.append('Re_err')
line.append('\t')
line.append('Im')
line.append('\t')
line.append('Im_err')
line.append('\n')
f.writelines(line)

for idz in range(len(z_ls_da)):
    line = []
    line.append(str(z_ls_da[idz]))
    line.append('\t')
    line.append(str(lam_ls[idz]))
    line.append('\t')
    line.append(str(a_re_avg[idz].mean))
    line.append('\t')
    line.append(str(a_re_avg[idz].sdev))
    line.append('\t')
    line.append(str(a_im_avg[idz].mean))
    line.append('\t')
    line.append(str(a_im_avg[idz].sdev))
    line.append('\n')
    f.writelines(line)
f.close


