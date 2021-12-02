

# %%
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
fig_size_lc = (fig_width * 0.8, fig_width * 0.8)
plt_axes = [0.1,0.12,0.85,0.8]
fs_p = {"fontsize": 13} # font size of text, label, ticks
ls_p = {"labelsize": 13}

errorb = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle

# labels
lambda_label = r"$\lambda = z P_z$"
x_label = r"$x$"
z_label = r"$z(fm)$"
t_label = r"$t$"
meff_local_label = r"$\ln(C(z=0, t) / C(z=0, t+1))$"
meff_non_local_label = r"$\ln(C(z, t) / C(z, t+1))$"
ratio_label = r"$C(z, t) / C(z=0, t)$"
hyb_ro_re_label = r'$Re[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
hyb_ro_im_label = r'$Im[e^{\frac{i z P_z}{2}} H_{\pi}(z)]$'
hyb_re_label = r'$Re[H_{\pi}(z)]$'
hyb_im_label = r'$Im[H_{\pi}(z)]$'

fig = plt.figure(figsize=fig_size)
ax = plt.axes(plt_axes)
ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in quasi_re_avg], [(val.mean - val.sdev) for val in quasi_re_avg], color='blue', alpha=0.7, label='quasi')
ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in lc_re_avg], [(val.mean - val.sdev) for val in lc_re_avg], color='red', alpha=0.7, label='lc')
ax.axhline(0, color='k', linestyle='--')
ax.legend(loc='upper right')
ax.set_xlabel(lambda_label, **fs_p)
ax.set_ylim([-0.5, 1.25])
ax.tick_params(direction='in', **ls_p)
plt.show()