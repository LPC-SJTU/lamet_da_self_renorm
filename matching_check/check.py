# %%
from meson_da_hyb_selfrenorm_class import *

meson = 'pion'
mom = 8
pz = mom * mom_to_pz


#! Matching in the momentum space
x_ls = np.arange(-2-0.01, 3.02, 0.01) # x after ft, for quasi before matching
print(len(x_ls))
y_ls = x_ls

# delta = 0.00001
# x_ls_mat = np.linspace(x_ls[0] + delta, x_ls[-1] - delta, 503)
# y_ls_mat = np.linspace(x_ls[0], x_ls[-1], 503)

x_ls_mat = x_ls
y_ls_mat = x_ls

inv_matching_mom = INV_MATCHING_MOM(x_ls, x_ls_mat, y_ls_mat)

lam_ls_ex = gv.load('cache/lam_ls_ex')
quasi_ext_ls = gv.load('cache/quasi_ext_ls')
quasi_mom_ls = gv.load('cache/quasi_mom_ls')

quasi_mom_avg = gv.dataset.avg_data(quasi_mom_ls, bstrap=True)
lc_mom_ls = inv_matching_mom.main(pz, quasi_mom_ls)
lc_mom_avg = gv.dataset.avg_data(lc_mom_ls, bstrap=True)


fig, ax = plt.figure(), plt.axes()
ax.fill_between(x_ls, gv.mean(lc_mom_avg) + gv.sdev(lc_mom_avg), gv.mean(lc_mom_avg) - gv.sdev(lc_mom_avg), color='green', alpha=0.5)
ax.plot(x_ls, gv.mean(lc_mom_avg), color='blue', label='Mom matching')
plt.tight_layout()
plt.xlim(-0.5, 1.5)
plt.show()


# %%
#! Matching in the coordinate space

lam_ls = gv.load('cache/lam_ls')
quasi_re_ls = gv.load('cache/quasi_re_ls')
quasi_im_ls = gv.load('cache/quasi_im_ls')

inv_matching_coor = INV_MATCHING_COOR(pz, lam_ls, quasi_re_ls, quasi_im_ls)
lc_re_ls, lc_im_ls = inv_matching_coor.main()

lc_re_avg = gv.dataset.avg_data(lc_re_ls, bstrap=True)
lc_im_avg = gv.dataset.avg_data(lc_im_ls, bstrap=True)

fig, ax = plt.figure(), plt.axes()
ax.fill_between(lam_ls, gv.mean(lc_re_avg) + gv.sdev(lc_re_avg), gv.mean(lc_re_avg) - gv.sdev(lc_re_avg), color='green', alpha=0.5)
ax.plot(lam_ls, gv.mean(lc_re_avg), color='blue', label='Coor matching')
plt.tight_layout()
plt.show()


# %%
