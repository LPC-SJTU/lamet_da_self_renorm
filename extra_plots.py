# %%
from head import *
from plot import *

def hyb_vs_RIMOM_plot():
    meson = 'pion'
    mom = 10

    rmom_a06 = []

    with open("RIMOM.txt", "r") as f:
        for line in f.readlines():
            line = line.strip()  #去掉列表中每一个元素的换行符
            line = line.split()
            if len(line) != 0 and line[0] == '0.057':
                rmom_a06.append(float(line[1]))

    rmom_z = []
    for i in range(len(rmom_a06)):
        rmom_z.append(round(0.06 * i, 2))

    rmom_z = np.array(rmom_z)

    da_a06_conf_ls = gv.load(meson+'/mom='+str(mom)+'/da_a06_conf_ls')

    def renorm(da_conf_ls, key): 
        da_re_ls = []
        da_im_ls = []
        for n_conf in range(len(da_conf_ls)):
            temp_re_ls = []
            temp_im_ls = []
            for idx in range(len(da_conf_ls[0])):
                temp_val = (da_conf_ls[n_conf][idx] / rmom_a06[idx+1])
                temp_re_ls.append(temp_val.real)
                temp_im_ls.append(temp_val.imag)

            da_re_ls.append(temp_re_ls)
            da_im_ls.append(temp_im_ls)

        return da_re_ls, da_im_ls

    lambda_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * 10 * gev_fm ) / gev_fm

    a06_hyb_re_ls, a06_hyb_im_ls = renorm(da_a06_conf_ls, key='a=0.0574')
    rmom_a06_lambda = turn_into_lambda(a06_hyb_re_ls, z_ls_extend, lambda_ls, 0.0574, 96, mom) 
    rmom_a06_lambda_im = turn_into_lambda(a06_hyb_im_ls, z_ls_extend, lambda_ls, 0.0574, 96, mom) 

    rmom_a06_re_ro, rmom_a06_im_ro = rotate(rmom_a06_lambda, rmom_a06_lambda_im, lambda_ls, back=False)

    rmom_a06_lambda_avg = gv.dataset.avg_data(rmom_a06_re_ro, bstrap=True)


    ### 

    a06_hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/a06_hyb_re_ls')
    a06_hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/a06_hyb_im_ls')

    pz = int(mom) * 0.215
    lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm

    a06_re_ro, a06_im_ro = rotate(a06_hyb_re_ls, a06_hyb_im_ls, lam_ls, back=False)

    a06_ro_avg = gv.dataset.avg_data(a06_re_ro, bstrap=True)


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(lam_ls, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='Hyrbid on self_R', fmt='o', **errorb)
    ax.errorbar(lam_ls, [val.mean for val in rmom_a06_lambda_avg], [val.sdev for val in rmom_a06_lambda_avg], color='red', label='Hybrid on RI/MOM', fmt='D', **errorb)
    ax.axhline(0, color='k', linestyle='--')
    ax.legend(loc='upper right')
    ax.set_xlabel(lambda_label, **fs_p)
    ax.set_ylim([-0.5, 1.25])
    ax.tick_params(direction='in', **ls_p)
    ax.set_title(hyb_ro_re_label, **fs_p)
    plt.savefig(meson+'/paper/renorm_comparison_mom_10_a06.pdf', transparent=True)
    plt.show()


# %%
from pdf_self_renorm import pdf_zR
zR_dic, m_pdf_dic = pdf_zR()

meson = 'pion'
mom = 8

da_a06_conf_ls = gv.load(meson+'/mom='+str(mom)+'/da_a06_conf_ls')
da_a09_conf_ls = gv.load(meson+'/mom='+str(mom)+'/da_a09_conf_ls')
da_a12_conf_ls = gv.load(meson+'/mom='+str(mom)+'/da_a12_conf_ls')

# z_ls_extend

def renorm(da_conf_ls, key): 
    da_re_ls = []
    da_im_ls = []
    for n_conf in range(len(da_conf_ls)):
        temp_re_ls = []
        temp_im_ls = []
        for idx in range(len(da_conf_ls[0])):
            temp_val = (da_conf_ls[n_conf][idx] / zR_dic[key][idx].mean) 
            temp_re_ls.append(temp_val.real)
            temp_im_ls.append(temp_val.imag)

        da_re_ls.append(temp_re_ls)
        da_im_ls.append(temp_im_ls)

    return da_re_ls, da_im_ls

a06_hyb_re_ls, a06_hyb_im_ls = renorm(da_a06_conf_ls, key='a=0.0574')
a09_hyb_re_ls, a09_hyb_im_ls = renorm(da_a09_conf_ls, key='a=0.0882')
a12_hyb_re_ls, a12_hyb_im_ls = renorm(da_a12_conf_ls, key='a=0.1213')

def lam(z_ls, a, L, mom):
    pz = 2*np.pi / (a * L) * mom * gev_fm # pz=2pi / (0.09*64) * mom=8 * 0.197
    return z_ls * pz / gev_fm

a06_re_ro, a06_im_ro = rotate(a06_hyb_re_ls, a06_hyb_im_ls, lam(z_ls_extend, 0.0574, 96, mom), back=False)
a09_re_ro, a09_im_ro = rotate(a09_hyb_re_ls, a09_hyb_im_ls, lam(z_ls_extend, 0.0882, 64, mom), back=False)
a12_re_ro, a12_im_ro = rotate(a12_hyb_re_ls, a12_hyb_im_ls, lam(z_ls_extend, 0.1213, 48, mom), back=False)

a06_ro_avg = gv.dataset.avg_data(a06_re_ro, bstrap=True)
a09_ro_avg = gv.dataset.avg_data(a09_re_ro, bstrap=True)
a12_ro_avg = gv.dataset.avg_data(a12_re_ro, bstrap=True)


fig = plt.figure(figsize=fig_size)
ax = plt.axes(plt_axes)
ax.errorbar(z_ls_extend, [val.mean for val in a06_ro_avg], [val.sdev for val in a06_ro_avg], color=color_list[0], label='a:0.06fm', fmt='o', **errorb)
ax.errorbar(z_ls_extend, [val.mean for val in a09_ro_avg], [val.sdev for val in a09_ro_avg], color=color_list[1], label='a:0.09fm', fmt='D', **errorb)
ax.errorbar(z_ls_extend, [val.mean for val in a12_ro_avg], [val.sdev for val in a12_ro_avg], color=color_list[2], label='a:0.12fm', fmt='s', **errorb)
ax.plot(z_ls_extend, [ZMS_da(val) for val in z_ls_extend], 'r-', label='ZMS-bar')

ax.axhline(0, color='k', linestyle='--')
ax.legend(loc='upper right')
ax.set_xlabel(lambda_label, **fs_p)
ax.set_ylim([-0.5, 2])
ax.tick_params(direction='in', **ls_p)
ax.set_title(hyb_ro_re_label, **fs_p)
plt.savefig(meson+'/paper/compare_with_ZMS-bar_mom_8.pdf', transparent=True)
plt.show()

# %%
