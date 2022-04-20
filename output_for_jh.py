# %%
# from head import *

# meson = 'pion'
# mom = 10

# hyb_re_ls = gv.load(meson+'/mom='+str(mom)+'/lc_re_ls') # shape = (N_conf, N_z)
# hyb_im_ls = gv.load(meson+'/mom='+str(mom)+'/lc_im_ls')

# a_re_avg = gv.dataset.avg_data(hyb_re_ls, bstrap=True)
# a_im_avg = gv.dataset.avg_data(hyb_im_ls, bstrap=True)


# pz = int(mom) * 0.215
# lam_ls = z_ls_da * ( 2*np.pi / (0.0574*96) * mom * gev_fm ) / gev_fm


# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)

# ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_re_avg], [(val.mean - val.sdev) for val in a_re_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-0.5, 1.25])
# ax.tick_params(direction='in', **ls_p)
# ax.set_title(hyb_re_label, **fs_p)
# plt.show()




# fig = plt.figure(figsize=fig_size)
# ax = plt.axes(plt_axes)

# ax.fill_between(lam_ls, [(val.mean + val.sdev) for val in a_im_avg], [(val.mean - val.sdev) for val in a_im_avg], color='pink', alpha=0.8, label=r'$a\ \to \ 0$')

# ax.axhline(0, color='k', linestyle='--')
# ax.legend(loc='upper right')
# ax.set_xlabel(lambda_label, **fs_p)
# ax.set_ylim([-1.1, 0.5])
# ax.tick_params(direction='in', **ls_p)
# ax.set_title(hyb_im_label, **fs_p)
# plt.show()


# %%
# f = open('./lc.txt', 'w')
# line = []
# line.append('z')
# line.append('\t')
# line.append('lambda')
# line.append('\t')
# line.append('Re')
# line.append('\t')
# line.append('Re_err')
# line.append('\t')
# line.append('Im')
# line.append('\t')
# line.append('Im_err')
# line.append('\n')
# f.writelines(line)

# for idz in range(len(z_ls_da)):
#     line = []
#     line.append(str(z_ls_da[idz]))
#     line.append('\t')
#     line.append(str(lam_ls[idz]))
#     line.append('\t')
#     line.append(str(a_re_avg[idz].mean))
#     line.append('\t')
#     line.append(str(a_re_avg[idz].sdev))
#     line.append('\t')
#     line.append(str(a_im_avg[idz].mean))
#     line.append('\t')
#     line.append(str(a_im_avg[idz].sdev))
#     line.append('\n')
#     f.writelines(line)
# f.close


# %%
from head import *
from meson_da_hyb_selfrenorm_class import *
from pdf_self_renorm import pdf_zR

meson = 'kaon'

### bare matrix element ###
bare_dic = {}
for mom in [6, 8, 10]:
    temp = gv.load(meson+'/mom='+str(mom)+'/da_an_ls')
    bare_dic['mom'+str(mom)] = temp
    bare_dic['mom'+str(mom)+'re'] = [[], [], []]
    bare_dic['mom'+str(mom)+'im'] = [[], [], []]
    for ida in range(3):
        for n_conf in range(600):
            bare_dic['mom'+str(mom)+'re'][ida].append([])
            bare_dic['mom'+str(mom)+'im'][ida].append([])
            for idz in range(25):
                bare_dic['mom'+str(mom)+'re'][ida][n_conf].append(temp[ida][n_conf][idz].real)
                bare_dic['mom'+str(mom)+'im'][ida][n_conf].append(temp[ida][n_conf][idz].imag)

print(np.shape(bare_dic['mom6im']))

bare_dic_avg = {}
for mom in [6, 8, 10]:
    bare_dic_avg['mom'+str(mom)+'re'] = []
    bare_dic_avg['mom'+str(mom)+'im'] = []
    for ida in range(3):
        bare_dic_avg['mom'+str(mom)+'re'].append(gv.dataset.avg_data(bare_dic['mom'+str(mom)+'re'][ida], bstrap=True))
        bare_dic_avg['mom'+str(mom)+'im'].append(gv.dataset.avg_data(bare_dic['mom'+str(mom)+'im'][ida], bstrap=True))

print(np.shape(bare_dic_avg['mom6im'])) # [a06, a09, a12], z



### renormalized matrix element ###
zR_dic, m_pdf_dic = pdf_zR()
self_renormalization = RENORMALIZATION(zR_dic)

renorm_dic = {}
for mom in [6, 8, 10]:
    renorm_dic['mom'+str(mom)+'re'], renorm_dic['mom'+str(mom)+'im'] = self_renormalization.main(bare_dic['mom'+str(mom)]) # shape = (N(a_str), N_conf, len(z_ls_extend))

renorm_dic_avg = {}
for mom in [6, 8, 10]:
    renorm_dic_avg['mom'+str(mom)+'re'] = []
    renorm_dic_avg['mom'+str(mom)+'im'] = []
    for ida in range(3):
        renorm_dic_avg['mom'+str(mom)+'re'].append(gv.dataset.avg_data(renorm_dic['mom'+str(mom)+'re'][ida], bstrap=True))
        renorm_dic_avg['mom'+str(mom)+'im'].append(gv.dataset.avg_data(renorm_dic['mom'+str(mom)+'im'][ida], bstrap=True))

print(np.shape(renorm_dic_avg['mom6im']))


# print(bare_dic_avg) # mom, a, z
# print(renorm_dic_avg)

a_str_ls = [ '0.0574', '0.0882', '0.1213' ]

# f = open(meson+'_bare.txt', 'w')
# line = []
# line.append('mom')
# line.append('\t')
# line.append('a')
# line.append('\t')
# line.append('z')
# line.append('\t')
# line.append('Re')
# line.append('\t')
# line.append('Re_err')
# line.append('\t')
# line.append('Im')
# line.append('\t')
# line.append('Im_err')
# line.append('\n')
# f.writelines(line)

# for mom in [6, 8, 10]:
#     for ida in range(3):
#         a_str = a_str_ls[ida]
#         for idz in range(25):
#             line = []
#             line.append(str(mom))
#             line.append('\t')
#             line.append(a_str)
#             line.append('\t')
#             line.append(str(z_ls_extend[idz]))
#             line.append('\t')
#             line.append(str(bare_dic_avg['mom'+str(mom)+'re'][ida][idz].mean))
#             line.append('\t')
#             line.append(str(bare_dic_avg['mom'+str(mom)+'re'][ida][idz].sdev))
#             line.append('\t')
#             line.append(str(bare_dic_avg['mom'+str(mom)+'im'][ida][idz].mean))
#             line.append('\t')
#             line.append(str(bare_dic_avg['mom'+str(mom)+'im'][ida][idz].sdev))
#             line.append('\n')
#             f.writelines(line)
# f.close

            
# f = open(meson+'_renorm.txt', 'w')
# line = []
# line.append('mom')
# line.append('\t')
# line.append('a')
# line.append('\t')
# line.append('z')
# line.append('\t')
# line.append('Re')
# line.append('\t')
# line.append('Re_err')
# line.append('\t')
# line.append('Im')
# line.append('\t')
# line.append('Im_err')
# line.append('\n')
# f.writelines(line)

# for mom in [6, 8, 10]:
#     for ida in range(3):
#         a_str = a_str_ls[ida]
#         for idz in range(25):
#             line = []
#             line.append(str(mom))
#             line.append('\t')
#             line.append(a_str)
#             line.append('\t')
#             line.append(str(z_ls_extend[idz]))
#             line.append('\t')
#             line.append(str(renorm_dic_avg['mom'+str(mom)+'re'][ida][idz].mean))
#             line.append('\t')
#             line.append(str(renorm_dic_avg['mom'+str(mom)+'re'][ida][idz].sdev))
#             line.append('\t')
#             line.append(str(renorm_dic_avg['mom'+str(mom)+'im'][ida][idz].mean))
#             line.append('\t')
#             line.append(str(renorm_dic_avg['mom'+str(mom)+'im'][ida][idz].sdev))
#             line.append('\n')
#             f.writelines(line)
# f.close






# %%
