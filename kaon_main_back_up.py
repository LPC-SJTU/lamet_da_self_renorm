# %%
from scipy.sparse import construct
from head import *
from pdf_self_renorm import pdf_zR
from meson_da_hyb_class import MESON_DA_HYB
from plot import *
np.set_printoptions(threshold=np.inf)

meson = 'kaon'

def kaon_main():
    zR_dic, m_pdf_dic = pdf_zR()

    x_ls = np.arange(-2-0.01, 3.02, 0.01) # x after ft, for quasi before matching

    delta = 0.01

    ls_1 = np.arange(-2-0.000001, 0.5-delta, delta)
    ls_2 = np.array([0.5])
    ls_3 = np.arange(0.5+delta+0.000001, 3+delta, delta)
    # ls_1 = np.arange(-2+0.0000001, 0.5, delta)
    # ls_2 = np.array([0.5])
    # ls_3 = np.arange(0.5+delta-0.0000001, 3, delta)

    x_ls_matching = np.hstack((ls_1, ls_2, ls_3)) # x for matching, quasi

    ls_4 = np.arange(-2+0.0000001, 0.5, delta)
    ls_5 = np.array([0.5])
    ls_6 = np.arange(0.5+delta-0.0000001, 3, delta)
    y_ls_matching = np.hstack((ls_4, ls_5, ls_6)) # x for matching, light-cone

    print(len(x_ls_matching))
    print(len(y_ls_matching))

    extend_point = {}
    extend_point['mom=6'] = -2
    extend_point['mom=8'] = -5
    extend_point['mom=10'] = -8
    extend_point['mom=12'] = -10


    extend_fit_start = {}
    extend_fit_start['mom=6'] = -8
    extend_fit_start['mom=8'] = -11
    extend_fit_start['mom=10'] = -13
    extend_fit_start['mom=12'] = -14
    
    t_dic = {}
    t_dic['mom=6'] = {}
    t_dic['mom=6']['a06'] = [12, 13, 14]
    t_dic['mom=6']['a09'] = [9, 10, 11]
    t_dic['mom=6']['a12'] = [5, 6, 7]
    t_dic['mom=8'] = {}
    t_dic['mom=8']['a06'] = [8, 9, 10]
    t_dic['mom=8']['a09'] = [6, 7, 8]
    t_dic['mom=8']['a12'] = [4, 5, 6]
    t_dic['mom=10'] = {}
    t_dic['mom=10']['a06'] = [7, 8, 9]
    t_dic['mom=10']['a09'] = [5, 6, 7]
    t_dic['mom=10']['a12'] = [3, 4, 5]
    t_dic['mom=12'] = {}
    t_dic['mom=12']['a06'] = [7, 8, 9]
    t_dic['mom=12']['a09'] = [5, 6, 7]
    t_dic['mom=12']['a12'] = [3, 4, 5]

    lcda_fit = True
    constant_fit = False

    ################################################
    lambda_ls = z_ls_da * 2*np.pi / (0.0574*96) * 6
    print('mom=6:')
    print( 'For lambda bigger than ' + str(lambda_ls[extend_point['mom=6']]) + ', use extrapolation function.' )
    print( 'Lambda starts from ' + str(lambda_ls[extend_fit_start['mom=6']]) + ' are included in the extrapolation fit.' )
    kaon_mom_6 = MESON_DA_HYB(meson, 6, x_ls, x_ls_matching, y_ls_matching, extend_point['mom=6'], extend_fit_start['mom=6'], t_dic['mom=6'], gs_extract=False, fit_1=False, fit_2=False, rotate=False, lcda_fit=lcda_fit, constant_fit=constant_fit)
    mom_6_kernel, mom_6_quasi_da, mom_6_y_ls, mom_6_lic_da = kaon_mom_6.main(zR_dic)
    gv.dump(mom_6_lic_da, 'kaon/mom=6/mom_6_lic_da')
    quasi_vs_lc_plot(x_ls, mom_6_y_ls, mom_6_quasi_da, mom_6_lic_da, pz=1.29, meson='kaon')
    y_ls = mom_6_y_ls

    ################################################
    lambda_ls = z_ls_da * 2*np.pi / (0.0574*96) * 8
    print('mom=8:')
    print( 'For lambda bigger than ' + str(lambda_ls[extend_point['mom=8']]) + ', use extrapolation function.' )
    print( 'Lambda starts from ' + str(lambda_ls[extend_fit_start['mom=8']]) + ' are included in the extrapolation fit.' )
    kaon_mom_8 = MESON_DA_HYB(meson, 8, x_ls, x_ls_matching, y_ls_matching, extend_point['mom=8'], extend_fit_start['mom=8'], t_dic['mom=8'], gs_extract=False, fit_1=False, fit_2=False, rotate=False, lcda_fit=lcda_fit, constant_fit=constant_fit)
    mom_8_kernel, mom_8_quasi_da, mom_8_y_ls, mom_8_lic_da = kaon_mom_8.main(zR_dic)
    gv.dump(mom_8_lic_da, 'kaon/mom=8/mom_8_lic_da')
    quasi_vs_lc_plot(x_ls, mom_8_y_ls, mom_8_quasi_da, mom_8_lic_da, pz=1.72, meson='kaon')
    y_ls = mom_8_y_ls

    ################################################
    lambda_ls = z_ls_da * 2*np.pi / (0.0574*96) * 10
    print('mom=10:')
    print( 'For lambda bigger than ' + str(lambda_ls[extend_point['mom=10']]) + ', use extrapolation function.' )
    print( 'Lambda starts from ' + str(lambda_ls[extend_fit_start['mom=10']]) + ' are included in the extrapolation fit.' )
    kaon_mom_10 = MESON_DA_HYB(meson, 10, x_ls, x_ls_matching, y_ls_matching, extend_point['mom=10'], extend_fit_start['mom=10'], t_dic['mom=10'], gs_extract=False, fit_1=False, fit_2=False, rotate=False, lcda_fit=lcda_fit, constant_fit=constant_fit)
    mom_10_kernel, mom_10_quasi_da, mom_10_y_ls, mom_10_lic_da = kaon_mom_10.main(zR_dic)
    gv.dump(mom_10_lic_da, 'kaon/mom=10/mom_10_lic_da')
    quasi_vs_lc_plot(x_ls, mom_10_y_ls, mom_10_quasi_da, mom_10_lic_da, pz=1.72, meson='kaon')
    y_ls = mom_10_y_ls

    ################################################
    lambda_ls = z_ls_da * 2*np.pi / (0.0574*96) * 12
    print('mom=12:')
    print( 'For lambda bigger than ' + str(lambda_ls[extend_point['mom=12']]) + ', use extrapolation function.' )
    print( 'Lambda starts from ' + str(lambda_ls[extend_fit_start['mom=12']]) + ' are included in the extrapolation fit.' )
    kaon_mom_12 = MESON_DA_HYB(meson, 12, x_ls, x_ls_matching, y_ls_matching, extend_point['mom=12'], extend_fit_start['mom=12'], t_dic['mom=12'], gs_extract=False, fit_1=False, fit_2=False, rotate=False, lcda_fit=lcda_fit, constant_fit=constant_fit)
    mom_12_kernel, mom_12_quasi_da, mom_12_y_ls, mom_12_lic_da = kaon_mom_12.main(zR_dic)
    gv.dump(mom_12_lic_da, 'kaon/mom=12/mom_12_lic_da')
    quasi_vs_lc_plot(x_ls, mom_12_y_ls, mom_12_quasi_da, mom_12_lic_da, pz=1.72, meson='kaon')
    y_ls = mom_12_y_ls


    # ################################################
    mom_6_lic_da = gv.load('kaon/mom=6/mom_6_lic_da')
    mom_8_lic_da = gv.load('kaon/mom=8/mom_8_lic_da')
    mom_10_lic_da = gv.load('kaon/mom=10/mom_10_lic_da')
    # mom_12_lic_da = gv.load('kaon/mom=12/mom_12_lic_da')

    ### large momentum limit ###
    ############################
    large_mom_lic_da = large_mom_limit(y_ls, mom_6_lic_da, mom_8_lic_da, mom_10_lic_da, mom_ls=[6, 8, 10], meson=meson)
    gv.dump(large_mom_lic_da, 'kaon/large_mom_lic_da')

    # paper_plot_discrete_effect(mom=6, meson=meson)
    # paper_plot_discrete_effect(mom=8, meson=meson)
    paper_plot_discrete_effect(mom=10, meson=meson)
    # paper_plot_discrete_effect(mom=12, meson=meson)

    # continuous_limit_pz_mix(meson, mom_ls=[6,8,10,12])
    lcda_mix_pz_plot(meson, y_ls)

    lcda_large_pz_plot(meson, y_ls, mom_10_lic_da, large_mom_lic_da)

if __name__ == '__main__':
    kaon_main()
# %%