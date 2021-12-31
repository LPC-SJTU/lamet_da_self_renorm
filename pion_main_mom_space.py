# %%
from numpy import add
from head import *
from meson_da_hyb_selfrenorm_class import *
from pdf_self_renorm import pdf_zR
from plot import *
np.set_printoptions(threshold=np.inf)

def pion_main():
    meson = 'pion'
    zR_dic, m_pdf_dic = pdf_zR()

    x_ls = np.arange(-2-0.01, 3.02, 0.01) # x after ft, for quasi before matching
    y_ls = x_ls

    delta = 0.00001
    x_ls_mat = np.linspace(x_ls[0] + delta, x_ls[-1] - delta, 500)
    y_ls_mat = np.linspace(x_ls[0], x_ls[-1], 500)

    extend_point = {}
    extend_point['mom=6'] = -5#-4
    extend_point['mom=8'] = -9#-7
    extend_point['mom=10'] = -12#-10
    extend_point['mom=12'] = -10


    extend_fit_start = {}
    extend_fit_start['mom=6'] = -11
    extend_fit_start['mom=8'] = -14#-13
    extend_fit_start['mom=10'] = -16#-16
    extend_fit_start['mom=12'] = -16

    
    t_dic = {}
    t_dic['mom=6'] = {} 
    t_dic['mom=6']['a06'] = [12, 13, 14] # 1 means t=1
    t_dic['mom=6']['a09'] = [8, 9, 10]
    t_dic['mom=6']['a12'] = [6, 7, 8]
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

    file_path = 'DA_new.hdf5'
    data_pre = DATA_PRE(file_path, meson, const_fit=False)
    self_renormalization = RENORMALIZATION(zR_dic) # zR_dic 
    continuum_limit = CONTINUUM_LIMIT(['a06', 'a09', 'a12'], rotate=True)
    extrapolation_ft = EXTRAPOLATION_FT(meson, x_ls)
    inv_matching_mom = INV_MATCHING_MOM(x_ls, x_ls_mat, y_ls_mat)
    

    mom_ls = [6,8,10]
    lc_mom_mix = []
    quasi_mom_mix = []
    for mom in mom_ls:
        pz = mom * mom_to_pz

        if False:
            ### read data and normalize ###
            da_an_dic = data_pre.main(mom, t_dic['mom='+str(mom)]) # keys = ['a06', 'a09', 'a12'], list with this order will be used below
            da_an_ls = [da_an_dic[key] for key in da_an_dic] # ['a06', 'a09', 'a12'], shape = (N(a_str), N_conf, len(z_ls_extend)) (complex number)
            gv.dump(da_an_ls, meson+'/mom='+str(mom)+'/da_an_ls')

        da_an_ls = gv.load(meson+'/mom='+str(mom)+'/da_an_ls')


        if False:
            ### renormalization ###
            renorm_da_re_ls, renorm_da_im_ls = self_renormalization.main(da_an_ls) # shape = (N(a_str), N_conf, len(z_ls_extend))

            ### continuum_limit with rotate ###
            re_lam_ls, im_lam_ls, lam_ls, quasi_re_ls, quasi_im_ls = continuum_limit.main(mom, renorm_da_re_ls, renorm_da_im_ls) # quasi in the coordinate space
            gv.dump(re_lam_ls, meson+'/mom='+str(mom)+'/re_lam_ls')
            gv.dump(im_lam_ls, meson+'/mom='+str(mom)+'/im_lam_ls')
            gv.dump(lam_ls, meson+'/mom='+str(mom)+'/lam_ls')
            gv.dump(quasi_re_ls, meson+'/mom='+str(mom)+'/quasi_re_ls')
            gv.dump(quasi_im_ls, meson+'/mom='+str(mom)+'/quasi_im_ls')

        re_lam_ls = gv.load(meson+'/mom='+str(mom)+'/re_lam_ls')
        im_lam_ls = gv.load(meson+'/mom='+str(mom)+'/im_lam_ls')
        lam_ls = gv.load(meson+'/mom='+str(mom)+'/lam_ls')
        quasi_re_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_re_ls')
        quasi_im_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_im_ls')



        print('>>> mom = '+str(mom))
        print( 'For lambda bigger than ' + str(lam_ls[extend_point['mom='+str(mom)]]) + ', use extrapolation function.' )
        print( 'Lambda starts from ' + str(lam_ls[extend_fit_start['mom='+str(mom)]]) + ' are included in the extrapolation fit.' )

        # extrapolation_check(meson+', mom='+str(mom), lam_ls, quasi_re_ls, lam_ls[extend_fit_start['mom='+str(mom)]], lam_ls[extend_point['mom='+str(mom)]])

        

        if False:
            ### extrapolation and FT of quasi ###
            print('>>> extrapolation and FT of quasi at mom='+str(mom))
            lam_ls_ex, quasi_ext_ls, quasi_mom_ls = extrapolation_ft.main(lam_ls, extend_fit_start['mom='+str(mom)]+2, extend_point['mom='+str(mom)]+2, quasi_re_ls, quasi_im_ls, mode='lcda')
            gv.dump(lam_ls_ex, meson+'/mom='+str(mom)+'/lam_ls_ex')
            gv.dump(quasi_ext_ls, meson+'/mom='+str(mom)+'/quasi_ext_ls')
            gv.dump(quasi_mom_ls, meson+'/mom='+str(mom)+'/quasi_mom_ls')

        lam_ls_ex = gv.load(meson+'/mom='+str(mom)+'/lam_ls_ex')
        quasi_ext_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_ext_ls')
        quasi_mom_ls = gv.load(meson+'/mom='+str(mom)+'/quasi_mom_ls')

        quasi_mom_avg = gv.dataset.avg_data(quasi_mom_ls, bstrap=True)
        quasi_mom_mix.append(quasi_mom_avg)


        if False:
            ### inverse matching ###
            lc_mom_ls = inv_matching_mom.main(pz, quasi_mom_ls)
            gv.dump(lc_mom_ls, meson+'/mom='+str(mom)+'/lc_mom_ls_mom_mat')


        lc_mom_ls = gv.load(meson+'/mom='+str(mom)+'/lc_mom_ls_mom_mat') # in the mom space

        lc_mom_avg = gv.dataset.avg_data(lc_mom_ls, bstrap=True)

        

        

        # ### extrapolation at the endpoints ###
        # lc_mom_gv = [add_sdev(lc, lc_mom_avg) for lc in lc_mom_ls]

        # lc_mom_ls = []
        # print('>>> fitting the lc endpoints of '+meson)
        # for n_conf in tqdm(range(len(lc_mom_gv))):
        #     y_ls, lc_new = endpoint_ext(x_ls, lc_mom_gv[n_conf], meson)
        #     lc_mom_ls.append(lc_new)
        # lc_mom_avg = gv.dataset.avg_data(lc_mom_ls, bstrap=True)




        lc_mom_mix.append(lc_mom_avg)

        #quasi_vs_lc_plot(x_ls, y_ls, quasi_mom_avg, lc_mom_avg, pz, meson)


    large_mom_da = large_mom_limit(y_ls, lc_mom_mix, mom_ls)

    print('>>> large mom limit a2:')
    a2 = calc_an(y_ls, large_mom_da, 2)
    print(a2)

    print('>>> large mom limit a4:')
    a4 = calc_an(y_ls, large_mom_da, 4)
    print(a4)

    mellin_moment(y_ls, large_mom_da, 2)
    mellin_moment(y_ls, large_mom_da, 4)

    lcda_large_pz_plot(meson, y_ls, lc_mom_mix[-1], large_mom_da)
    # lcda_mix_pz_plot(meson, y_ls)
    

if __name__ == '__main__':
    pion_main()

# %%
