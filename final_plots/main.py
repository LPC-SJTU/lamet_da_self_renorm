# %%
import gvar as gv
from plot_func import *

#! Pion using momentum space matching

pion_mom_space = gv.load('dump/pion_final_plot_mom_space.pkl')

lcda_large_pz_plot('pion', pion_mom_space['y_ls'], pion_mom_space['mom_n_lic_da'], pion_mom_space['large_mom_lic_da'])


# %%
#! Kaon using momentum space matching

kaon_mom_space = gv.load('dump/kaon_final_plot_mom_space.pkl')

lcda_large_pz_plot('kaon', kaon_mom_space['y_ls'], kaon_mom_space['mom_n_lic_da'], kaon_mom_space['large_mom_lic_da'])


# %%
