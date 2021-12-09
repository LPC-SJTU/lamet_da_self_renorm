# %%

from head import *

def val(px):
    return 0.2 + (px - 495)/(844 - 495) * 0.1


fig = plt.figure(figsize=fig_size_lc)
plt_axes_ = [0.25,0.15,0.7,0.8]
ax = plt.axes(plt_axes_)

ax.errorbar([val(790)], [1], xerr=[val(980) - val(790)], fmt='D', **errorb)
ax.errorbar([val(774)], [2], xerr=[val(844) - val(774)], fmt='D', **errorb)
ax.errorbar([val(620)], [3], xerr=[val(646) - val(620)], fmt='D', **errorb)
ax.errorbar([val(648)], [4], xerr=[val(754) - val(648)], fmt='D', **errorb)
ax.errorbar([val(531)], [5], xerr=[val(576) - val(531)], fmt='D', color='brown', **errorb)
ax.errorbar([val(531)], [5], xerr=[val(691) - val(531)], fmt='D', color='brown', **errorb)
ax.errorbar([val(531)], [5.5], xerr=[val(664) - val(531)], color='violet', fmt='D', **errorb)
ax.errorbar([val(531)], [5.5], xerr=[val(576) - val(531)], fmt='D', color='violet', **errorb)
ax.errorbar([0.303], [6.5], xerr=[0.029], fmt='D', color='purple', **errorb)


ax.axvline(0.303-0.029, color='purple', linestyle='--', lw=0.5)
ax.axvline(0.303+0.029, color='purple', linestyle='--', lw=0.5)
ax.axvline(0.303-2*0.029, color='purple', linestyle='dashdot', lw=0.5)

ax.set_xlabel(x_label, **fs_p)
ax.set_ylim([0.5, 7])
ax.set_xlim([0.15, 0.35])
ax.tick_params(direction='in', labelsize=12)
ax.set_xticks([0.2, 0.3])
ax.set_xlabel(r'$<\xi^2>(\mu = 2 \rm{GeV})$')
plt.yticks([1,2,3,4,5,5.5, 6.5], ['Del Debbio et al.\n(2003)', 'Arthur et al.\n(2011)', 'Bali et al.\n(2019)', 'Zhang et al.\n(2020)', 'HOPE Mom', 'HOPE TMR', 'This work'], fontsize=10)
plt.show()