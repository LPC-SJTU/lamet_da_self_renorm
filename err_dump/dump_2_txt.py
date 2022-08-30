# %%
import gvar as gv

f = open('./pion.txt', 'w')
line = []
line.append('x')
line.append('\t')
line.append('mean')
line.append('\t')
line.append('stat')
line.append('\t')
line.append('large_mom_limit')
line.append('\t')
line.append('continuum_limit')
line.append('\t')
line.append('extrapolation')
line.append('\t')
line.append('mu')
line.append('\n')
f.writelines(line)

x_ls = gv.load('x_ls_pion')
mean = gv.load('mean_pion')
stat = gv.load('stat_pion')
large_mom = gv.load('large_mom_limit_pion')
cont = gv.load('sys_continuum_limit_pion')
ext = gv.load('sys_extrapolation_pion')
mu = gv.load('sys_mu_pion')

for i in range(len(x_ls)):
    line = []
    line.append(str(x_ls[i]))
    line.append('\t')
    line.append(str(mean[i]))
    line.append('\t')
    line.append(str(stat[i]))
    line.append('\t')
    line.append(str(large_mom[i]))
    line.append('\t')
    line.append(str(cont[i]))
    line.append('\t')
    line.append(str(ext[i]))
    line.append('\t')
    line.append(str(mu[i]))
    line.append('\n')
    f.writelines(line)
f.close

# %%
f = open('./kaon.txt', 'w')
line = []
line.append('x')
line.append('\t')
line.append('mean')
line.append('\t')
line.append('stat')
line.append('\t')
line.append('large_mom_limit')
line.append('\t')
line.append('continuum_limit')
line.append('\t')
line.append('extrapolation')
line.append('\t')
line.append('mu')
line.append('\n')
f.writelines(line)

x_ls = gv.load('x_ls_kaon')
mean = gv.load('mean_kaon')
stat = gv.load('stat_kaon')
large_mom = gv.load('large_mom_limit_kaon')
cont = gv.load('sys_continuum_limit_kaon')
ext = gv.load('sys_extrapolation_kaon')
mu = gv.load('sys_mu_kaon')

for i in range(len(x_ls)):
    line = []
    line.append(str(x_ls[i]))
    line.append('\t')
    line.append(str(mean[i]))
    line.append('\t')
    line.append(str(stat[i]))
    line.append('\t')
    line.append(str(large_mom[i]))
    line.append('\t')
    line.append(str(cont[i]))
    line.append('\t')
    line.append(str(ext[i]))
    line.append('\t')
    line.append(str(mu[i]))
    line.append('\n')
    f.writelines(line)
f.close
# %%
