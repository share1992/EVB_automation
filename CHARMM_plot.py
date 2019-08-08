import pandas as pd
import matplotlib.pyplot as plt

# path = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating/100_sqa_52A_nvt_DYNA.txt'
# path_a = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating_continued/100_sqa_52A_nvt_cont_DYNA.txt'
# path_b = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating_continued_2/100_sqa_52A_nvt_cont_2_DYNA.txt'
# path_c = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating_continued_3/100_sqa_52A_nvt_cont_3_DYNA.txt'
# path_d = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating_continued_4/100_sqa_52A_nvt_cont_4_DYNA.txt'
# path_e = '~/Documents/CHAMPS/Squalane_project/52Abox_nve_post_nvt/100_sqa_52A_nve_DYNA.txt'
# path_f = '~/Documents/CHAMPS/Squalane_project/52Abox_nve_post_nvt_cont/100_sqa_52A_nve_cont_DYNA.txt'

# path = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_DYNA.txt'
# path_a = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_2_DYNA.txt'
# path_b = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_3_DYNA.txt'
# path_c = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_4_DYNA.txt'
# path_d = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_5_DYNA.txt'
# path_e = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_6_DYNA.txt'
# path_f = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_7_DYNA.txt'
# path_g = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_8_DYNA.txt'
# path_h = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_9_ext_DYNA.txt'
# path_i = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_10_ext_DYNA.txt'
# path_j = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_11_ext_DYNA.txt'
# path_k = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_12_ext_DYNA.txt'
# path_l = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_13_ext_DYNA.txt'
# path_m = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_14_ext_DYNA.txt'
# path_n = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_15_ext_DYNA.txt'
# path_o = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_16_ext_DYNA.txt'
# path_p = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_44p3A_nvt_DAVES_17_ext_DYNA_a.txt'

path = '~/Documents/CHAMPS/Squalane_project/DYNA_txts/100_sqa_1cn_nve_vv2_2_ext_DYNA.txt'


# path_d = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating_continued_4/100_sqa_52A_nvt_cont_4_DYNA.txt'
# path_e = '~/Documents/CHAMPS/Squalane_project/52Abox_nve_post_nvt/100_sqa_52A_nve_DYNA.txt'
# path_f = '~/Documents/CHAMPS/Squalane_project/52Abox_nve_post_nvt_cont/100_sqa_52A_nve_cont_DYNA.txt'


# NVT
dyna_data = pd.read_csv(path, delim_whitespace=True, header=None)
dyna_data.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_a = pd.read_csv(path_a, delim_whitespace=True, header=None)
# dyna_data_a.iloc[:,1] = dyna_data_a.iloc[:,1] + dyna_data['Time'].iloc[-1]
# dyna_data_a.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_b = pd.read_csv(path_b, delim_whitespace=True, header=None)
# dyna_data_b.iloc[:,1] = dyna_data_b.iloc[:,1] + dyna_data_a['Time'].iloc[-1]
# dyna_data_b.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_c = pd.read_csv(path_c, delim_whitespace=True, header=None)
# dyna_data_c.iloc[:,1] = dyna_data_c.iloc[:,1] + dyna_data_b['Time'].iloc[-1]
# dyna_data_c.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_d = pd.read_csv(path_d, delim_whitespace=True, header=None)
# dyna_data_d.iloc[:,1] = dyna_data_d.iloc[:,1] + dyna_data_c['Time'].iloc[-1]
# dyna_data_d.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_e = pd.read_csv(path_e, delim_whitespace=True, header=None)
# dyna_data_e.iloc[:,1] = dyna_data_e.iloc[:,1] + dyna_data_d['Time'].iloc[-1]
# dyna_data_e.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_f = pd.read_csv(path_f, delim_whitespace=True, header=None)
# dyna_data_f.iloc[:,1] = dyna_data_f.iloc[:,1] + dyna_data_e['Time'].iloc[-1]
# dyna_data_f.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_g = pd.read_csv(path_g, delim_whitespace=True, header=None)
# dyna_data_g.iloc[:,1] = dyna_data_g.iloc[:,1] + dyna_data_f['Time'].iloc[-1]
# dyna_data_g.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']


# EXPANDED BOX
# dyna_data_h = pd.read_csv(path_h, delim_whitespace=True, header=None)
# dyna_data_h.iloc[:,1] = dyna_data_h.iloc[:,1] + dyna_data_g['Time'].iloc[-1]
# dyna_data_h.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_i = pd.read_csv(path_i, delim_whitespace=True, header=None)
# dyna_data_i.iloc[:,1] = dyna_data_i.iloc[:,1] + dyna_data_h['Time'].iloc[-1]
# dyna_data_i.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_j = pd.read_csv(path_j, delim_whitespace=True, header=None)
# dyna_data_j.iloc[:,1] = dyna_data_j.iloc[:,1] + dyna_data_i['Time'].iloc[-1]
# dyna_data_j.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_k = pd.read_csv(path_k, delim_whitespace=True, header=None)
# dyna_data_k.iloc[:,1] = dyna_data_k.iloc[:,1] + dyna_data_j['Time'].iloc[-1]
# dyna_data_k.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_l = pd.read_csv(path_l, delim_whitespace=True, header=None)
# dyna_data_l.iloc[:,1] = dyna_data_l.iloc[:,1] + dyna_data_k['Time'].iloc[-1]
# dyna_data_l.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_m = pd.read_csv(path_m, delim_whitespace=True, header=None)
# dyna_data_m.iloc[:,1] = dyna_data_m.iloc[:,1] + dyna_data_l['Time'].iloc[-1]
# dyna_data_m.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_n = pd.read_csv(path_n, delim_whitespace=True, header=None)
# dyna_data_n.iloc[:,1] = dyna_data_n.iloc[:,1] + dyna_data_m['Time'].iloc[-1]
# dyna_data_n.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_o = pd.read_csv(path_o, delim_whitespace=True, header=None)
# dyna_data_o.iloc[:,1] = dyna_data_o.iloc[:,1] + dyna_data_n['Time'].iloc[-1]
# dyna_data_o.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
# dyna_data_p = pd.read_csv(path_p, delim_whitespace=True, header=None)
# dyna_data_p.iloc[:,1] = dyna_data_p.iloc[:,1] + dyna_data_o['Time'].iloc[-1]
# dyna_data_p.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']

# dyna_data = dyna_data_a.drop(dyna_data_a.index[0:10], axis=0)
# dyna_data = dyna_data.append(dyna_data_b)
# dyna_data = dyna_data.append(dyna_data_c)
# dyna_data = dyna_data.append(dyna_data_d)
# dyna_data = dyna_data.append(dyna_data_e)
# dyna_data = dyna_data.append(dyna_data_f)
# dyna_data = dyna_data.append(dyna_data_g)


# #
# dyna_data_exp = dyna_data_h.drop(dyna_data_h.index[0:10], axis=0)
# dyna_data_exp = dyna_data_exp.append(dyna_data_i)
# dyna_data_exp = dyna_data_exp.append(dyna_data_j)
# dyna_data_exp = dyna_data_exp.append(dyna_data_k)
# dyna_data_exp = dyna_data_exp.append(dyna_data_l)
# dyna_data_exp = dyna_data_exp.append(dyna_data_m)
# dyna_data_exp = dyna_data_exp.append(dyna_data_n)
# dyna_data_exp = dyna_data_exp.append(dyna_data_o)
# dyna_data_exp = dyna_data_exp.append(dyna_data_p)


first_energy = dyna_data['TOTEner'].iloc[0]
dyna_data['Relative Energy'] = dyna_data['TOTEner'] - first_energy
# dyna_data_exp['Relative Energy'] = dyna_data_exp['TOTEner'] - first_energy
# dyna_data_exp['Time'] = dyna_data_exp['Time'] + dyna_data['Time'].iloc[-1]

fig = plt.figure(figsize=(8,8))

ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
# ax3 = plt.subplot(2,2,3)
# ax4 = plt.subplot(2,2,4)

ax1.plot(dyna_data['Time'], dyna_data['TEMPerature'], color='C0')
# ax3.plot(dyna_data_exp['Time'], dyna_data_exp['TEMPerature'], color='C4')

ax2.plot(dyna_data['Time'], dyna_data['Relative Energy'], color='C0')
# ax4.plot(dyna_data_exp['Time'], dyna_data_exp['Relative Energy'], color='C4')

ax1.set_xlabel('Time (ps)')
ax1.set_ylabel('Temperature (K)')
ax1.set_title('Temperature Equilibration')
ax2.set_xlabel('Time (ps)')
ax2.set_ylabel('Relative Energy (kcal/mol)')
ax2.set_title('Energy Equilibration')
# ax3.set_xlabel('Time (ps)')
# ax3.set_ylabel('Temperature (K)')
# ax3.set_title('Temperature Equilibration')
# ax4.set_xlabel('Time (ps)')
# ax4.set_ylabel('Relative Energy (kcal/mol)')
# ax4.set_title('Energy Equilibration')

smooth_data_T = dyna_data['TEMPerature'].rolling(50).mean()
smooth_data_E = dyna_data['Relative Energy'].rolling(50).mean()

# smooth_data_T_exp = dyna_data_exp['TEMPerature'].rolling(50).mean()
# smooth_data_E_exp = dyna_data_exp['Relative Energy'].rolling(50).mean()

ax1.plot(dyna_data['Time'], smooth_data_T, color='C1')
# ax3.plot(dyna_data_exp['Time'], smooth_data_T_exp, color='yellow')

ax2.plot(dyna_data['Time'], smooth_data_E, color='C1')
# ax4.plot(dyna_data_exp['Time'], smooth_data_E_exp, color='yellow')

# ax1.axvline(1000, color='k', linestyle=':')
# ax2.axvline(1000, color='k', linestyle=':')

plt.tight_layout()
plt.savefig("/Users/ec18006/Documents/CHAMPS/Squalane_project/CHARMM_equil_plot_cube_slab.png")
plt.show()

