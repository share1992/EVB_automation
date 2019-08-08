import pandas as pd
import matplotlib.pyplot as plt

path = '~/Documents/CHAMPS/Squalane_project/52Abox_nvt_post_heating/100_sqa_52A_nvt_DYNA.txt'

dyna_data = pd.read_csv(path, delim_whitespace=True, header=None)
dyna_data.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']
first_energy = dyna_data['TOTEner'][0]
dyna_data['Relative Energy'] = dyna_data['TOTEner'] - first_energy

ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)
ax1.plot(dyna_data['Step']/1000, dyna_data['TEMPerature'])
ax2.plot(dyna_data['Step']/1000, dyna_data['Relative Energy'])

ax1.set_xlabel('Time (ps)')
ax1.set_ylabel('Temperature (K)')
ax1.set_title('Temperature Equilibration')
ax2.set_xlabel('Time (ps)')
ax2.set_ylabel('Relative Energy (kcal/mol)')
ax2.set_title('Energy Equilibration')

smooth_data_T = dyna_data['TEMPerature'].rolling(50).mean()
smooth_data_E = dyna_data['Relative Energy'].rolling(50).mean()

ax1.plot(dyna_data['Step']/1000, smooth_data_T)
ax2.plot(dyna_data['Step']/1000, smooth_data_E)


plt.tight_layout()
plt.show()
