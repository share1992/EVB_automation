import pandas as pd
import matplotlib.pyplot as plt

def plot_T_E(path):

    dyna_data = pd.read_csv(path, delim_whitespace=True, header=None)
    dyna_data.columns = ['DYNA', 'Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature']

    first_energy = dyna_data['TOTEner'].iloc[0]
    dyna_data['Relative Energy'] = dyna_data['TOTEner'] - first_energy
    new_dyna_data = dyna_data.tail(20000)

    fig = plt.figure(figsize=(8, 8))

    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    ax1.plot(dyna_data['Time'], dyna_data['TEMPerature'], color='C0')
    # ax2.plot(dyna_data['Time'], dyna_data['Relative Energy'], color='C0')
    ax2.plot(dyna_data['Time'], dyna_data['TOTEner'], color='C0')
    ax3.plot(new_dyna_data['Time'], new_dyna_data['TEMPerature'], color='C0')
    ax4.plot(new_dyna_data['Time'], new_dyna_data['Relative Energy'], color='C0')

    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Temperature')
    ax2.set_xlabel('Time (ps)')
    # ax2.set_ylabel('Relative Energy (kcal/mol)')
    # ax2.set_title('Relative Total Energy')
    ax2.set_ylabel('Total Energy (kcal/mol)')
    ax2.set_title('Total Energy')

    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Temperature (K)')
    ax3.set_title('Temperature')

    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Relative Energy (kcal/mol)')
    ax4.set_title('Relative Energy')

    smooth_data_T = dyna_data['TEMPerature'].rolling(100).mean()
    smooth_data_T2 = new_dyna_data['TEMPerature'].rolling(100).mean()
    # smooth_data_E = dyna_data['Relative Energy'].rolling(50).mean()
    smooth_data_E = dyna_data['TOTEner'].rolling(100).mean()
    smooth_data_E2 = new_dyna_data['Relative Energy'].rolling(100).mean()


    ax1.plot(dyna_data['Time'], smooth_data_T, color='C1')
    ax2.plot(dyna_data['Time'], smooth_data_E, color='C1')
    ax3.plot(new_dyna_data['Time'], smooth_data_T2, color='C1')
    ax4.plot(new_dyna_data['Time'], smooth_data_E2, color='C1')

    plt.tight_layout()
    # plt.savefig("CHARMM_equil_plot_sqa_cn.png")
    plt.show()

    return dyna_data, new_dyna_data

