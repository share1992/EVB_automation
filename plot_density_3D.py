import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import squalane_CN_contact as sqcn
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import cm


def plot_3d_hist(trajectory, slices=50, step=10, azim=330, elev=50):

    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(111, projection='3d')

    times = np.arange(0, len(trajectory), step)
    color = iter(cm.winter(np.linspace(0, 1, len(times))))

    for i in times:
        c = next(color)
        hist, bins = np.histogram(trajectory[i]['Z'], bins=slices)
        xs = (bins[:-1] + bins[1:]) / 2
        ax.bar(xs, hist, zs=i, zdir='y', color=c, alpha=0.6)
        ax.plot(xs, hist, zs=i, zdir='y', color=c, linewidth=3)

    ax.set_xlabel('Slice')
    ax.set_ylabel('Time')
    ax.set_zlabel('Number of Atoms (density)')
    ax.view_init(elev=elev, azim=azim)
    # ax.set_zlim(0, 100)
    # plt.savefig('density_slab_%s_jet.png' % azim)
    plt.show()

path = '/Users/ec18006/OneDrive - University of Bristol/CHAMPS/Research_Topics/Squalane_project/100SQA_equilibration/300K_slab_equil.pdb'
frames = sqcn.read_pdb(path)
plot_3d_hist(frames, slices=30, step=100)

# For gif making
# for n in range(0, 360, 15):
#     plot_3d_hist(frames, step=100, azim=n)

