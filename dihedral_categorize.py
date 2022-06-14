from md_to_pd_tools import dihedral_into_pd
from identify_rotatable_bonds import identify_main_rotatable_bonds
import numpy as np
import scipy.stats as scistat
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import pandas as pd

struc = 'epctcng.pdb'
traj = 'crest_rotamers.xyz'

dihedral_df = dihedral_into_pd(identify_main_rotatable_bonds(struc), struc, traj)

#print(dihedral_df.tail())

col_1 = dihedral_df['dihedral_1']

#eval_points = np.linspace(np.min(col_1), np.max(col_1))
#kde = scistat.gaussian_kde(col_1,bw_method='scott')
#y = kde.pdf(eval_points)
#plt.plot(eval_points,y)
#plt.savefig('dihedral_1.kdescistat.png')
#plt.clf()

#model = KernelDensity()
#model.fit(col_1)
#log_dens = model.score_samples(eval_points)
#plt.fill(eval_points, np.exp(log_dens), c='cyan')
#plt.savefig('dihedral_1.scikit.png')
#plt.clf()

#scat = col_1.plot.scatter()
#figscat = scat.get_figure()
#figscat.savefig('dihedral_1.scatter.png')
#figscat.clf()

#ax = col_1.plot.kde()
#kdy = ax.get_lines()[0].get_xydata()

#mi, ma = argrelextrema(ax.get_lines()[0].get_xydata(), np.less)[0], argrelextrema(ax.get_lines()[0].get_xydata(), np.greater)[0]
#print(mi, ma)
#fig = ax.get_figure()
#fig.savefig('dihedral_1.kde.plot.png')
#fig.clf()

### DO KDE for regular dataset

#a = np.array(col_1).reshape(-1,1)
#kde = KernelDensity(kernel='gaussian', bandwidth=3).fit(a)
#s = np.linspace(np.min(col_1),np.max(col_1))
#e = kde.score_samples(s.reshape(-1,1))
#plt.plot(s,e)
#plt.savefig('dihedral_1.sklearn.test.png')
#plt.clf()

#mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
#min_index = np.where(e[mi]==e[mi].min())
#print(s[mi][min_index])

#print(col_1.loc[2.16679613 > col_1])

#sorted(s[mi], key=lambda x: e[x])

def apply_kde_to_col(col, plot=False, shift=False, clusters=False):
    categories = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    a = np.array(col).reshape(-1,1)
    kde = KernelDensity(kernel='gaussian', bandwidth=3).fit(a)
    s = np.linspace(np.min(col),np.max(col))
    e = kde.score_samples(s.reshape(-1,1))
    if plot:
        plt.plot(s,e)
        title = col.name + '.kde.plot.png'
        plt.savefig(title)
        plt.clf()
    if shift:
        mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
        min_index = np.where(e[mi]==e[mi].min())
        new_180 = s[mi][min_index]
        pos_shift = 180 - new_180
        neg_shift = new_180 - 180
        shifted_col = col.map(lambda x: (x + pos_shift) if (x + pos_shift)<180 else (x - pos_shift))
        shifted_a = np.array(shifted_col).reshape(-1,1)
        shifted_kde = KernelDensity(kernel='gaussian', bandwidth=3).fit(shifted_a)
        shifted_s = np.linspace(np.min(shifted_col),np.max(shifted_col))
        shifted_e = shifted_kde.score_samples(shifted_s.reshape(-1,1))
        #plt.plot(shifted_s,shifted_e)
        #plt.savefig('test_shifted_7.png')
        #plt.clf()
        shifted_mi, shifted_ma = argrelextrema(shifted_e, np.less)[0], argrelextrema(shifted_e, np.greater)[0]
    if clusters:
        mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
        #print(s[mi])
        bins = [-1000]+list(s[mi])+[1000]
        labs = categories[:len(bins)-2] + [categories[0]]
        categorized_col = pd.cut(col,bins=bins,right=False,labels=labs,ordered=False)
