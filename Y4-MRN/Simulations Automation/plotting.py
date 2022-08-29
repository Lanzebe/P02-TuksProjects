import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


RawDataSet = np.load('Dataset5.npy', allow_pickle=True)



OmegaND1 = np.log10(RawDataSet[:,6]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND2 = np.log10(RawDataSet[:,7]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND3 = np.log10(RawDataSet[:,8]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND4 = np.log10(RawDataSet[:,9]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND5 = np.log10(RawDataSet[:,10]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND6 = np.log10(RawDataSet[:,11]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND7 = np.log10(RawDataSet[:,12]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
OmegaND8 = np.log10(RawDataSet[:,13]*RawDataSet[:,0]*np.sqrt(RawDataSet[:,4]/RawDataSet[:,3]))
AR_A = RawDataSet[:,1]/RawDataSet[:,0]
AR_B = RawDataSet[:,2]/RawDataSet[:,0]
NU = RawDataSet[:,5]


# Create Figure

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
 
# Create Plot

scatter_plot = ax.scatter3D(AR_A,AR_B,OmegaND8, c= NU)
 

plt.colorbar(scatter_plot)
# Show plot

plt.show()