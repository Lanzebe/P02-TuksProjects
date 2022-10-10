import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


RawDataSet = np.load('CleanedDataSet.npy', allow_pickle=True)



# Create Figure

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")

# Create Plot

scatter_plot = ax.scatter3D(RawDataSet[:,1],RawDataSet[:,2],RawDataSet[:,10], c= RawDataSet[:,5])
#ax.axes.set_zlim3d(bottom=-4, top=-1)
ax.axes.set_xlabel('l/t')
ax.axes.set_ylabel('b/t')
ax.axes.set_zlabel('log(omega x t x (rho/E)^0.5')
plt.colorbar(scatter_plot)
# Show plot

plt.show()