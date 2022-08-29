from ElementaryFlows import *
import numpy as np

n = 500
xs = np.linspace(-2,4,n)
ys = np.linspace(-2,2,n)

V_inf = 10
Lambda = 10

flow = PotentialFlow(xs,ys)
flow.AddUniformFlow(V_inf)
flow.AddSourceFlow(Lambda,1,0)
flow.SolveFlow()
flow.PlotFlow(PlotInstantly = False, PressureCoefficient=False, StreamFunction=False,StreamLines=True,LevelLabels=True,PotentialFunction=False,Grid=True, VelocityLevels=[Lambda/2])


vmin = 0
vmax = 10

for idx,i in enumerate(flow.PressureCoefficient):
    for jdx,j in enumerate(i):       
        if (j < vmin):
            flow.PressureCoefficient[idx,jdx] = 0 
plt.contourf(flow.Xs,flow.Ys,flow.PressureCoefficient)
plt.colorbar()
#plt.clim(vmin=vmin,vmax=vmax)
plt.show()