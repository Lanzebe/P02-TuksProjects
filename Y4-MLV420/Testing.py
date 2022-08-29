from ElementaryFlows import *
import numpy as np

V_inf = 30
r = 0.3

Gamma_crit = 4*np.pi*V_inf*r

n = 500
xs = np.linspace(-2,2,n)
ys = np.linspace(-2,2,n)

flow1 = PotentialFlow(xs,ys)
flow1.AddUniformFlow(V_inf)
flow1.AddDoubletFlow(2*np.pi*V_inf*(r**2))
#flow1.AddSourceFlow(100)
#flow1.AddVortexFlow(Gamma_crit/2,RadiusZeroPoint = r)
flow1.SolveFlow()
""" 
for idx,i in enumerate(flow1.Vel):
    for jdx,j in enumerate(i):       
        if (flow1.Xs[idx,jdx]**2 + flow1.Ys[idx,jdx]**2)**0.5 < r:
            flow1.Vel[idx,jdx] = 0 """

flow1.PlotFlow(StreamFunction=True,LevelLabels=False,StreamLines=True,PotentialFunction=False,PlotInstantly=True)