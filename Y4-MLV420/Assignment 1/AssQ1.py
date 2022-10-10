from ElementaryFlows import *
import numpy as np

n = 1000
xs = np.linspace(-2,2,n)
ys = np.linspace(-2,2,n)

flow = PotentialFlow(xs,ys)
flow.AddDoubletFlow(10)
flow.SolveFlow()
flow.PlotFlow(StreamFunction=True,StreamLines=False,LevelLabels=True,PotentialFunction=False,Grid=True)


