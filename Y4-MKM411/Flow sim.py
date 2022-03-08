from re import A
from xml.etree.ElementTree import ElementTree
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


class TwoDAirElement(object):
    def __init__(self,x1,x2,x3,x4,y1,y2,y3,y4):
        self.Xs = [x1,x2,x3,x4]
        self.Ys = [y1,y2,y3,y4]

        self.CG = [np.average(self.Xs),np.average(self.Ys)]

        self.Volume = 0.5*np.abs(np.dot(self.Xs,np.roll(self.Ys,1))-np.dot(self.Ys,np.roll(self.Xs,1)))

        self.VelocityVector = [0,0,0]




xs = np.linspace(0,10,5)
ys = np.linspace(0,5,4)

Elements = [ [0] * (len(ys)-1) ] * (len(xs)-1) 

a = 1


for idx, i in enumerate(Elements):
    for jdx, j in enumerate(i):
        Elements[idx][jdx] = TwoDAirElement(xs[idx],xs[idx+1],xs[idx+1],xs[idx],ys[jdx],ys[jdx],ys[jdx+1],ys[jdx+1])
        print(Elements[idx][jdx].CG,Elements[idx][jdx].Volume)

