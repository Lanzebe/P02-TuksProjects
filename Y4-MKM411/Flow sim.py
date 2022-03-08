from queue import PriorityQueue
import numpy as np

class TwoDAirElement(object):
    def __init__(self,x1,x2,x3,x4,y1,y2,y3,y4):
        self.Xs = [x1,x2,x3,x4]
        self.Ys = [y1,y2,y3,y4]

        self.CG = [np.average(self.Xs),np.average(self.Ys)]

        self.Volume = 0.5*np.abs(np.dot(self.Xs,np.roll(self.Ys,1))-np.dot(self.Ys,np.roll(self.Xs,1)))

        self.VelocityVector = [0,0,0]

        self.BCVelocity = ['None','None','None','None']




xs = np.linspace(0,10,5)
ys = np.linspace(0,5,4)

Elements =  [0] * ((len(ys)-1)*(len(xs)-1)) 

a = 1

print(np.shape(Elements))


for i in range(len(Elements)):
    idx = i%(len(xs)-1)
    jdx = i//(len(xs)-1)

    #print(idx,jdx)

    Elements[i] = TwoDAirElement(xs[idx],xs[idx+1],xs[idx+1],xs[idx],ys[jdx],ys[jdx],ys[jdx+1],ys[jdx+1])
    if idx == 0:
        Elements[i].BCVelocity[2] = 0
    if idx == len(xs)-2:
        Elements[i].BCVelocity[0] = 0
    if jdx == 0:
        Elements[i].BCVelocity[3] = 0
    if jdx == len(ys)-2:
        Elements[i].BCVelocity[1] = 0
    #print(idx,jdx)
    #print(Elements[i].CG,Elements[i].BCVelocity)

        

        
#### inlett and outlett
Elements[1].BCVelocity = ['None','None','None',1]
Elements[-2].BCVelocity = ['None','None','None','None']


""" for i in Elements:
    print(i.CG,i.BCVelocity) """
