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
flow.PlotFlow(PlotInstantly = False, StreamFunction=False,StreamLines=True,LevelLabels=True,PotentialFunction=False,Grid=True, VelocityLevels=[Lambda/2])

xs = np.linspace(0,1,100)

theta1 = np.pi/12
theta2 = -np.pi/12 + 2*np.pi

A = [[0,0,0,0,1],[1,1,1,1,1],[4,3,2,1,0],[0,0,0,1,0],[12*((0.5)**2),6*0.5, 2 ,0 ,0]]
B = [theta1,theta2,0,0,0]

a,b,c,d,e = np.linalg.solve(A,B)

Thetas = a*xs**4 + b*xs**3 + c*xs**2 + d*xs + e

Radii = (Lambda/(2*V_inf*np.sin(Thetas)))*(1-(Thetas/np.pi))

Xss = np.multiply(Radii,np.cos(Thetas)) + flow.ElementaryFlows[1].Position[0]
Yss = np.multiply(Radii,np.sin(Thetas)) + flow.ElementaryFlows[1].Position[1]


PressC = flow.PressureCoefficientAt(Xss,Yss)

scale = 3.5
headlenght = 0.04

dXss = np.multiply(scale,np.multiply(np.abs(PressC),np.cos(Thetas)))
dYss = np.multiply(scale,np.multiply(np.abs(PressC),np.sin(Thetas)))

for i in range(len(Xss)):
    dx = scale*abs(PressC[i])*np.cos(Thetas[i])
    dy = scale*abs(PressC[i])*np.sin(Thetas[i])
    if PressC[i] <0:
            
        plt.arrow(Xss[i],Yss[i], dx,dy,color = 'blue',length_includes_head = True,head_length = headlenght, head_width = 0.5*headlenght)
    else:
        plt.arrow(Xss[i],Yss[i], dx,dy,color = 'red',length_includes_head = True,head_length = headlenght, head_width = 0.5*headlenght)

plt.plot(np.add(Xss,dXss),np.add(Yss,dYss),color = '0.3', linewidth = 0.5)

plt.show()