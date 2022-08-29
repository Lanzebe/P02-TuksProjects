from ElementaryFlows import *
import numpy as np
import matplotlib.pyplot as plt

V_inf = 30
r = 0.3

Gamma_crit = 4*np.pi*V_inf*r
Gamma_1 = Gamma_crit/2
Gamma_3 = Gamma_crit*2

n = 800
xs = np.linspace(-2,4,n)
ys = np.linspace(-2,2,n)

flow1 = PotentialFlow(xs,ys)
flow1.AddUniformFlow(V_inf)
flow1.AddDoubletFlow(2*np.pi*V_inf*(r**2))
flow1.AddVortexFlow(Gamma_1,RadiusZeroPoint = r)
flow1.SolveFlow()

### Calculation of Stagnation Points
flowThetaPoint1 = [np.arcsin(-Gamma_1/(4*np.pi*V_inf*r)),-np.pi-np.arcsin(-Gamma_1/(4*np.pi*V_inf*r))]
flowRadiiPoint1 = [r,r]

StagnationPoints1X = flowRadiiPoint1*np.cos(flowThetaPoint1)
StagnationPoints1Y = flowRadiiPoint1*np.sin(flowThetaPoint1)

flow1.PlotFlow(PotentialFunction=False,PlotInstantly=False,VelocityMagnitude=False,LevelLabels=True)
plt.scatter(StagnationPoints1X,StagnationPoints1Y,color = 'red')
annotations = ['P2','P1']
for i, label in enumerate(annotations):
    plt.annotate(label, (StagnationPoints1X[i], StagnationPoints1Y[i]),xytext =(StagnationPoints1X[i]-0.05, StagnationPoints1Y[i]-0.2),color = 'blue')

### 2nd flow
flow2 = PotentialFlow(xs,ys)
flow2.AddUniformFlow(V_inf)
flow2.AddDoubletFlow(2*np.pi*V_inf*(r**2))
flow2.AddVortexFlow(Gamma_crit,RadiusZeroPoint = r)
flow2.SolveFlow()
### Calculation of Stagnation Points
flowThetaPoint2 = [-np.pi/2]
flowRadiiPoint2 = [r]

StagnationPoints2X = flowRadiiPoint2*np.cos(flowThetaPoint2)
StagnationPoints2Y = flowRadiiPoint2*np.sin(flowThetaPoint2)

for idx,i in enumerate(flow2.PressureCoefficient):
    for jdx,j in enumerate(i):       
        if (flow2.Xs[idx,jdx]**2 + flow2.Ys[idx,jdx]**2)**0.5 < r:
            flow2.PressureCoefficient[idx,jdx] = 0 
            
flow2.PlotFlow(PotentialFunction=False,PlotInstantly=False,VelocityMagnitude=False,LevelLabels=True, PressureCoefficient=True)
plt.title('C_p Plot around cylinder')
plt.scatter(StagnationPoints2X,StagnationPoints2Y,color = 'red')
annotations = ['P3']
for i, label in enumerate(annotations):
    plt.annotate(label, (StagnationPoints2X[i], StagnationPoints2Y[i]),xytext =(StagnationPoints2X[i]-0.05, StagnationPoints2Y[i]-0.2),color = 'blue')

### 3rd flow
flow3 = PotentialFlow(xs,ys)
flow3.AddUniformFlow(V_inf)
flow3.AddDoubletFlow(2*np.pi*V_inf*(r**2))
flow3.AddVortexFlow(Gamma_3,RadiusZeroPoint = r)
flow3.SolveFlow()
### Calculation of Stagnation Points
flowThetaPoint3 = [-np.pi/2,-np.pi/2]
flowRadiiPoint3 = [(Gamma_3/(4*np.pi*V_inf))+np.sqrt((Gamma_3/(4*np.pi*V_inf))**2 - r**2),(Gamma_3/(4*np.pi*V_inf))-np.sqrt((Gamma_3/(4*np.pi*V_inf))**2 - r**2)]

StagnationPoints3X = flowRadiiPoint3*np.cos(flowThetaPoint3)
StagnationPoints3Y = flowRadiiPoint3*np.sin(flowThetaPoint3)

flow3.PlotFlow(PotentialFunction=False,VelocityMagnitude=False,PlotInstantly=False,LevelLabels=True)
plt.scatter(StagnationPoints3X,StagnationPoints3Y,color = 'red')
annotations = ['P4','P5']
for i, label in enumerate(annotations):
    plt.annotate(label, (StagnationPoints3X[i], StagnationPoints3Y[i]),xytext =(StagnationPoints3X[i]-0.05, StagnationPoints3Y[i]-0.2),color = 'blue')
### Plotting of Pressure coefiecient

def PressureCoefficientPlot(Object):

    Thetas = np.linspace(0,2*np.pi,50)
    Radii = np.multiply(r,np.ones(np.shape(Thetas)))

    Xss = np.multiply(Radii,np.cos(Thetas))
    Yss = np.multiply(Radii,np.sin(Thetas))

    PressC = Object.PressureCoefficientAt(Xss,Yss)

    scale = 0.09

    dXss = np.multiply(scale,np.multiply(np.abs(PressC),np.cos(Thetas)))
    dYss = np.multiply(scale,np.multiply(np.abs(PressC),np.sin(Thetas)))

    fig,ax = plt.subplots()
    ax.set_aspect('equal','box')
    ax.set_xlim(Object.xlimt)
    ax.set_ylim(Object.ylimt)

    ax.contour(Object.Xs,Object.Ys,Object.psi,colors = 'k',linestyles = 'solid', levels = [0])
    ax.grid()

    ax.streamplot(Object.Xs,Object.Ys,Object.u_vel,Object.v_vel,density=2, color='0.8', linewidth= 0.5)
    N_std = 2
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')

    headlenght = 0.04
    for i in range(len(Xss)):
        dx = scale*abs(PressC[i])*np.cos(Thetas[i])
        dy = scale*abs(PressC[i])*np.sin(Thetas[i])
        if PressC[i] <0:
            
            ax.arrow(Xss[i],Yss[i], dx,dy,color = 'blue',length_includes_head = True,head_length = headlenght, head_width = 0.5*headlenght)
        else:
            ax.arrow(Xss[i],Yss[i], dx,dy,color = 'red',length_includes_head = True,head_length = headlenght, head_width = 0.5*headlenght)

    plt.plot(np.add(Xss,dXss),np.add(Yss,dYss),color = '0.3', linewidth = 0.5)

    
def PressureCoefficientPlot2(Objects,Labels):

    Thetas = np.linspace(0,2*np.pi,50)
    Radii = np.multiply(r,np.ones(np.shape(Thetas)))

    Xss = np.multiply(Radii,np.cos(Thetas))
    Yss = np.multiply(Radii,np.sin(Thetas))

    fig,ax = plt.subplots()
    ax.set_xlim([0,2*np.pi])
    for i in range(len(Objects)):
        PressC = Objects[i].PressureCoefficientAt(Xss,Yss)
        ax.plot(Thetas,PressC, label = Labels[i] )

    ax.grid()
    ax.title.set_text('Pressure Coefficient along the arc of the Cylinder')
    plt.legend()
    plt.xlabel('Theta [rad]')
    plt.ylabel('C_p [-]')
    


PressureCoefficientPlot(flow1)
PressureCoefficientPlot(flow2)
PressureCoefficientPlot(flow3)

PressureCoefficientPlot2([flow1,flow2,flow3],['GammaSmaller','GammaCritical','GammaLarger'])

####Pressure Coefficient


#### final plotting
plt.show()