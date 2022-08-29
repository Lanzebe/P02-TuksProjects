import matplotlib.pyplot as plt
import numpy as np

class PotentialFlow:
    def __init__(self,xs,ys, AmbientPressure = 1.01325E5 , AmbientDensity = 1.2250) -> None:
        self.xlimt = (np.min(xs),np.max(xs))
        self.ylimt = (np.min(ys),np.max(ys))
        self.Xs, self.Ys = np.meshgrid(xs,ys)
        self.ElementaryFlows = []
        self.V_inf = 0.001
        self.psi = np.empty(np.shape(self.Xs))
        self.phi = np.empty(np.shape(self.Xs))
        self.u_vel = np.empty(np.shape(self.Xs))
        self.v_vel = np.empty(np.shape(self.Xs))
        self.Vel = np.empty(np.shape(self.Xs))
        self.PressureCoefficient = np.empty(np.shape(self.Xs))
        self.AmbientPressure = AmbientPressure
        self.AmbientDensity = AmbientDensity

    def SolveStreamLineFunction(self):
        for i in self.ElementaryFlows:
            self.psi = np.add(self.psi,i.CalculateStreamLine())

    def SolveVelocityPotentialFunction(self):
        for i in self.ElementaryFlows:
            self.phi = np.add(self.phi, i.CalculateVelocityPotential())

    def SolveVelocities(self):
        for i in self.ElementaryFlows:
            self.u_vel = np.add(self.u_vel, i.CalculateU_Velocity())
            self.v_vel = np.add(self.v_vel, i.CalculateV_Velocity())

        self.Vel = np.sqrt(np.square(self.u_vel) + np.square(self.v_vel))

    def SolvePressureCoefficient(self):
        self.PressureCoefficient = np.subtract(1 , np.square(np.divide(self.Vel,self.V_inf)))

    def PressureCoefficientAt(self, X , Y):
        store = [self.Xs,self.Ys]
        self.Xs,self.Ys = np.array(X),np.array(Y)

        u_vel = np.zeros(np.shape(X))
        v_vel = np.zeros(np.shape(X))

        for i in self.ElementaryFlows:
            #ubuf = i.CalculateU_Velocity()
            #vbuf = i.CalculateV_Velocity()
            #print(i,ubuf)
            #print(i,vbuf)
            u_vel = np.add(u_vel, i.CalculateU_Velocity())
            v_vel = np.add(v_vel, i.CalculateV_Velocity())

            #print("looking",u_vel)
            #print("looking",v_vel)
        
        Vel = np.sqrt(np.add(np.square(u_vel) , np.square(v_vel)))

        self.Xs,self.Ys = store[0],store[1]

        return np.subtract(1 , np.square(np.divide(Vel,self.V_inf)))

    def AddUniformFlow(self,V_inf,Theta = 0):
        self.V_inf += V_inf
        self.ElementaryFlows.append(UniformFlow(V_inf,FlowField=self,Theta=Theta))

    def AddSourceFlow(self,Lambda,Xposition = 0,Yposition = 0):
        self.ElementaryFlows.append(SourceFlow(Lambda,Xposition,Yposition,FlowField=self))

    def AddDoubletFlow(self,Kappa,Xposition = 0,Yposition = 0):
        self.ElementaryFlows.append(DoubletFlow(Kappa,Xposition,Yposition,FlowField=self))

    def AddVortexFlow(self,Gamma,Xposition = 0,Yposition = 0, RadiusZeroPoint = 1):
        self.ElementaryFlows.append(VortexFlow(Gamma,Xposition,Yposition,FlowField=self,RadiusZeroPoint = RadiusZeroPoint))

    def SolveFlow(self):
        self.SolveStreamLineFunction()
        self.SolveVelocityPotentialFunction()
        self.SolveVelocities()
        self.SolvePressureCoefficient()

        ShapeSize = np.shape(self.Xs)

        FlowExtreames = [self.psi[:,0],self.psi[:,ShapeSize[1]-1],self.psi[0,:],self.psi[ShapeSize[0]-1,:]]
        PoteExtreames = [self.phi[:,0],self.phi[:,ShapeSize[1]-1],self.phi[0,:],self.phi[ShapeSize[0]-1,:]]

        self.FlowLines = np.linspace(np.min(FlowExtreames),np.max(FlowExtreames),16+1)
        self.PotentialLines = np.linspace(np.min(PoteExtreames),np.max(PoteExtreames),16+1)

    def PlotFlow(self, StreamFunction = False, StreamLines = True, PotentialFunction = False, PressureCoefficient = False, VelocityMagnitude = False, LevelLabels = False, VelocityLevels = [0], Grid = False, PlotInstantly = True):
        fig,ax = plt.subplots()
        ax.set_aspect('equal','box')
        ax.set_xlim(self.xlimt)
        ax.set_ylim(self.ylimt)
        ax.contour(self.Xs,self.Ys,self.psi,colors = 'k',linestyles = 'solid', levels = VelocityLevels)

        if Grid:
            ax.grid()
        
        if VelocityMagnitude:
            colorbar = ax.contourf(self.Xs,self.Ys,self.Vel)

        if PressureCoefficient:
            colorbar = ax.contourf(self.Xs,self.Ys,self.PressureCoefficient)
            
        if StreamLines:
            ax.streamplot(self.Xs,self.Ys,self.u_vel,self.v_vel,density=2, color='0.2', linewidth= 0.5)
            

        if StreamFunction:
            mean = np.mean(self.psi)
            N_std = 2
            levels = np.linspace(np.subtract(mean,N_std*np.std(self.psi)),np.add(mean,N_std*np.std(self.psi)),27)
            flow = ax.contour(self.Xs,self.Ys,self.psi,colors = '0.2',linestyles = 'solid', levels = levels)
            if LevelLabels:
                ax.clabel(flow,flow.levels)

        if PotentialFunction:
            pote = ax.contour(self.Xs,self.Ys,self.phi,colors = 'r',linestyles = 'solid', levels = self.PotentialLines)
            if LevelLabels:
                ax.clabel(pote,pote.levels)


        plt.xlabel('x [m]')
        plt.ylabel('y [m]')

        if VelocityMagnitude or PressureCoefficient:
            plt.colorbar(colorbar)
        
        if PlotInstantly:
            plt.show()

class UniformFlow:
    def __init__(self,V_inf, FlowField,Theta = 0) -> None:
        self.V_inf = V_inf
        self.Theta = Theta
        self.FlowField = FlowField
    def CalculateStreamLine(self):
        return np.multiply(self.FlowField.Ys,self.V_inf)
    def CalculateVelocityPotential(self):
        return np.multiply(self.FlowField.Xs,self.V_inf)
    def CalculateU_Velocity(self):
        return np.multiply(np.divide(self.V_inf,self.FlowField.AmbientDensity),np.ones(np.shape(self.FlowField.Xs)))
    def CalculateV_Velocity(self):
        return np.multiply(np.multiply(0,self.FlowField.Xs),np.ones(np.shape(self.FlowField.Xs)))

class SourceFlow:
    def __init__(self,Lambda,Xposition,Yposition, FlowField) -> None:
        self.Lambda = Lambda
        self.Position = [Xposition,Yposition]
        self.FlowField = FlowField
        self.radii = np.sqrt(np.add(np.square(np.subtract(self.FlowField.Xs,Xposition)) , np.square(np.subtract(self.FlowField.Ys,Yposition))))
        self.Theta = np.arctan2(np.subtract(self.FlowField.Ys,Yposition),np.subtract(self.FlowField.Xs,Xposition))
        for idx,i in enumerate(self.Theta):
            for jdx,j in enumerate(i):       
                if j < 0:
                    self.Theta[idx,jdx] = j + 2 * np.pi
        #print(np.min(self.Theta))
        #print(np.max(self.Theta))
    def CalculateStreamLine(self):
        return (self.Lambda/(2*np.pi))*self.Theta
    def CalculateVelocityPotential(self):
        return (self.Lambda/(2*np.pi))*np.log(self.radii)
    def CalculateU_Velocity(self):
        return np.multiply((self.Lambda/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(np.subtract(self.FlowField.Xs,self.Position[0]),np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys)))))
    def CalculateV_Velocity(self):
        return np.multiply((-self.Lambda/(2*np.pi*self.FlowField.AmbientDensity)),np.divide((np.subtract(self.Position[1] , self.FlowField.Ys)),np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys)))))

class DoubletFlow:
    def __init__(self,Kappa,Xposition,Yposition, FlowField) -> None:
        self.Kappa = Kappa
        self.Position = [Xposition,Yposition]
        self.FlowField = FlowField
        self.radii = np.sqrt(np.add(np.square(np.subtract(self.FlowField.Xs,Xposition)) , np.square(np.subtract(self.FlowField.Ys,Yposition))))
        self.Theta = np.arctan2(np.subtract(self.FlowField.Ys,Yposition),np.subtract(self.FlowField.Xs,Xposition))
        for idx,i in enumerate(self.Theta):
            for jdx,j in enumerate(i):       
                if j < 0:
                    self.Theta[idx,jdx] = j + 2 * np.pi

        #print(np.min(self.Theta))
        #print(np.max(self.Theta))
    def CalculateStreamLine(self):
        return -(self.Kappa/(2*np.pi))*(np.divide(np.sin(self.Theta),self.radii))
    def CalculateVelocityPotential(self):
        return (self.Kappa/(2*np.pi))*(np.divide(np.cos(self.Theta),self.radii))
    def CalculateU_Velocity(self):
        #print(np.multiply((-self.Kappa/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(np.subtract(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys))),np.square(np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys)))))))
        return np.multiply((-self.Kappa/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(np.subtract(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys))),np.square(np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys))))))
    def CalculateV_Velocity(self):
        return np.multiply((-self.Kappa/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(2*np.multiply(np.subtract(self.Position[0] , self.FlowField.Xs),np.subtract(self.Position[1] , self.FlowField.Ys)),np.square(np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys))))))

class VortexFlow:
    def __init__(self,Gamma,Xposition,Yposition, FlowField, RadiusZeroPoint = 1) -> None:
        self.Gamma = Gamma
        self.Position = [Xposition,Yposition]
        self.FlowField = FlowField
        self.RadiusZeroPoint = RadiusZeroPoint
        self.radii = np.sqrt(np.add(np.square(np.subtract(self.FlowField.Xs,Xposition)) , np.square(np.subtract(self.FlowField.Ys,Yposition))))
        self.Theta = np.arctan2(np.subtract(self.FlowField.Ys,Yposition),np.subtract(self.FlowField.Xs,Xposition))
        for idx,i in enumerate(self.Theta):
            for jdx,j in enumerate(i):       
                if j < 0:
                    self.Theta[idx,jdx] = j + 2 * np.pi
    def CalculateStreamLine(self):
        return (self.Gamma/(2*np.pi))*np.log(np.divide(self.radii,self.RadiusZeroPoint))
    def CalculateVelocityPotential(self):
        return -(self.Gamma/(2*np.pi))*self.Theta

    def CalculateU_Velocity(self):
        #print(np.multiply((self.Gamma/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(np.subtract(self.FlowField.Ys,self.Position[1]),np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys))))))
        return np.multiply((self.Gamma/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(np.subtract(self.FlowField.Ys,self.Position[1]),np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys)))))
    def CalculateV_Velocity(self):
        return np.multiply((-self.Gamma/(2*np.pi*self.FlowField.AmbientDensity)),np.divide(np.subtract(self.FlowField.Xs,self.Position[0]),np.add(np.square(np.subtract(self.Position[0] , self.FlowField.Xs)) , np.square(np.subtract(self.Position[1] , self.FlowField.Ys)))))


""" n = 1000
xs = np.linspace(-20,20,n)
ys = np.linspace(-10,10,n)

#print(ys)
#print(xs)

flow = PotentialFlow(xs,ys)
flow.AddUniformFlow(10)
#flow.AddSourceFlow(-100,5,0)
flow.AddSourceFlow(100,-15,-5)
flow.AddSourceFlow(-100,-10,-5)
flow.AddDoubletFlow(1000,0,0)
#flow.AddSourceFlow(100)
#flow.AddDoubletFlow(1000,5,0)
flow.AddVortexFlow(100,-10,5)
flow.AddVortexFlow(-100,-5,5)
flow.SolveFlow()
flow.PlotFlow(PotentialFunction=False,VelocityMagnitude=False) """
