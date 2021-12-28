from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

class DD_KineBody(object):

    def __init__(self,CenterOfGravityPos = [0,0], Mass = 1, RotationalInertiaZZ = 1 ):
        self.CGLocal = CenterOfGravityPos
        self.Mass = Mass
        self.RotationalInertiaZZ = RotationalInertiaZZ
        self.PointsLocal = []
        self.VelLocal = []



        ## Global Parameters
        self.Pos = np.array([0.1,0.1])
        self.PosD = np.array([0.1,0.1])
        self.PosDD = np.array([0.1,0.1])
        self.Theta = 0.1
        self.ThetaD = 0.1
        self.TheraDD = 0.1

        self.PointsGlobal = np.array([[]])
        self.VelGlobal = np.array([[]])


    def AddLocalPointsToBody(self,ListOfPoints):
        for i in ListOfPoints:
            self.PointsLocal.append(i)
            self.VelLocal.append([0,0])


    def UpdateGlobalPointPositions(self):
        A = []
        for i in self.PointsLocal:
            x = i[0]
            y = i[1]
            #l = (x**2 + y**2)**0.5
            A.append([self.Pos[0] + x*np.cos(self.Theta) - y*np.sin(self.Theta),self.Pos[1] + x*np.sin(self.Theta) + y*np.cos(self.Theta)])

        self.PointsGlobal = np.array(A)

    def UpdateGlobalPointVelocities(self):
        A = []
        for i in self.PointsLocal:
            x = i[0]
            y = i[1]
            #l = (x**2 + y**2)**0.5
            A.append([self.PosD[0] - self.ThetaD*x*np.sin(self.Theta) - self.ThetaD*y*np.cos(self.Theta), self.PosD[1] + self.ThetaD*x*np.cos(self.Theta) - self.ThetaD*y*np.sin(self.Theta)])

        self.VelGlobal = np.array(A)


    # Lets handle duplicates later ok
    """     for i in self.PointsLocal:
            if any(i in self.PointsLocal[i,:]): """



class KinematicSystem(object):
    def __init__(self):
        self.ComponentBodies = [] 
        self.FixedPoints = []
        self.Connections = []
        self.BCBodyAngle = []
        self.BCBodyAngleVelocity = []


    def AddComponnentsToSystem(self, ComponentBodies):
        for i in ComponentBodies:
            self.ComponentBodies.append(i)

    def SetFixedPoints(self,Body,point,PointPosition = [0,0]):
        self.FixedPoints.append([Body,point,PointPosition])

    def SetConnectComponents(self,Comp1,Comp2,point1,point2):
        self.Connections.append([Comp1,Comp2,point1,point2])

    def SetAngleComponent(self, Body, angle):
        self.BCBodyAngle.append([Body,angle])

    def SetAngularVelocityComponent(self, Body, angle):
        self.BCBodyAngleVelocity.append([Body,angle])

    def PositionalCalc(self,init_guess):
        #print(init_guess)


        for i in range(len(self.ComponentBodies)):  
            self.ComponentBodies[i].Theta = init_guess[3*i]
            self.ComponentBodies[i].Pos[0] = init_guess[3*i+1]
            self.ComponentBodies[i].Pos[1] = init_guess[3*i+2]
            self.ComponentBodies[i].UpdateGlobalPointPositions()


        xs = []
        ys = []

        for i in self.FixedPoints:
            xs.append(abs(i[0].PointsGlobal[i[1]][0] - i[2][0]))
            ys.append(abs(i[0].PointsGlobal[i[1]][1] - i[2][1]))

        conx = []
        cony = []

        for i in self.Connections:

            conx.append(abs(i[0].PointsGlobal[i[2]][0]-i[1].PointsGlobal[i[3]][0]))
            cony.append(abs(i[0].PointsGlobal[i[2]][1]-i[1].PointsGlobal[i[3]][1]))

        angs = []
        for i in self.BCBodyAngle:
            angs.append(abs(i[0].Theta- i[1]))
            

        ans = xs+ys+conx+cony+angs
        #print(ans)

        return np.array(ans)

    def SolvePositions(self):
    
        InitialGuess  = []

        for i in self.ComponentBodies:
            InitialGuess.append(i.Theta)
            InitialGuess.append(i.Pos[0])
            InitialGuess.append(i.Pos[1])


        sol = fsolve(self.PositionalCalc, InitialGuess)
        #print(sol)
        
        
        
        for i in range(len(self.ComponentBodies)):
            self.ComponentBodies[i].Theta = sol[3*i]
            self.ComponentBodies[i].Pos[0] = sol[3*i+1]
            self.ComponentBodies[i].Pos[1] = sol[3*i+2]
            self.ComponentBodies[i].UpdateGlobalPointPositions()

    def SpeedCalc(self,init_guess):
        #print(init_guess)


        for i in range(len(self.ComponentBodies)):  
            self.ComponentBodies[i].ThetaD = init_guess[3*i]
            self.ComponentBodies[i].PosD[0] = init_guess[3*i+1]
            self.ComponentBodies[i].PosD[1] = init_guess[3*i+2]
            self.ComponentBodies[i].UpdateGlobalPointVelocities()


        xs = []
        ys = []

        for i in self.FixedPoints:
            xs.append(abs(i[0].VelGlobal[i[1]][0] - 0))
            ys.append(abs(i[0].VelGlobal[i[1]][1] - 0))

        conx = []
        cony = []

        for i in self.Connections:

            conx.append(abs(i[0].VelGlobal[i[2]][0]-i[1].VelGlobal[i[3]][0]))
            cony.append(abs(i[0].VelGlobal[i[2]][1]-i[1].VelGlobal[i[3]][1]))

        angs = []
        for i in self.BCBodyAngleVelocity:
            angs.append(abs(i[0].ThetaD - i[1]))
            

        ans = xs+ys+conx+cony+angs
        #print(ans)

        return np.array(ans)

    def SolveVelocities(self):
    
        InitialGuess  = []

        for i in self.ComponentBodies:
            InitialGuess.append(i.ThetaD)
            InitialGuess.append(i.PosD[0])
            InitialGuess.append(i.PosD[1])

        sol = fsolve(self.SpeedCalc, InitialGuess)
        #print(sol)
        
        
        
        for i in range(len(self.ComponentBodies)):
            self.ComponentBodies[i].ThetaD = sol[3*i]
            self.ComponentBodies[i].PosD[0] = sol[3*i+1]
            self.ComponentBodies[i].PosD[1] = sol[3*i+2]
            self.ComponentBodies[i].UpdateGlobalPointVelocities()

    def BodyPlotStack(self):
        for i in self.ComponentBodies:
            plt.plot(i.PointsGlobal[:,0],i.PointsGlobal[:,1])



        







""" Body1 = DD_KineBody(CenterOfGravityPos=[1,2],Mass=10,RotationalInertiaZZ=10)
Body1.AddLocalPointsToBody([[0,0],[1,0]])
Body2 = DD_KineBody(CenterOfGravityPos=[1,2],Mass=20,RotationalInertiaZZ=22)
Body2.AddLocalPointsToBody([[0,0],[2,0]])






KSystem = KinematicSystem()



KSystem.AddComponnentsToSystem([Body1,Body2])
KSystem.SetConnectComponents(Body1,Body2,1,0)
KSystem.SetFixedPoints(Body1,0)
KSystem.SetFixedPoints(Body2,1,PointPosition=[2,2])



KSystem.SolvePositions()
plt.figure()

KSystem.BodyPlot()
plt.axis("equal")
plt.grid()
plt.show()


print("Body points:")
print(Body1.PointsGlobal)
print(Body2.PointsGlobal) """

""" print(KSystem.ComponentBodies)
print(KSystem.Connections)
print(KSystem.FixedPoints) """


