import tkinter as tk
import numpy as np
from scipy import optimize
    

class Mesh:
    def __init__(self):
        self.Elements = []
        self.Residual = 0
        
    def InitializeElements(self):
        for i in self.Elements:
            i.InitilizeElementValues()

    def CalculateResidualReturnFunction(self):
        self.InitializeElements()

        Residual = []
        for i in self.Elements:
            Residual.append(i.CostFunction())

        self.Residual = np.sum(Residual)


    def GenerateInputVector(self):
        A = []
        for i in self.Elements:
            if not i.BCxVelocity:
                A.append(0)
            if not i.BCyVelocity:
                A.append(0)
            if not i.BCPressure:
                A.append(0)
        return A

    def Solve(self):
        A = self.GenerateInputVector()

        
        sol = optimize.minimize(self.MinimizeFunction,A)
        print(sol)
        

    def MinimizeFunction(self,A):
        count = 0

        res = []
        for i in self.Elements:
            i.InitilizeElementValues()
            if not i.BCxVelocity:
                i.VelocityX = A[count]
                count += 1
            if not i.BCyVelocity:
                i.VelocityY = A[count]
                count += 1
            if not i.BCPressure:
                i.Pressure = A[count]
                count += 1

        self.CalculateResidualReturnFunction()
        #print(res)
        return self.Residual


class TwoDAirElement:
    CountElements = 1
    mu = 2*(10**(-5))
    rho = 1.2 
    
    def __init__(self,x,y,ConnectionPositionsX,ConnectionPositionsY,ElementType = 'Normal'):
        """ 
        Element Type is either Normal, BL_Corner, TL_Corner, BR_Corner, TR_Corner, L_Side, R_Side, T_Side, or B_Side.
        """
        super(TwoDAirElement,self)
        self.ElementID = self.CountElements
        self.ElementType = ElementType
        TwoDAirElement.CountElements += 1

        self.CG = [x,y]
        self.VelocityX = 0
        self.VelocityY = 0
        self.Presssure = 0
        self.BCxVelocity = None
        self.BCyVelocity = None
        self.BCPressure = None

        self.ConnectionPositionsX = ConnectionPositionsX
        self.ConnectionPositionsY = ConnectionPositionsY

        self.ElementL = None
        self.ElementR = None
        self.ElementT = None
        self.ElementB = None

        self.Residual = 0

    def CostFunction(self):
        x = self.CG[0]
        y = self.CG[1]

        u = self.VelocityX
        v = self.VelocityY
        p = self.Presssure



        ### conditions if elements around current element does not exist:
        if not self.ElementL:
            x = (x + self.ElementR.CG[0])/2
            x_l = self.CG[0]
            y_l = self.CG[1]
            p_l = self.Presssure
            u_l = self.VelocityX
            v_l = self.VelocityY
        else:
            x_l = self.ElementL.CG[0]
            y_l = self.ElementL.CG[1]
            p_l = self.ElementL.Presssure
            u_l = self.ElementL.VelocityX
            v_l = self.ElementL.VelocityY

        if not self.ElementR:
            x = (x + self.ElementL.CG[0])/2
            x_r = self.CG[0]
            y_r = self.CG[1]
            p_r = self.Presssure
            u_r = self.VelocityX
            v_r = self.VelocityY
        else:
            x_r = self.ElementR.CG[0]
            y_r = self.ElementR.CG[1]
            p_r = self.ElementR.Presssure
            u_r = self.ElementR.VelocityX
            v_r = self.ElementR.VelocityY

        if not self.ElementB:
            y = (y+ self.ElementT.CG[1])/2
            x_b = self.CG[0]
            y_b = self.CG[1]
            p_b = self.Presssure
            u_b = self.VelocityX
            v_b = self.VelocityY
        else:
            x_b = self.ElementB.CG[0]
            y_b = self.ElementB.CG[1]
            p_b = self.ElementB.Presssure
            u_b = self.ElementB.VelocityX
            v_b = self.ElementB.VelocityY

        if not self.ElementT:
            y = (y + self.ElementB.CG[1])/2
            x_t = self.CG[0]
            y_t = self.CG[1]
            p_t = self.Presssure
            u_t = self.VelocityX
            v_t = self.VelocityY
        else:
            x_t = self.ElementT.CG[0]
            y_t = self.ElementT.CG[1]
            p_t = self.ElementT.Presssure
            u_t = self.ElementT.VelocityX
            v_t = self.ElementT.VelocityY

        
        dudx = ((u-u_l)/(x-x_l))+((((u_r-u)/(x_r-x))-((u-u_l)/(x-x_l)))*((x-x_l)/(x_r-x_l)))
        dudy = ((u-u_b)/(y-y_b))+((((u_t-u)/(y_t-y))-((u-u_b)/(y-y_b)))*((y-y_b)/(y_t-y_b)))
        dvdx = ((v-v_l)/(x-x_l))+((((v_r-v)/(x_r-x))-((v-v_l)/(x-x_l)))*((x-x_l)/(x_r-x_l)))
        dvdy = ((v-v_b)/(y-y_b))+((((v_t-v)/(y_t-y))-((v-v_b)/(y-y_b)))*((y-y_b)/(y_t-y_b)))
        dpdx = ((p-p_l)/(x-x_l))+((((p_r-p)/(x_r-x))-((p-p_l)/(x-x_l)))*((x-x_l)/(x_r-x_l)))
        dpdy = ((p-p_b)/(y-y_b))+((((p_t-p)/(y_t-y))-((p-p_b)/(y-y_b)))*((y-y_b)/(y_t-y_b)))
        dduddxx = (((u_r-u)/(x_r-x))-((u-u_l)/(x-x_l)))/((x_r-x_l)**2)
        dduddyy = (((u_t-u)/(y_t-y))-((u-u_b)/(y-y_b)))/((y_t-y_b)**2)
        ddvddxx = (((v_r-v)/(x_r-x))-((v-v_l)/(x-x_l)))/((x_r-x_l)**2)
        ddvddyy = (((v_t-v)/(y_t-y))-((v-v_b)/(y-y_b)))/((y_t-y_b)**2)



            
        #print(((p-p_l)/(x-x_l)))


        mu = TwoDAirElement.mu
        rho = TwoDAirElement.rho


        Continuity = dudx + dvdy
        MommentumX = rho*(u*dudx + v*dudy) - dpdx - mu*(dduddxx + dduddyy)
        MommentumY = rho*(u*dvdx + v*dvdy) - dpdy - mu*(ddvddxx + ddvddyy)
        COST = Continuity**2  +  MommentumX**2 +  MommentumY**2 

        #print(COST)
        if COST == None:
            print("This is element cost:",COST)
            print(Continuity)
            print(MommentumX)
            print(MommentumY)
            print('x',x);            print('x_l',x_l);            print('x_r',x_r);            print('x_b',x_b);            print('x_t',x_t)
            print('y',y);            print('y_l',y_l);            print('y_r',y_r);            print('y_b',y_b);            print('y_t',y_t)
            print('u',u);            print('u_l',u_l);            print('u_r',u_r);            print('u_b',u_b);            print('u_t',u_t)
            print('v',v);            print('v_l',v_l);            print('v_r',v_r);            print('v_b',v_b);            print('v_t',v_t)
            print('p',p);            print('p_l',p_l);            print('p_r',p_r);            print('p_b',p_b);            print('p_t',p_t)
            quit()

        self.Residual = COST
        return COST

    def InitilizeElementValues(self):
        if self.BCxVelocity != None:
            self.VelocityX = self.BCxVelocity


        if self.BCyVelocity != None:
            self.VelocityY = self.BCyVelocity


        if self.BCPressure != None:
            self.Presssure = self.BCPressure





#Generating nodes
xs = np.linspace(0,10,11)
ys = np.linspace(0,5,5)
Mesh = Mesh()


for idx, x in enumerate(xs):
    for jdx, y in enumerate(ys):
        if idx == 0 and jdx == 0 :
            type = 'BL_Corner'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = ['False',xs[idx+1],'False',x]
            ConnectionPositionsY = ['False',y,'False',ys[jdx+1]]
        elif idx == len(xs)-1 and jdx == 0 :
            type = 'BR_Corner'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = [xs[idx-1],'False','False',x]
            ConnectionPositionsY = [y,'False','False',ys[jdx+1]]
        elif idx == 0 and jdx == len(ys)-1 :
            type = 'TL_Corner'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = ['False',xs[idx+1],x,'False']
            ConnectionPositionsY = ['False',y,ys[jdx-1],'False']
        elif idx == len(xs)-1 and jdx == len(ys)-1 :
            type = 'TR_Corner'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = [xs[idx-1],'False',x,'False']
            ConnectionPositionsY = [y,'False',ys[jdx-1],'False']
        elif idx == 0 and jdx != 0 and jdx != len(ys)-1 :
            type = 'L_side'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = ['False',xs[idx+1],x,x]
            ConnectionPositionsY = ['False',y,ys[jdx-1],ys[jdx+1]]
        elif idx == len(xs)-1 and jdx != 0 and jdx != len(ys)-1 :
            type = 'R_side'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = [xs[idx-1],'False',x,x]
            ConnectionPositionsY = [y,'False',ys[jdx-1],ys[jdx+1]]
        elif idx != len(xs)-1 and idx != 0 and jdx == 0 :
            type = 'B_side'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = [xs[idx-1],xs[idx+1],'False',x]
            ConnectionPositionsY = [y,y,'False',ys[jdx+1]]
        elif idx != len(xs)-1 and idx != 0 and jdx == len(ys)-1 :
            type = 'T_side'
            XVel = 0
            YVel = 0
            Pres = None
            ConnectionPositionsX = [xs[idx-1],xs[idx+1],x,'False']
            ConnectionPositionsY = [y,y,ys[jdx-1],'False']
        else:
            type = 'Normal'
            XVel = None
            YVel = None
            Pres = None
            ConnectionPositionsX = [xs[idx-1],xs[idx+1],x,x]
            ConnectionPositionsY = [y,y,ys[jdx-1],ys[jdx+1]]
        
        Mesh.Elements.append(TwoDAirElement(x,y,ConnectionPositionsX,ConnectionPositionsY,type))
        Mesh.Elements[-1].BCxVelocity = XVel
        Mesh.Elements[-1].BCyVelocity = YVel
        Mesh.Elements[-1].BCPressure = Pres

###Meshing
for i in Mesh.Elements:    
    CoordsL = [i.ConnectionPositionsX[0],i.ConnectionPositionsY[0]]
    CoordsR = [i.ConnectionPositionsX[1],i.ConnectionPositionsY[1]]
    CoordsB = [i.ConnectionPositionsX[2],i.ConnectionPositionsY[2]]
    CoordsT = [i.ConnectionPositionsX[3],i.ConnectionPositionsY[3]]
    for j in Mesh.Elements:
        Coords = j.CG
        if Coords == CoordsL and CoordsL[0]!='False':
            i.ElementL = j
        elif Coords == CoordsR and CoordsR[0]!='False':
            i.ElementR = j
        elif Coords == CoordsB and CoordsB[0]!='False':
            i.ElementB = j
        elif Coords == CoordsT and CoordsT[0]!='False':
            i.ElementT = j



Mesh.Elements[2].BCxVelocity = 1
Mesh.Elements[2].BCyVelocity = 0
Mesh.Elements[2].BCPressure = 100E3

Mesh.Elements[-2].BCxVelocity = 1
Mesh.Elements[-2].BCyVelocity = 0
Mesh.Elements[-2].BCPressure = 100E3


#Mesh.CalculateResidualReturnFunction()
Mesh.Solve()