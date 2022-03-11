from traceback import print_tb
from xml.etree.ElementTree import ElementTree
import numpy as np
    

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
        self.BCxVelocity = 'none'
        self.BCyVelocity = 'none'
        self.BCPressure = 'none'

        self.ConnectionPositionsX = ConnectionPositionsX
        self.ConnectionPositionsY = ConnectionPositionsY

        self.ElementL = None
        self.ElementR = None
        self.ElementT = None
        self.ElementB = None

    def CostFunction(self):
        x = self.CG[0]
        x_l = self.ElementL.CG[0]
        x_r = self.ElementR.CG[0]
        #x_b = self.ElementB.CG[0]
        #x_t = self.ElementT.CG[0]

        y = self.CG[1]
        #y_l = self.ElementL.CG[1]
        #y_r = self.ElementR.CG[1]
        y_b = self.ElementB.CG[1]
        y_t = self.ElementT.CG[1]

        p = self.VelocityX
        p_l = self.ElementL.VelocityX
        p_r = self.ElementR.VelocityX
        p_b = self.ElementB.VelocityX
        p_t = self.ElementT.VelocityX

        u = self.VelocityX
        u_l = self.ElementL.VelocityX
        u_r = self.ElementR.VelocityX
        u_b = self.ElementB.VelocityX
        u_t = self.ElementT.VelocityX

        v = self.VelocityY
        v_l = self.ElementL.VelocityY
        v_r = self.ElementR.VelocityY
        v_t = self.ElementT.VelocityY
        v_b = self.ElementB.VelocityY

        mu = TwoDAirElement.mu
        rho = TwoDAirElement.rho


        return (-mu*(((-u + u_t)/(-y + y_t) - (u - u_b)/(y - y_b))/(-y_b + y_t)**2 + ((-u + u_r)/(-x + x_r) - (u - u_l)/(x - x_l))/(-x_l + x_r)**2) + rho*(u*((u - u_l)/(x - x_l) + (x - x_l)*((-u + u_r)/(-x + x_r) - (u - u_l)/(x - x_l))/(-x_l + x_r)) + v*((u - u_b)/(y - y_b) + (y - y_b)*((-u + u_t)/(-y + y_t) - (u - u_b)/(y - y_b))/(-y_b + y_t))) - (p - p_l)/(x - x_l) - (x - x_l)*((-p + p_r)/(-x + x_r) - (p - p_l)/(x - x_l))/(-x_l + x_r))**2 + (-mu*(((-v + v_t)/(-y + y_t) - (v - v_b)/(y - y_b))/(-y_b + y_t)**2 + ((-v + v_r)/(-x + x_r) - (v - v_l)/(x - x_l))/(-x_l + x_r)**2) + rho*(u*((v - v_l)/(x - x_l) + (x - x_l)*((-v + v_r)/(-x + x_r) - (v - v_l)/(x - x_l))/(-x_l + x_r)) + v*((v - v_b)/(y - y_b) + (y - y_b)*((-v + v_t)/(-y + y_t) - (v - v_b)/(y - y_b))/(-y_b + y_t))) - (p - p_b)/(y - y_b) - (y - y_b)*((-p + p_t)/(-y + y_t) - (p - p_b)/(y - y_b))/(-y_b + y_t))**2 + ((u - u_l)/(x - x_l) + (v - v_b)/(y - y_b) + (x - x_l)*((-u + u_r)/(-x + x_r) - (u - u_l)/(x - x_l))/(-x_l + x_r) + (y - y_b)*((-v + v_t)/(-y + y_t) - (v - v_b)/(y - y_b))/(-y_b + y_t))**2

#Generating nodes
xs = np.linspace(0,10,5)
ys = np.linspace(0,5,4)
Elements = []
for idx, x in enumerate(xs):
    for jdx, y in enumerate(ys):
        if idx == 0 and jdx == 0 :
            type = 'BL_Corner'
            ConnectionPositionsX = ['False',xs[idx+1],'False',x]
            ConnectionPositionsY = ['False',y,'False',ys[jdx+1]]
        elif idx == len(xs)-1 and jdx == 0 :
            type = 'BR_Corner'
            ConnectionPositionsX = [xs[idx-1],'False','False',x]
            ConnectionPositionsY = [y,'False','False',ys[jdx+1]]
        elif idx == 0 and jdx == len(ys)-1 :
            type = 'TL_Corner'
            ConnectionPositionsX = ['False',xs[idx+1],x,'False']
            ConnectionPositionsY = ['False',y,ys[jdx-1],'False']
        elif idx == len(xs)-1 and jdx == len(ys)-1 :
            type = 'TR_Corner'
            ConnectionPositionsX = [xs[idx-1],'False',x,'False']
            ConnectionPositionsY = [y,'False',ys[jdx-1],'False']
        elif idx == 0 and jdx != 0 and jdx != len(ys)-1 :
            type = 'L_side'
            ConnectionPositionsX = ['False',xs[idx+1],x,x]
            ConnectionPositionsY = ['False',y,ys[jdx-1],ys[jdx+1]]
        elif idx == len(xs)-1 and jdx != 0 and jdx != len(ys)-1 :
            type = 'R_side'
            ConnectionPositionsX = [xs[idx-1],'False',x,x]
            ConnectionPositionsY = [y,'False',ys[jdx-1],ys[jdx+1]]
        elif idx != len(xs)-1 and idx != 0 and jdx == 0 :
            type = 'B_side'
            ConnectionPositionsX = [xs[idx-1],xs[idx+1],'False',x]
            ConnectionPositionsY = [y,y,'False',ys[jdx+1]]
        elif idx != len(xs)-1 and idx != 0 and jdx == len(ys)-1 :
            type = 'T_side'
            ConnectionPositionsX = [xs[idx-1],xs[idx+1],x,'False']
            ConnectionPositionsY = [y,y,ys[jdx-1],'False']
        else:
            type = 'Normal'
            ConnectionPositionsX = [xs[idx-1],xs[idx+1],x,x]
            ConnectionPositionsY = [y,y,ys[jdx-1],ys[jdx+1]]
        
        Elements.append(TwoDAirElement(x,y,ConnectionPositionsX,ConnectionPositionsY,type))
###Meshing
for i in Elements:    
    CoordsL = [i.ConnectionPositionsX[0],i.ConnectionPositionsY[0]]
    CoordsR = [i.ConnectionPositionsX[1],i.ConnectionPositionsY[1]]
    CoordsB = [i.ConnectionPositionsX[2],i.ConnectionPositionsY[2]]
    CoordsT = [i.ConnectionPositionsX[3],i.ConnectionPositionsY[3]]
    for j in Elements:
        Coords = j.CG
        if Coords == CoordsL and CoordsL[0]!='False':
            i.ElementL = j
        elif Coords == CoordsR and CoordsR[0]!='False':
            i.ElementR = j
        elif Coords == CoordsB and CoordsB[0]!='False':
            i.ElementB = j
        elif Coords == CoordsT and CoordsT[0]!='False':
            i.ElementT = j

for i in Elements:
    print(i.CostFunction())