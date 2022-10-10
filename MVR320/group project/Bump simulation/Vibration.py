#from LSODESolver import LSODECC
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class MultiDOF(object):
    '''
    Zeroth node is the relative fixed node.
    '''
    def __init__(self,Masses,SpringElementsK,DamperElementsC,SpringElementsConnectivity,DamperElementsConnectivity):
        self.Masses = Masses
        self.SpringElementsK = SpringElementsK
        self.DamperElementsC = DamperElementsC
        self.SpringElementsConnectivity = SpringElementsConnectivity
        self.SpringCompOnlyType = []
        self.DamperElementsConnectivity = DamperElementsConnectivity
        self.LiveK = []
        self.C = []
        self.M = []
        self.initX = []
        self.initXp = []
        self.DisplacementFunctions = []
        self.ForcingFunctions = []
        self.y = 0
        self.once = True
        

        ###init calcs
        self.n = len(self.Masses)+1

        self.__BuildK_Matrix()
        self.__BuildC_Matrix()
        self.__BuildM_Matrix()

    def __BuildK_Matrix(self):
        
        self.LiveK = np.zeros([self.n,self.n])
        for idx,i in enumerate(self.SpringElementsK):
            ans = np.zeros([self.n,self.n])
            con = self.SpringElementsConnectivity[idx]
            ans[con[0],con[0]] = i
            ans[con[0],con[1]] = -i
            ans[con[1],con[0]] = -i
            ans[con[1],con[1]] = i
            self.LiveK  = self.LiveK + ans

    def BuildLiveK(self,X):
        self.LiveK = np.zeros([self.n,self.n])
        for idx,i in enumerate(self.SpringElementsK):
            if(self.SpringCompOnlyType[idx]==0)or((X[self.SpringElementsConnectivity[idx][1]]-X[self.SpringElementsConnectivity[idx][0]])<0):
                ans = np.zeros([self.n,self.n])
                con = self.SpringElementsConnectivity[idx]
                ans[con[0],con[0]] = i
                ans[con[0],con[1]] = -i
                ans[con[1],con[0]] = -i
                ans[con[1],con[1]] = i
                self.LiveK = self.LiveK + ans

    def __BuildC_Matrix(self):
        self.C = np.zeros([self.n,self.n])
        for idx,i in enumerate(self.DamperElementsC):
            ans = np.zeros([self.n,self.n])
            con = self.DamperElementsConnectivity[idx]
            ans[con[0],con[0]] = i
            ans[con[0],con[1]] = -i
            ans[con[1],con[0]] = -i
            ans[con[1],con[1]] = i
            self.C = self.C + ans

    def __BuildM_Matrix(self):
        self.M = np.zeros([self.n,self.n])
        self.M[0][0] = 1
        for idx,i in enumerate(self.Masses):
            self.M[idx+1][idx+1] = i

    def DerivativeX(self,X,t):
        Xdis = X[0:self.n]
        Xvel = X[self.n:]
        XdisOld = Xdis
        XvelOld = Xvel
        if self.once:
            print('Xdis,Xvel',Xdis,Xvel)
        
        Xdis,Xvel = self.OverWriteDisplacement(Xdis,Xvel,t)
        if self.once:
            print('Xdis,Xvel',Xdis,Xvel)

            
        #self.BuildLiveK(Xdis)


        xd = Xvel
        
        term1 = np.matmul(-np.linalg.inv(self.M), self.LiveK)		
        term2 = np.matmul(-np.linalg.inv(self.M), self.C)
        term3 = np.matmul(np.linalg.inv(self.M), self.ApplyForcingFuncitons(t))
        
        xdd = np.matmul(term1,Xdis) + np.matmul(term2,Xvel) + term3


        #print('M \n ', self.M)
        #print('K \n', self.LiveK)
        #print('C \n', self.C)
        #print('invM',np.linalg.inv(self.M))
        #print('term1:\n',term1)
        #print('term2:\n',term2)
        #print('applyforcingfunctoins',self.ApplyForcingFuncitons(t))
        #print('term3:\n',term3)
        while self.once:
            print('M \n ', self.M)
            print('K \n', self.LiveK)
            print('C \n', self.C)
            print('invM',np.linalg.inv(self.M))
            print('term1:\n',term1)
            print('term2:\n',term2)
            #print('applyforcingfunctoins',self.ApplyForcingFuncitons(t))
            print('term3:\n',term3)
            #print('term1 OLD:\n',np.matmul(term1,XdisOld) )
            print('term1 next:\n',np.matmul(term1,Xdis) )
            print('term2 next:\n',np.matmul(term2,Xvel))
            print('term3 next:\n',term3)
            self.once = False
        #print('term1*Xdis',np.matmul(term1,Xdis.T))
        #print('xd',xd)
        #print('xdd',xdd)
        return np.append(xd,xdd)

    def SetSpringsToCompressionOnly(self,S):
        self.SpringCompOnlyType = S

    def SetInitialDisplacements(self,X):
        self.initX = np.array(X)



    def SetInitialVelocities(self,Xp):
        self.initXp = np.array(Xp)

    def SetConstraints(self,FuncDescribingX,nodenumber):
        '''
        FuncDescribingX must be a input function that spits out x,v(displacement,velocity) when given t
        '''
        self.DisplacementFunctions = [0]*self.n
        self.DisplacementFunctions[nodenumber] = FuncDescribingX

    def SetForcingFunction(self, FuncDescribingF,NodeNumber):
        '''
        FuncDescribingF must take one input t and output a force for one node NodeNumber
        '''
        self.ForcingFunctions = [0]*self.n
        self.ForcingFunctions[NodeNumber] = FuncDescribingF

    def OverWriteDisplacement(self,X,Xp,t):
        #print(self.DisplacementFunctions)
        for idx,i in enumerate(self.DisplacementFunctions):
            if i!=0:
                
                X[idx],Xp[idx] = i(t)
            #print(X,',',Xp)
        return X,Xp

    def ApplyForcingFuncitons(self,t):
        F = [0]*self.n
        for idx,i in enumerate(self.ForcingFunctions):
            
            if i!=0:
                F[idx]= i(t)
        return F


    def Labels(self):
        l = []
        for i in range(self.n):
            l.append('x' + str(i))
        return l


    def SolveSystem(self,t0,t1,res=1000):
        self.res = res
        t = np.linspace(t0,t1,res)
        InitialConditions = np.append(self.initX,self.initXp)
        #print('initial con:',InitialConditions)
        self.y = odeint(self.DerivativeX,InitialConditions,t)


    def plot(self,t0,t1,FigureNum=1):
        t = np.linspace(t0,t1,self.res)
        plt.figure(FigureNum)
        displacements = self.y[:,0:self.n]
        plt.plot(t,displacements, label = self.Labels())
        plt.xlabel('time')
        plt.ylabel('x(t)')
        plt.legend()
        plt.grid()
        #plt.show()


            


