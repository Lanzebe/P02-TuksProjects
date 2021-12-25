import sys
sys.path.append(r'E:\OneDrive\ENGINEERING\01-ENGINEERING TOOLS\engineering_library')
from FourierApprox import *
from Vibration import *
import numpy as np
import matplotlib.pyplot as plt
import cmath


t0 = 0
t1 = 0.25

def approx_fun(t):    
    if t<0.1:
        return (20/0.1)*t
    else:
        return 20
 
series = FourierSeries(t0,t1,N_max_series=40,Integration_n=1000)
series.calc_cs(approx_fun)

t = np.linspace(t0,t1*2,1000)

y = []

for i in t:
    y.append(approx_fun(i))


def gen_XY(T,obj):
    X = []
    Y = []
    for i in T:
        ans = obj.output(i)
        X.append(ans.real)
        Y.append(ans.imag)
    return X,Y

X2,Y2 = gen_XY(t,series)

plt.figure(1)
#plt.plot(t,X1)
plt.plot(t,X2)
plt.plot(t,y)
plt.grid()
#plt.show()


for idx,i in enumerate(series.C):
    text = 'Constant:'+ str(idx+1)
    #print(text,i.real)

def DisplacementFunction1(t):
    return 0, 0

def ForcingFunction1(t):
    return series.output(t).real


Masses = [2]
SpringElementsK = [1500]
DamperElementsC = [10]
SpringElementsConnectivity = [[0,1]]
DamperElementsConnectivity = [[0,1]]

sys = MultiDOF(Masses,SpringElementsK,DamperElementsC,SpringElementsConnectivity,DamperElementsConnectivity)
sys.SetInitialDisplacements([0,0])
sys.SetInitialVelocities([0,0])
sys.SetConstraints(DisplacementFunction1,0)
sys.SetForcingFunction(ForcingFunction1,1)
print('masses:',sys.Masses)
print('n:',sys.n)
sys.SolveSystem(t0,t1*2)
sys.plot(t0,t1*2,2)
plt.show()