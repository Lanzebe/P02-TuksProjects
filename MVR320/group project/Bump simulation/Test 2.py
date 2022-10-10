from Vibration import *


### This is a functiond describing the base exitiation.

def DisplacementFunction1(t):
    if (t<3)and(t>2):
        dis = -(t-2)*(t-3)
        vel  = -(t-2) - (t-3)

    else:
        dis = 0
        vel  = 0
    return dis,vel



Masses = [10]
SpringElementsK = [200000]
DamperElementsC = [10]

SpringElementsConnectivity = [[0,1]]
DamperElementsConnectivity = [[0,1]]


sys = MultiDOF(Masses,SpringElementsK,DamperElementsC,SpringElementsConnectivity,DamperElementsConnectivity)

sys.SetInitialDisplacements([0,0])

sys.SetInitialVelocities([0,0])

sys.SetConstraints(DisplacementFunction1,0)


sys.SolveSystem(0,5,res=10000)
sys.plot(0,5)

plt.show()