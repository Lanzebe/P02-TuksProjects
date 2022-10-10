from Vibration import *


### This is a functiond describing the base exitiation.

def DisplacementFunction1(t):
    if (t<5)and(t>4):
        dis = -(t-5)*(t-4)
        vel  = -(t-5) - (t-4)

    else:
        dis = 0
        vel  = 0
    return dis,vel



Masses = [1,1]
SpringElementsK = [4*np.pi**2,4*np.pi**2]
DamperElementsC = [1,2]

SpringElementsConnectivity = [[0,1],[1,2]]
DamperElementsConnectivity = [[0,1],[1,2]]


sys = MultiDOF(Masses,SpringElementsK,DamperElementsC,SpringElementsConnectivity,DamperElementsConnectivity)

sys.SetInitialDisplacements([0,0.25,0])

sys.SetInitialVelocities([0,0,0])

sys.SetConstraints(DisplacementFunction1,0)


sys.SolveSystem(0,10)
sys.plot(0,10)

plt.show()