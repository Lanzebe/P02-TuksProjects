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

#Mases of node 1,2,...etc.
Masses = [50,50,500,400,56.8]
SpringElementsK = [65000,65000,32000,32000,75500]
DamperElementsC = [250,250,7500,7500,3840]

SpringElementsConnectivity = [[0,1,0,0],[0,2],[1,3],[2,3]]
DamperElementsConnectivity = [[0,1],[0,2],[1,3],[2,3]]


sys = MultiDOF(Masses,SpringElementsK,DamperElementsC,SpringElementsConnectivity,DamperElementsConnectivity)

sys.SetInitialDisplacements([0,0,0])

sys.SetInitialVelocities([0,0,0])

sys.SetConstraints(DisplacementFunction1,0)


sys.SolveSystem(0,5,res=10000)
sys.plot(0,5)

plt.show()