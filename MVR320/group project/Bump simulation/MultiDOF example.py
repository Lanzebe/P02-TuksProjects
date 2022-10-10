from Vibration import *


### This is a functiond describing the base exitiation.
def DisplacementFunction1(t):
    a = 2
    b = 3
    if (t<b)and(t>a):
        dis = -(t-a)*(t-b)
        vel  = -(t-a) - (t-b)

    else:
        dis = 0
        vel  = 0
    return dis,vel

### This is a function describing a force time function applied to a specified node.
def ForcingFunction1(t):
    omega = 150
    M = 200
    if np.sin(omega*t)>0.8:

        return M
    else: return 0


#Mases of node 1,2,...etc.
Masses = [100,500,56.8]
SpringElementsK = [65000*2,32000*2,754500]
DamperElementsC = [250*2,7500*2,3840]

SpringElementsConnectivity = [[0,1],[1,2],[2,3]]
DamperElementsConnectivity = [[0,1],[1,2],[2,3]]


sys = MultiDOF(Masses,SpringElementsK,DamperElementsC,SpringElementsConnectivity,DamperElementsConnectivity)

k1 = 65000*2
k2 = 32000*2
k3 = 754500

c1 = 250*2
c2 = 7500*2
c3 = 3840

sys.LiveK = [
[k1,-k1,0,0],
[-k1,k1+k2,-k2,0],
[0,-k2,k2+k3,-k3],
[0,0,-k3,k3],
]

sys.C = [
[c1,-c1,0,0],
[-c1,c1+c2,-c2,0],
[0,-c2,c2+c3,-c3],
[0,0,-c3,c3],
]

sys.SetInitialDisplacements([0,0,0,0])

sys.SetInitialVelocities([0,0,0,0])

sys.SetConstraints(DisplacementFunction1,0)


sys.SolveSystem(0,10,res=10000)
sys.plot(0,10)

""" print(sys.M)
print(sys.LiveK)
print(sys.C) """

plt.show()