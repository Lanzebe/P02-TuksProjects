from Vibration import *

### This is a functiond describing the base exitiation.
raw_data = np.load('X-Y random road.npy')
""" plt.figure()
plt.plot(raw_data[0],raw_data[1])
plt.show() """

speed = 40 #km.h
v = speed*1000/3600
t2 = 1000/v
def BaseExitation1(t):
    
    findx = t*v
    x = raw_data[0]
    y = raw_data[1]
    index1 = 0
    for idx, i in enumerate(x):
        if i >= findx:
            index1 = idx
            break

    dis = (y[index1-1])+((findx-x[index1-1])/(x[index1]-x[index1-1]))*(y[index1]-y[index1-1])

    vel = v*(y[index1]-y[index1-1])/(x[index1]-x[index1-1])

    return dis,vel

def BaseExitation2(t):
    
    findx = (t*v)-1.5
    x = raw_data[0]
    y = raw_data[1]
    index1 = 0
    for idx, i in enumerate(x):
        if i >= findx:
            index1 = idx
            break

    dis = (y[index1-1])+((findx-x[index1-1])/(x[index1]-x[index1-1]))*(y[index1]-y[index1-1])

    vel = v*(y[index1]-y[index1-1])/(x[index1]-x[index1-1])

    return dis,vel

Masses = [1,1,50,50,500,400,56.8]
a = 1.2
b = 1.2


sys = MultiDOF(Masses)

k1 = 65000
k2= 65000
k3 = 32000
k4 = 32000
k5 = 754500

c1 = 250
c2 = 250
c3 = 7500
c4 = 7500
c5 = 3840

sys.LiveK = [
[k1,0,-k1,0,0,0,0],
[0,k2,0,-k2,0,0,0],
[-k1,0,k1+k3,0,-k3,a*k3,0],
[0,-k2,0,k2+k4,-k4,-b*k4,0],
[0,0,-k3,-k4,k3+k4+k5,a*k3-b*k4,-k5],
[0,0,a*k3,-b*k4,a*k3-b*k4,(a**2)*k3+(b**2)*k4,0],
[0,0,0,0,-k5,0,k5]
]

sys.C = [
[c1,0,-c1,0,0,0,0],
[0,c2,0,-c2,0,0,0],
[-c1,0,c1+c3,0,-c3,a*c3,0],
[0,-c2,0,c2+c4,-c4,-b*c4,0],
[0,0,-c3,-c4,c3+c4+c5,a*c3-b*c4,-c5],
[0,0,a*c3,-b*c4,a*c3-b*c4,(a**2)*c3+(b**2)*c4,0],
[0,0,0,0,-c5,0,c5]
]

sys.SetInitialDisplacements([0,0,0,0,0,0,0])

sys.SetInitialVelocities([0,0,0,0,0,0,0])

sys.SetConstraints(BaseExitation1,1)
sys.SetConstraints(BaseExitation2,0)

t0 = 0
t1 = 10

sys.SolveSystem(t0,t1,res=10000)
#sys.plot(0,5)

t = np.linspace(t0,t1,sys.res)
plt.figure(1)
displacements = sys.y[:,0:sys.n]
plt.title('Resonse of vehicle with offsets')
plt.plot(t,displacements[:,0], label = 'Bump rear')
plt.plot(t,displacements[:,1], label = 'Bump front')
plt.plot(t,displacements[:,2]+0.1, label = 'Rear wheel')
plt.plot(t,displacements[:,3]+0.1, label = 'Front wheel')
plt.plot(t,displacements[:,4]+0.2, label = 'Car Body')
plt.plot(t,displacements[:,5], label = 'Car Body rotation(rad)')
plt.plot(t,displacements[:,6]+0.3, label = 'Driver head')
plt.xlabel('time')
plt.ylabel('x(t)')
plt.legend()
plt.grid()

plt.figure(2)
displacements = sys.y[:,0:sys.n]
plt.title('Resonse of vehicle with no offsets')
plt.plot(t,displacements[:,0], label = 'Bump rear')
plt.plot(t,displacements[:,1], label = 'Bump front')
plt.plot(t,displacements[:,2], label = 'Rear wheel')
plt.plot(t,displacements[:,3], label = 'Front wheel')
plt.plot(t,displacements[:,4], label = 'Car Body')
plt.plot(t,displacements[:,5], label = 'Car Body rotation(rad)')
plt.plot(t,displacements[:,6], label = 'Driver head')
plt.xlabel('time')
plt.ylabel('x(t)')
plt.legend()
plt.grid()

k = np.array(sys.C)
#print(k-k.T)
plt.show()
