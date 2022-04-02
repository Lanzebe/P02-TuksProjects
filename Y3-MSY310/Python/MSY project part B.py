import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#Knowns
P = 80 * 10**3  #W
N = 2000        #rpm
d = 0.067       #m
D_B = 0.2       #m
D_C = 0.25      #m
a = 0.4         #m
b = 1           #m
c = 0.3         #m
sigma_allow = 80 * 10**6    #Pa
tau_allow = 40 * 10**6      #Pa
delta_allow = 0.002         #m
theta_allow = 3             #deg

E = 200 * 10**9 #Pa
I = np.pi*(d**4)/64 #m^4

#Initail CALCS
P_new = 1.5 * P #W
T = 30*P_new/(np.pi*N)
F_B = 10*T/(4*D_B)
F_C= 10*T/(4*D_C)

print('torque in shaft ', T)
#required functions
def mac(x,a):
    if x>=a:
        return x-a
    else:return 0

def macv(x,a):
    if x>a:
        return 1
    else:return 0

## initial conditions
Fay = F_B*(1.2-((1.2*a)/(a+b+c)))
Fdy = ((F_B*1.2*a)/(a+b+c))

Faz = F_C*(((a+b)*1.2/(a+b+c))-1.2)
Fdz = -(((a+b)*1.2*F_C)/(a+b+c))


print('Fay',Fay)
print('Fdy',Fdy)
print('Faz',Faz)
print('Fdz',Fdz)
Y  = np.array([0,0,0])   #np.array([F_Ay ,F_Dy ,F_Ey ,C1 ,C2 ,C3])
Z  = np.array([0,0,0])   #np.array([F_Az ,F_Dz ,F_Ez ,k1 ,k2 ,k3])


###XY_Plane Functions
def shearF_y(x0, *args):
    x = args[0]
    #print(x)
    residual = args[1]
    #print(residual)
    return (((Fay*(macv(x,0)))-(1.2*F_B*((macv(x,a))))+(Fdy*(macv(x,a+b+c))))-residual)

def M_z (x0, *args):
    x = args[0]
    residual = args[1]
    return ((Fay*(x))-(1.2*F_B*mac(x,a))+(Fdy*mac(x,a+b+c))+(x0[0])-residual)

def theta_y(x0, *args):
    x = args[0]
    #print(x)
    residual = args[1]
    #print(residual)
    return (((Fay*(x**2)/2)-(0.6*F_B*((mac(x,a)**2)))+(Fdy*(mac(x,a+b+c)**2)/2)+(x0[0]*(x**1)/1)+(x0[1]))-residual)

def deflection_y(x0, *args):
    x = args[0]
    #print(x)
    residual = args[1]
    #print(residual)
    return (((Fay*(x**3)/6)-(0.2*F_B*((mac(x,a)**3)))+(Fdy*(mac(x,a+b+c)**3)/6)+(x0[0]*(x**2)/2)+(x0[1]*(x))+x0[2])-residual)

def sum_forces_y(x0):
    return Fay-(1.2*F_B)+Fdy

def sum_moments_z(x0):
    return (-1.2*F_B*a)+((a+b+c)*Fdy)

###XZ_Plane Functions
def shearF_z(x0, *args):
    x = args[0]
    #print(x)
    residual = args[1]
    #print(residual)
    return (((Faz*(macv(x,0)))+(1.2*F_C*((macv(x,a+b))))+(Fdz*(macv(x,a+b+c))))-residual)

def M_y (x0, *args):
    x = args[0]
    residual = args[1]
    return ((Faz*(x))+(1.2*F_C*mac(x,a+b))+(Fdz*mac(x,a+b+c))+(x0[0])-residual)

def theta_z(x0, *args):
    x = args[0]
    #print(x)
    residual = args[1]
    #print(residual)
    return (((Faz*(x**2)/2)+(0.6*F_C*((mac(x,a+b)**2)))+(Fdz*(mac(x,a+b+c)**2)/2)+(x0[0]*(x**1)/1)+(x0[1]))-residual)

def deflection_z(x0, *args):
    x = args[0]
    #print(x)
    residual = args[1]
    #print(residual)
    return (((Faz*(x**3)/6)+(0.2*F_C*((mac(x,a+b)**3)))+(Fdz*(mac(x,a+b+c)**3)/6)+(x0[0]*(x**2)/2)+(x0[1]*(x))+x0[2])-residual)

def sum_forces_z(x0):
    return Faz+(1.2*F_C)+Fdz

# Solve forces for x-y plane
bb = (-10**9, 10**9)
bnds = (bb, bb, bb)
con1 = {'type': 'eq', 'fun': deflection_y, 'args':(a+b+c,0)}
con2 = {'type': 'eq', 'fun': deflection_y, 'args':(0,0)}
con3 = {'type': 'eq', 'fun': M_z,'args':(0,0)}
cons = ([con1,con2,con3])
myoptions={'disp':True}

solution = minimize(deflection_y ,Y ,method='SLSQP',options = myoptions,bounds=bnds ,constraints=cons,args=(0,0))
Y = solution.x


# print solution
print('SolutionY')
print('x = ' + str(Y))
print('delta_a = ', deflection_y(Y,0,0))
print('delta_d = ', deflection_y(Y,a+b+c,0))
print('M_za = ', M_z(Y,0,0))
print('M_zd = ', M_z(Y,a+b+c,0))
print('Sum of forces = ', sum_forces_y(Y))
print('Sum of moments in z = ', sum_moments_z(Y))


# Solve forces for x-z plane
con1 = {'type': 'eq', 'fun': deflection_z, 'args':(a+b+c,0)}
con2 = {'type': 'eq', 'fun': deflection_z, 'args':(0,0)}
con3 = {'type': 'eq', 'fun': M_y,'args':(0,0)}
cons = ([con1,con2,con3])

solution = minimize(M_z ,Z ,method='SLSQP',options = myoptions,bounds=bnds ,constraints=cons,args=(0,0))
Z = solution.x

# print solution
print('SolutionZ')
print('x = ' + str(Z))
print('delta_a = ', deflection_z(Z,0,0))
print('delta_d = ', deflection_z(Z,a+b+c,0))
print('M_za = ', M_y(Z,0,0))
print('M_zd = ', M_y(Z,a+b+c,0))
print('Sum of forces = ', sum_forces_z(Z))
#print('Sum of moments in z = ', sum_moments_z(Z))


x = np.linspace(0, a+b+c, 1000)

#### XY plane
def V_y_plot(x):
    y = []
    for i in x:
        y.append(shearF_y(Y,i,0))
    return y

def mom_z_plot(x):
    y = []
    for i in x:
        y.append(M_z(Y,i,0))
    return y

def thetay_plot(x):
    y = []
    for i in x:
        y.append(theta_y(Y,i,0))
    return y

def def_y_plot(x):
    y = []
    for i in x:
        y.append(deflection_y(Y,i,0))
    return y

V_y = V_y_plot(x)
Mz = mom_z_plot(x)
theta_y = thetay_plot(x)
delta_y = def_y_plot(x)

### XZ Plane
def V_z_plot(x):
    y = []
    for i in x:
        y.append(shearF_z(Z,i,0))
    return y

def mom_y_plot(x):
    y = []
    for i in x:
        y.append(M_y(Z,i,0))
    return y

def thetaz_plot(x):
    y = []
    for i in x:
        y.append(theta_z(Z,i,0))
    return y

def def_z_plot(x):
    y = []
    for i in x:
        y.append(deflection_z(Z,i,0))
    return y

V_z = V_z_plot(x)
My = mom_y_plot(x)
theta_z = thetaz_plot(x)
delta_z = def_z_plot(x)

###
plt.figure(1)
rows = 4
columns = 3
grid = plt.GridSpec(rows, columns, wspace = .25, hspace = .5)
### 
plt.subplot(grid[0])
plt.plot(x,V_y,'b-')
plt.title('Shear force(N) in XY-Plane')
plt.grid()

plt.subplot(grid[3])
plt.plot(x,Mz,'b-')
plt.title('Moments(Nm) about in XY-Plane')
plt.grid()

plt.subplot(grid[6])
plt.plot(x,theta_y,'b-')
plt.title('EIv_prime(x) in XY-Plane')
plt.grid()

plt.subplot(grid[9])
plt.plot(x,delta_y,'b-')
plt.title('EIv(x) in XY-Plane')
plt.grid()
###
plt.subplot(grid[1])
plt.plot(x,V_z,'g-')
plt.title('Shear force(N) in XZ-Plane')
plt.grid()

plt.subplot(grid[4])
plt.plot(x,My,'g-')
plt.title('Moments(Nm) about in XZ-Plane')
plt.grid()

plt.subplot(grid[7])
plt.plot(x,theta_z,'g-')
plt.title('EIv_prime(x) in XZ-Plane')
plt.grid()

plt.subplot(grid[10])
plt.plot(x,delta_z,'g-')
plt.title('EIv(x) in XZ-Plane')
plt.grid()
###
plt.subplot(grid[2])
plt.plot(x,np.power((np.power(V_z,2)+np.power(V_y,2)),0.5),'c-')
plt.title('Superpositioned combined Shear force(N)')
plt.grid()

plt.subplot(grid[5])
plt.plot(x,np.power((np.power(My,2)+np.power(Mz,2)),0.5),'c-')
plt.title('Superpositioned combined Moments(Nm)')
plt.grid()

plt.subplot(grid[8])
plt.plot(x,np.power((np.power(theta_y,2)+np.power(theta_z,2)),0.5),'c-')
plt.title('Superpositioned combined EIv_prime(x)')
plt.grid()

plt.subplot(grid[11])
plt.plot(x,np.power((np.power(delta_y,2)+np.power(delta_z,2)),0.5),'c-')
plt.title('Superpositioned combined EIv(x)')
plt.grid()

plt.show()