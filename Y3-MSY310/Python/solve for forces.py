import numpy as np
from scipy.optimize import minimize
import one_D_beam as BM

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


#Initail CALCS
P_new = 1.5 * P #W
T = 30*P_new/(np.pi*N)
F_B = 10*T/(4*D_B)
F_C= 10*T/(4*D_C)


#print(T, F_B, F_C)


## initial conditions
h = 0.15
x0  = np.array([10,10,10,10,10,10])   #np.array([F_Ay ,F_Dy ,F_Ey ,C1 ,C2 ,C3])

def delta_d (x0):
    x = 1.7
    return ((x0[0]*(x**3)/6)-(1.2*F_B*((x-0.4)**3)/6)+(x0[2]*((x-0.4-h)**3)/6)+(x0[1]*((x-1.7)**3)/6)+(x0[3]*(x**2)/2)+(x0[4]*(x))+x0[5])

def delta_a (x0):
    x = 0
    return ((x0[0]*(x**3)/6)+(x0[3]*(x**2)/2)+(x0[4]*(x))+x0[5])

def delta_e (x0):
    x = 0.4+h
    return ((x0[0]*(x**3)/6)-(1.2*F_B*((x-0.4)**3)/6)+(x0[3]*(x**2)/2)+(x0[4]*(x))+x0[5])

def M_a (x0):
    x = 0
    return ((x0[3]))

def M_d (x0):
    x = 1.7
    return ((x0[0]*(x))-(1.2*F_B*(x-0.4))+(x0[2]*(x-0.4-h))+(x0[1]*(x-1.7))+(x0[3]))

def sum_forces_y(x0):
    return x0[0]-(1.2*F_B)+x0[1]+x0[2]

# show initial objective
print('Initial Objective: ' + str(delta_d(x0)))


# optimize
b = (-10**9, 10**9)
bnds = (b, b, b, b, b, b)
con1 = {'type': 'eq', 'fun': delta_d}
con2 = {'type': 'eq', 'fun': delta_a}
con3 = {'type': 'eq', 'fun': delta_e}
con4 = {'type': 'eq', 'fun': M_a}
con5 = {'type': 'eq', 'fun': M_d}
con6 = {'type': 'eq', 'fun': sum_forces_y}
cons = ([con1,con2,con3,con4,con5,con6])
myoptions={'disp':True}
solution = minimize(delta_d ,x0 ,method='SLSQP',options = myoptions,bounds=bnds ,constraints=cons)

x = solution.x

# show final objective
print('Final Objective: ' + str(delta_d(x)))

# print solution
print('Solution')
print('x = ' + str(x))

if (abs(delta_a(x))<10**(-6))and(abs(delta_d(x))<10**(-6))and(abs(delta_e(x))<10**(-6))and(abs(M_a(x))<10**(-6))and(abs(M_d(x))<10**(-6))and(abs(sum_forces_y(x))<10**(-6)):
    print('Good solution')
else:
    print('delta_a=', delta_a(x))
    print('delta_d=', delta_d(x))
    print('delta_e=', delta_e(x))
    print('M_a=',M_a(x))
    print('M_d=',M_d(x))
    print('sum_forces_y=',sum_forces_y(x))

print('sum_moments_yA=',((0.4+h)*x[2])+(1.7*x[1])-(0.48*F_B))

Shaft_y = BM.Beam(1.7)

Shaft_y.add_point_load(x[0],0)
Shaft_y.add_point_load(x[1],1.7)
Shaft_y.add_point_load(x[2],0.4+h)
Shaft_y.add_point_load(-1.2*F_B,0.4)

Shaft_y.plot_beam_data()