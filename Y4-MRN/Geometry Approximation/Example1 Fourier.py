import cmath
from math import pi
import matplotlib.pyplot as plt
import numpy as np
from FourierApprox import *




t0 = 0
t1 = 1

def approx_fun(t):
    ans = 0
    #return t**2
    #return np.sin(t)
    #return np.sin(3*(1-np.cos(t)))
    
    """ if t>=0 and t<1:
        ans = np.sin(t)
    else:
        ans = -0.5*t**2
    
    return ans*(-np.cos(np.pi*t)+1) """

    par = -(t-0)*(t-2)
    return np.sin(20*par)
    



series = FourierSeries(t0,t1,N_max_series=10,Integration_n=5000)
series.calc_cs(approx_fun)


t = np.linspace(t0,t1,10000)

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

plt.figure()
#plt.plot(t,X1)
plt.plot(t,X2)
plt.plot(t,y)
plt.grid()
""" for i in range(10):
    t = t0+ i*((t1-t0)/10)
    ans = series.ouput(t)
    plt.scatter(ans.real,ans.imag) """
plt.show()