import cmath
import numpy as np
import math


def print_polar(n):
    print(cmath.polar(n)[0],'/__',360*cmath.polar(n)[1]/(2*cmath.pi),'deg')

def parallel(Zlist):
    ans = 0
    for i in Zlist:
        ans = ans + i**-1
    return ans**-1

def ComplexPolar(z,thetaDeg):
    return z*np.cos(math.radians(thetaDeg)) + z*np.sin(math.radians(thetaDeg))*1j


Z1 = 13-5j
Z2 = 10
Z3 = 8+6j
Z4 = 12

V = ComplexPolar(90,30)

Zth = parallel([Z1,Z2]) + parallel([Z3,Z4])

Vab = V*((Z2/(Z2+Z1))-(Z3/(Z3+Z4)))

print('Rth:',Zth)
print('Vth:',Vab)
print_polar(Vab)

In = Vab/Zth

print_polar(In)