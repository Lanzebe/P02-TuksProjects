import cmath
import math

def print_polar(n):
    print(cmath.polar(n)[0],'/__',360*cmath.polar(n)[1]/(2*cmath.pi),'deg')


omega = 60
C1 = 100E-6
C2 = 300E-6
L1 = 400E-3
L2 = 200E-3
R1 = 12
R2 = 20

Zc1 = 1/(omega*C1*(1j))
Zc2 = 1/(omega*C2*(1j))
Zl1 = omega*L1*(1j)
Zl2 = omega*L2*(1j)
Zr1 = R1
Zr2 = R2

Vs1 = -7.7646 - 28.978j
Vs2 = 19.319 + 5.1764j

#print_polar(Vs1)
#print_polar(Vs2)


V2 = ((Vs2/Zr2)-(Vs1/(Zc1+Zr1)))/((1/Zc1)+(1/Zl2)+(1/Zc2)+(1/Zr2))

V1 = ((Vs1/(Zc1+Zr1)) + (V2/Zl2))/((1/(Zc1+Zr1))+(1/Zl2)+(1/Zl1))

M1 = 10/((1/(Zc1+Zr1))+(1/Zl2)+(1/Zl1))
M2 = 20/(-1/Zl2)



k1 = (-1/Zl2)*M1
k2 = (Vs1/(Zc1+Zr1))*M1
k3 = ((1/Zl2)+(1/Zc2)+(1/Zr2))*M2
k4 = (Vs2/Zr2)*M2


""" print('Zl1',Zl1)
print('Zl2',Zl2)
print('Zc1',Zc1)
print('Zc2',Zc2)


print_polar(k1)
print_polar(k2)
print_polar(k3)
print_polar(k4) """



Z1 = Zr1 + Zc1
I1 = Vs1/Z1

#print(Z1)
#print(I1)


I4 = (Zr2*Vs2*(((1/Zr2)+(1/Zc2))**-1))/(Zl2+((((1/Zr2)+(1/Zc2))**-1)))
Z2 = (Zl2+((((1/Zr2)+(1/Zc2))**-1)))

V3 = (I1*I4)*(((1/Z1)+(1/Z2))**-1)
Z3 = (((1/Z1)+(1/Z2))**-1)


V0 = V3*Zl1/(Zl1+Z3)
print(V0)
#print(Z3)