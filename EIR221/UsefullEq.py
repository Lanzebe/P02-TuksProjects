import math
import numpy as np
import cmath

# conversion of degrees to rad
theta = math.radians(21.8) 



### Linear complex equatoins
A = np.array([[1+1.5j,2.5j],
              [11,15]])

B = np.array([[20],
              [0]])


def cramer_solve2(A,B):
    C1 = np.array(A,copy=True)
    C2 = np.array(A,copy=True)
    print('A =\n',A,'\nB =\n',B)
    for i in range(0,len(A)):
        C1[i,0] = B[i,0] 
        C2[i,1] = B[i,0]
    print('C1 =\n',C1,'\nC2 =\n',C2)
    det_A = np.linalg.det(A)
    print('det A =',np.round(det_A,3))
    det_C1 = np.linalg.det(C1)
    print('det_C1 =',np.round(det_C1,3))
    det_C2 = np.linalg.det(C2)
    print('det_C2 =',np.round(det_C2,3))
    print('\nx1 = {} \nx2 = {}'.format(np.round(det_C1/det_A,3),np.round(det_C2/det_A,3)))
    
    A1,theta1 =cmath.polar(det_C1/det_A)
    print('\nx1 = {} /_{}deg'.format(np.round(A1,3),np.round(math.degrees(theta1),3)))
    A2,theta2 =cmath.polar(det_C2/det_A)
    print('x2 = {} /_{}deg'.format(np.round(A2,3),np.round(math.degrees(theta2),3)))

cramer_solve2(A,B)

print("\nCHECKING:")
x1,x2 = np.linalg.solve(A,B)
print('x1 = {} \nx2 = {}'.format(x1[0],x2[0]))




###
def print_polar(n):
    print(cmath.polar(n)[0],'/__',360*cmath.polar(n)[1]/(2*cmath.pi),'deg')


def parallel(Zlist):
    ans = 0
    for i in Zlist:
        ans = ans + i**-1
    return ans**-1

def ComplexPolar(z,thetaDeg):
    return z*np.cos(math.radians(thetaDeg)) + z*np.sin(math.radians(thetaDeg))*1j