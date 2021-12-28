import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from IPython.display import display


def plot2d(u,v,x0,x1,y0,y1,n=50):
    uL = sp.lambdify([x,y],u,'numpy')
    vL = sp.lambdify([x,y],v,'numpy')
    X = np.linspace(x0,x1,n)
    Y = np.linspace(y0,y1,n)

    Xp,Yp = np.meshgrid(X,Y)

    plt.figure()
    plt.quiver(Xp,Yp,uL(Xp,Yp),vL(Xp,Yp))
    plt.show()

def plot3d(u,v,w,x0,x1,y0,y1,z0,z1,n=50):
    uL = sp.lambdify([x,y,z],u,'numpy')
    vL = sp.lambdify([x,y,z],v,'numpy')
    wL = sp.lambdify([x,y,z],w,'numpy')
    X = np.linspace(x0,x1,n)
    Y = np.linspace(y0,y1,n)
    Z = np.linspace(z0,z1,n)

    Xp,Yp,Zp = np.meshgrid(X,Y,Z)

    plt.figure().add_subplot(projection='3d')
    scale = 50
    plt.quiver(Xp,Yp,Zp,uL(Xp,Yp,Zp)*scale,vL(Xp,Yp,Zp)*scale,wL(Xp,Yp,Zp)*scale)
    plt.show()

def StrainMatrixFromDisplacement3d(u,v,w,x,y,z):
    displacement_vector = [u,v,w] 
    spatial_coordinates = [x, y, z]
    displacement_gradient = sp.derive_by_array(displacement_vector, spatial_coordinates)
    return sp.simplify(0.5*(displacement_gradient.transpose() + displacement_gradient))

def StressTensorFromStrain3d(StrainMatrix ,D = 1):
    if D == 1:
        D = E/((1 + nu)*(1 - 2*nu))*sp.Matrix([[1-nu,nu,nu,0,0,0],[nu,1-nu,nu,0,0,0],[nu,nu,1-nu,0,0,0],[0,0,0,0.5-nu,0,0],[0,0,0,0,0.5-nu,0],[0,0,0,0,0,0.5-nu]])
    strain_vector = sp.Matrix([StrainMatrix[0,0],StrainMatrix[1,1],StrainMatrix[2,2],2*StrainMatrix[1,2],2*StrainMatrix[0,2],2*StrainMatrix[0,1]])
    
    return sp.simplify(D*strain_vector)

def StrainMatrixFromDisplacement2d(u,v,x,y):
    displacement_vector = [u,v] 
    spatial_coordinates = [x, y]
    displacement_gradient = sp.derive_by_array(displacement_vector, spatial_coordinates)
    return sp.simplify(0.5*(displacement_gradient.transpose() + displacement_gradient))

def StressTensorFromStrain2d(StrainMatrix ,D = 1):
    if D == 1:
        D = E/((1 + nu)*(1 - 2*nu))*sp.Matrix([[1-nu,nu,0],[nu,1-nu,0],[0,0,0.5-nu]])
    strain_vector = sp.Matrix([StrainMatrix[0,0],StrainMatrix[1,1],2*StrainMatrix[0,1]])
    return sp.simplify(D*strain_vector)




#Testing StrainMatrixFromDisplacement
rho, g, nu, E, W, B, H, x, y, z = sp.symbols('rho, g, nu, E, W, B, H, x, y, z')

u = -nu*rho*g*x*z/E
v = -nu*rho*g*y*z/E
w = rho*g*z**2/(2*E) + nu*rho*g*(x**2+y**2)/(2*E) - rho*g*H**2/(2*E)

D = E/((1 + nu)*(1 - 2*nu))*sp.Matrix([[1-nu,nu,nu,0,0,0],[nu,1-nu,nu,0,0,0],[nu,nu,1-nu,0,0,0],[0,0,0,0.5-nu,0,0],[0,0,0,0,0.5-nu,0],[0,0,0,0,0,0.5-nu]])

StrainMatrix = StrainMatrixFromDisplacement3d(u,v,w,x,y,z)
print('Strain Matrix:')
display(StrainMatrix)

StressMatrix = StressTensorFromStrain3d(StrainMatrix,D = D)
print('Stress Matrix:')
display(StressMatrix)

problemspec = {rho: 7600, g: 9.81, nu: 0.3, E: 210E9, W: 0.1, B: 0.1, H: 100}

uProb = u.subs(problemspec)
vProb = v.subs(problemspec)
wProb = w.subs(problemspec)

x0 = problemspec[B]
x1 = -problemspec[B]
y0 = problemspec[W]
y1 = -problemspec[W]
z0 = 0
z1 = problemspec[H]


plot3d(uProb,vProb,wProb,x0,x1,y0,y1,z0,z1,n=10)
