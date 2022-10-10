import numpy as np
import sympy as sp
from sympy.plotting import plot
import matplotlib.pyplot as plt


x_PluckLocation = np.linspace(0.001,0.5,10)
y_MaxDisplacement = 1

N_Approximations = 2


q = []

x = sp.symbols('x')

""" for x_Pluck in x_PluckLocation:
    DisplacementFun =  sp.Piecewise( ((1/x_Pluck)*x , x<x_Pluck), (((-1/(1-x_Pluck)*x) + (1/(1-x_Pluck)) , x >=x_Pluck)) )

    q_vec = []
    for n in range(N_Approximations):
        #print(sp.integrate((DisplacementFun*sp.sin((n+1)*sp.pi*x)), (x, 0, 1)).evalf())
        q_vec.append(sp.integrate((DisplacementFun*sp.sin((n+1)*sp.pi*x)), (x, 0, 1)))


    q.append(q_vec)

print(q) """


def GenerateDisplacementFun(PluckLocation,N):
    DisplacementFun =  sp.Piecewise( ((1/PluckLocation)*x , x<PluckLocation), (((-1/(1-PluckLocation)*x) + (1/(1-PluckLocation)) , x >=PluckLocation)) )
    Func = 0
    for n in range(N):
        Func = Func + (sp.integrate((DisplacementFun*sp.sin((n+1)*sp.pi*x)), (x, 0, 1))/sp.integrate((sp.sin((n+1)*sp.pi*x))**2, (x, 0, 1))) * sp.sin((n+1)*sp.pi*x)
    return Func


def EigenValues(PluckLocation,N):
    DisplacementFun =  sp.Piecewise( ((1/PluckLocation)*x , x<PluckLocation), (((-1/(1-PluckLocation)*x) + (1/(1-PluckLocation)) , x >=PluckLocation)) )
    EigValues = []
    for n in range(N):
        EigValues.append(sp.integrate((DisplacementFun*sp.sin((n+1)*sp.pi*x)), (x, 0, 1))/sp.integrate((sp.sin((n+1)*sp.pi*x))**2, (x, 0, 1)))
    return EigValues

Pluck = 0.2
ys = sp.Piecewise( ((1/0.2)*x , x<0.2), (((-1/(1-0.2)*x) + (1/(1-0.2)) , x >=0.2)) )
ys1 = sp.Piecewise( ((1/Pluck)*x , x<Pluck), (((-1/(1-Pluck)*x) + (1/(1-Pluck)) , x >=Pluck)) )#GenerateDisplacementFun(0.1, 10)
ys2 = GenerateDisplacementFun(0.2, 10)
ys3 = GenerateDisplacementFun(0.3, 10)
ys4 = GenerateDisplacementFun(0.4, 10)
ys5 = GenerateDisplacementFun(0.5, 10)

EigPlot = EigenValues()

p1 = plot(ys1,(x,0,1),show=False)
p2 = plot(ys2,(x,0,1),show=False)
p3 = plot(ys3,(x,0,1),show=False)
p4 = plot(ys4,(x,0,1),show=False)
p5 = plot(ys5,(x,0,1),show=False)

pp5 = plot()

p1.extend(p2)
p1.extend(p3)
p1.extend(p4)
p1.extend(p5)
p1.show()