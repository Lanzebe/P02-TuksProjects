import numpy as np
import matplotlib.pyplot as plt
from  AircraftSpecs import *
from random import random


C_Dmin = 0.028
C_L_TO = 1.1
C_D_TO = 0.05

e = 1.78*(1-0.045*(AR**0.68)) - 0.64
k = 1/(np.pi * AR * e)
n = 1.2
mu = 0.04

WtoS = np.linspace(2,80,1000)


TtoW_Bank = []
TtoW_Cruise = []
TtoW_TO = []
TtoW_Climb = []

PtoW_Bank = []
PtoW_Cruise = []
PtoW_TO = []
PtoW_Climb = []



for i in AverageMinSpeed:
    TtoW_b = (0.5*rho*(i**2))*((C_Dmin/WtoS)+(k*((n/(0.5*rho*(i**2)))**2)*WtoS))
    TtoW_c = (0.5*rho*(i**2))*C_Dmin*(1/WtoS)+(k*(1/(0.5*rho*(i**2)))* WtoS)
    TtoW_to = ((i**2)/(2*g*15))+((C_D_TO*0.5*rho*i**2)/WtoS) + mu*(1-((C_L_TO*0.5*rho*i**2)/(WtoS)))
    TtoW_climb = (0.4/(0.75*i))+((C_Dmin*0.5*rho*i**2)/WtoS)+(k*(1/(0.5*rho*(i**2)))* WtoS)


    TtoW_Bank.append(TtoW_b)
    TtoW_Cruise.append(TtoW_c)
    TtoW_TO.append(TtoW_to)
    TtoW_Climb.append(TtoW_climb)

    PtoW_Bank.append(735.5*(i/(eta_engine*eta_prop*550))*TtoW_b)  ### W per lb
    PtoW_Cruise.append(735.5*(i/(eta_engine*eta_prop*550))*TtoW_c)  ### W per lb
    PtoW_TO.append(735.5*(i/(eta_engine*eta_prop*550))*TtoW_to)  ### W per lb
    PtoW_Climb.append(735.5*(i/(eta_engine*eta_prop*550))*TtoW_climb)  ### W per lb

    C_L_max = (1/(0.5*rho*(i*0.7)**2))*(50)
    print(C_L_max)


##Cruise Condition:

#PtoW_cruise = 

plt.figure()

plotlaps = [3]

for i in plotlaps:

    lableplot = 'Laps: ' + str(i)
    color = [random(),random(),random()]
    plt.title('4 Lap Constraint Diagam')
    """ plt.plot(WtoS,PtoW_Bank[i-1],color=color, label = lableplot)
    plt.plot(WtoS,PtoW_Cruise[i-1],color=color, linestyle='dotted')
    plt.plot(WtoS,PtoW_TO[i-1],color=color, linestyle='dashed')
    plt.plot(WtoS,PtoW_Climb[i-1],color=color, linestyle='dashdot') """

    color_r = [1,0,0]
    color_b = [0,0,1]
    plt.plot(WtoS,PtoW_Bank[i-1],color=color_b, linestyle='solid', label = 'Bank Power')
    plt.plot(WtoS,PtoW_Cruise[i-1],color=color_b, linestyle='dotted', label = 'Cruise Power')
    plt.plot(WtoS,PtoW_TO[i-1],color=color_b, linestyle='dashed', label = 'Take Off Power')
    plt.plot(WtoS,PtoW_Climb[i-1],color=color_b, linestyle='dashdot', label = 'Climb Power')

plt.xlabel('W/S [N/m^2]')
plt.ylabel('P/W [W/N]')
plt.grid()
plt.legend()
plt.show()