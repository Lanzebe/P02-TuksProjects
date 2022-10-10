""" 
Weather for irene:
https://weatherspark.com/y/148604/Average-Weather-at-Pretoria-Irene-South-Africa-Year-Round
"""


import numpy as np
import matplotlib.pyplot as plt
from  AircraftSpecs import *





print("VoltageFull:             [V]:",VoltageFull)
print("VoltageDepeted           [V]:",VoltageDepeted)
print("BatteryEnergy            [J]:",BatteryEnergy)
print("Ranges                   [m]:",Range)
print("Average Minimum Speed  [m/s]:",AverageMinSpeed)
#print("Ideal Speed            [m/s]:",IdealSpeed)


Mass1 = (eta_engine*eta_prop * C_lod * BatteryEnergy)/(Range*g)
Mass2 = (eta_engine*eta_prop * C_lod * PowerActual)/(CruiseVSAveVelocity * AverageMinSpeed * g)

PowerToWeight = MaxEnginePower / Mass2    #(AverageMinSpeed*g)/(eta_engine*eta_prop*C_lod)

print("Mass1 if all energy is used    [kg]:", Mass1)
print("Mass2 if all max power limeted [kg]:", Mass2)
print("Power to Weight Ratio        [W/kg]:", PowerToWeight)


Sref = (2*Mass2*g)/(rho * CruiseVSAveVelocity * AverageMinSpeed**2 * ClAtAlpha)

b = (AR*Sref)**0.5
c = b/AR

print("Estimated Sref area: [m^2]",Sref)
print('Estimated Wingspan [m]', b)
print('Estimated Chord [m]', c)

plt.figure(1)
plt.plot(Laps,PowerToWeight)
plt.plot(Laps,110 * np.ones(len(Laps)))
plt.show()

