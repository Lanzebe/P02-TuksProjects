import numpy as np

### Physical Constants:
g = 9.81
rho = 1.0687

### Battery Specifications:
EnergyCapacityBattery = 2200  #mAh
VoltageFull = 12.6  ### 4.2 max cell voltage dws 4.2*3 = 12.6 V max
MinCellVol = 3.5  ### 3.2 min typical 3.5-3.7
VoltageDepeted = 3 * MinCellVol

VoltageUsable = VoltageFull - VoltageDepeted
VoltageAverage = (VoltageFull + VoltageDepeted)/2
BatteryEnergy = 60 * 60 * VoltageAverage * EnergyCapacityBattery / 1000

### Motor Specification
MaxEnginePower = 155 #155.4   ###max contineous power
PowerActual = 100

### Airplane Properties:
eta_engine = 0.5  ##0.65
eta_prop = 0.5    ##0.7
C_lod = 18.5  
ClAtAlpha = 0.866
CruiseVSAveVelocity = 1.25
AR = 10

### Mission specs:
TurningRadius = 50 # m worst case
Laps = np.array([3,4])
Range = (150 + 150 + np.pi* TurningRadius) * Laps
CompetitionTime = 5*60              ###s
PercentageTimeNotFlying = 0.3



AverageMinSpeed = Range/(CompetitionTime*(1-PercentageTimeNotFlying))
IdealSpeed = Range*PowerActual/BatteryEnergy