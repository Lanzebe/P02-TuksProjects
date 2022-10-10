import sys
sys.path.append(r'C:\Users\ebers\Documents\GitHub Repositories\L02-PythonEngineeringLibrary')
from EngLib.Beam import *
import matplotlib.pyplot as plt
import numpy as np
import joblib

#### Problem Definition
b = 1.53
mass = 2.1
G_Loading = 6
g = 9.81
SF = 1.5


### Calculate a simulated distributed loading
r = 0.8

A = [[(b**3)/12,b],[(b**2)/4,(1-r)]]
B = [mass*g*G_Loading*SF,0]
a,c = np.linalg.solve(A,B)

def DistributedLoading(x):
    return a*x**2 + c


UpdatedIValue = joblib.load("C:\\Users\\ebers\\Documents\\GitHub Repositories\\P02-TuksProjects\\Y4-MLV420\\SparDesign\\IvalueData.txt")

### Loadcase 1



UPW_LC1 = Beam(b)
UPW_LC1.add_fixed_support('a',b/2)

for i in UpdatedIValue:
    UPW_LC1.define_I(i[1][0],i[1][1],i[0])


NumElemDist = 100


for i in range(NumElemDist):
    x0 = (i/NumElemDist)*b
    x1 = ((i+1)/NumElemDist)*b
    UPW_LC1.distributed_loading('dist_'+ str(i+1),[DistributedLoading(x0-b/2),DistributedLoading(x1-b/2)],[x0,x1])


UPW_LC1.solve_beam()

#UPW_LC1.plot_beam_diagrams( save=False)


### Loadcase 2:
LandingGearPos = 0.2
HeightDrop = 0.5
DropEnergy = g*HeightDrop*mass*SF
ForceLand = 5000

UPW_LC2 = Beam(b)
UPW_LC2.add_fixed_support('a',b/2)
UPW_LC2.force('Left_Land',ForceLand,b/2 - LandingGearPos)
UPW_LC2.force('Right_Land',ForceLand,b/2 + LandingGearPos)




for i in UpdatedIValue:
    UPW_LC2.define_I(i[1][0],i[1][1],i[0])



UPW_LC2.solve_beam()
UPW_LC2.plot_beam_diagrams(save=False)

Devlection1 = UPW_LC2.v(b/2 - LandingGearPos)
Devlection2 = UPW_LC2.v(b/2 + LandingGearPos)

print("Energy required to absorb:", DropEnergy)
print("Energy actualy absorbed", 0.5*ForceLand*Devlection1 + 0.5*ForceLand*Devlection2)

### Beam sizing


Sigma_yeild = 290E6     #https://www.matweb.com/search/DataSheet.aspx?MatGUID=57483b4d782940faaf12964a1821fb61&ckck=1
Height_ave = 0.0145
TaperRatio = 1.86

Height_r = 2*Height_ave*TaperRatio/(1+TaperRatio)
Height_t = 2*Height_ave/(1+TaperRatio)


ThicknessGrad = (Height_r-Height_t)/(b/2)

Thickness = 0.0004

SpanThicknessPos = np.linspace(0,b,100)

#print(SpanThicknessPos)

HeightArr = []
Width = []
IValue = []



###[i-value, [position start,end position]]
IvalueData = []
for idx,i in enumerate(SpanThicknessPos):
    
    if i <= b/2:
        Height = i*(ThicknessGrad) + Height_t
    if i > b/2:
        Height = ((SpanThicknessPos[idx])-b/2)*(-ThicknessGrad) + Height_r
        
    HeightArr.append(Height)
    Height_y = Height/2
    WidthsForLoadcases = [(12*Height_y*UPW_LC1.EIv_prime2(i)/Sigma_yeild)/(Height**3 - (Height - 2*Thickness)**3),(12*Height_y*UPW_LC2.EIv_prime2(i)/Sigma_yeild)/(Height**3 - (Height - 2*Thickness)**3)]

    Width_p = max(WidthsForLoadcases)
    if Width_p < 0.005:
        Width_p = 0.005

    IValue_p = ((Width_p * Height**3)/12)-((Width_p * (Height-2*Thickness)**3)/12)
    Width.append(Width_p)
    IValue.append(IValue_p)

for idx,i in enumerate(SpanThicknessPos):
    if idx !=0:
        if i <= b/2:
            IValue_p = ((Width[idx] * HeightArr[idx]**3)/12)-((Width[idx] * (HeightArr[idx]-2*Thickness)**3)/12)
        if i > b/2:
            IValue_p = ((Width[idx-1] * HeightArr[idx-1]**3)/12)-((Width[idx-1] * (HeightArr[idx-1]-2*Thickness)**3)/12)

        IvalueData.append([IValue_p,[SpanThicknessPos[idx-1],i]])

        
""" 
print('###########')
print(IvalueData)
print('###########') """

# write to txt
joblib.dump(IvalueData, "C:\\Users\\ebers\\Documents\\GitHub Repositories\\P02-TuksProjects\\Y4-MLV420\\SparDesign\\IvalueData.txt")

plt.figure(1)
plt.plot(SpanThicknessPos,Width)
plt.grid()

plt.figure(2)
plt.plot(SpanThicknessPos,IValue)
plt.ylabel('I-value [m^4]')
plt.grid()

""" plt.figure(3)
plt.plot(SpanThicknessPos,HeightArr)
plt.ylabel('Airfoil Thickness [m]')
plt.grid() """

plt.show()