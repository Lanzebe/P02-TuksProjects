from ansys.mapdl.core import launch_mapdl
import re
import numpy as np

path = 'c:\\Users\\ebers\\Simulations\\MRN412\\Testing Scripting'



lenght = 0.1077
width = 0.108
thickness = 0.00203
YoungsModulus = 73.1E9
Density = 2780
PoisonRatio = 0.33



mapdl = launch_mapdl(run_location=path, jobname='TestingScripting',override=True)

Core_data = [lenght,width,thickness,YoungsModulus,Density,PoisonRatio]

## General Setup
mapdl.clear()
mapdl.prep7()
mapdl.units("SI")  # SI - International system (m, kg, s, K).

##Geometry
mapdl.block(0,lenght,0,width,0,thickness)
#mapdl.vplot('all')


## Material properties
mapdl.lesize('ALL',thickness, layer1=1) ### check seeds.

mapdl.mp('ex',1,YoungsModulus)
mapdl.mp('nuxy',1,PoisonRatio)
mapdl.mp('dens',1,Density)

mapdl.et(1, 186)
mapdl.esize(thickness)
mapdl.vsweep("all")


mapdl.mshape(1,'3D')
mapdl.mshkey(0)
mapdl.vmesh('all')
#mapdl.eplot()

mapdl.run("/SOLU")
mapdl.antype("MODAL")
output = mapdl.modal_analysis(nmode=8, freqb=1)



mapdl.post1()
for i in range(8):
    resp = mapdl.set(lstep='1',sbstep=str(i+1))
    Core_data.append(float(resp[153:159]))


print(Core_data)
mapdl.exit()