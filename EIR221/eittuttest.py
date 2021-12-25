from engineering_lib import *

#3.49
C1 = RLCSF_Series_Circuit(4,0.45,10*10**(-3))
C1.initial_conditions([0,0],[0,-320/3])
C1.plot(0,2)