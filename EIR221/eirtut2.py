from engineering_lib import *

#3.49
C1 = RLCS_Parallel_Circuit(5,5,50*10**(-3),9)
C1.initial_conditions([0,12],[0,0])
C1.plot(0,5)

#8.14

#R = 10
#L = 4
#C = 0.25

#DE3 = LSODECC(1,R/L,1/(L*C))
#DE3.initial_conditions([0,2],[0,-5])
#DE3.plot(0,1)


#8.9
#DE2 = LSODECC(1,10,25)
#DE2.initial_conditions([0,10],[0,0])
#DE2.plot(0,1)


#8.7

R = 20*10**3
L = 0.2*10**(-3)
C = 5*10**(-6)

#DE1 = LSODECC(1,R/L,1/(L*C))
#DE1.initial_conditions([0,1],[0,1])
#DE1.plot(0,1)