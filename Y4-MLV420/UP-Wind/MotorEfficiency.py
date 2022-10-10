""" 
For Motor efficiency go to :
https://www.radiocontrolinfo.com/brushless-motor-efficiency/
https://www.bavaria-direct.co.za/constants/#:~:text=Rm%20(or%20Ri)%20is%20the,value%20would%20be%20the%20Rm.
"""

import numpy as np
import matplotlib.pyplot as plt

#I,V = sp.symbols('I,V')
I = 14 #np.linspace(1,19,100)
V = 9.9

I0 = (0.7/8) * 9.9  ## 0.7 a @ 8V
Rm = 0.2            ## Equivalent to Ri


PowerIn = V * I
CuLosses = Rm * I**2
FeLosses = V * I0
PowerOut = (PowerIn - CuLosses - FeLosses) 
Eff = PowerOut / PowerIn


print('PowerIn:',PowerIn)
print('Eff:',Eff)
print('CuLosses:',CuLosses)
print('FeLosses:',FeLosses)
print('PowerOut:',PowerOut)

""" plt.figure()
plt.plot(I,Eff)
plt.grid()
plt.show() """