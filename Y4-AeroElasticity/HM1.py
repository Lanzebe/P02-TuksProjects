import numpy as np
import matplotlib.pyplot as plt

wn = 2*np.pi
w = wn*np.array([0.5,1.00001,2])


t = np.linspace(0,20,1000)

A1 = ((-1)/(wn**2 - w**2))
#A2 = 0

freq = 1
x1 = A1[freq]*np.exp(wn*t*(1j)) + (1/(wn**2 - w[freq]**2))*np.exp(-w[freq]*t*(1j))


plt.figure()
plt.plot(t,np.real(x1))
plt.show()