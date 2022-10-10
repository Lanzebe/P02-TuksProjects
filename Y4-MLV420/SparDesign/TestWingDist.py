import matplotlib.pyplot as plt
import numpy as np


b = 1.95
r = 1
mass = 1.9
g = 4
SF = 1.5

A = [[(b**3)/12,b],[(b**2)/4,(1-r)]]
B = [mass*g*SF,0]
a,c = np.linalg.solve(A,B)


xs = np.linspace(-b/2,b/2,1000)

plt.figure()
plt.grid()
plt.ylim([0,c*1.2])
plt.plot(xs, a*xs**2 + c)
plt.show()