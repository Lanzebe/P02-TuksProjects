from engineering_lib import Beam
import numpy as np

a = 0.4
b = 1
c = 0.3
d = 0.067
P = 80*10**3
N = 2000 #rpm
Db = 0.2
Dc = 0.25


L = a+b+c
Pnew = 1.5*P
T = 30*Pnew/(N*np.pi)
Fb = T*10/(4*Db)
Fc = T*10/(4*Dc)
I = (np.pi*d**4)/64

print('L:',L,'T:',T,'Fb:',Fb,'Fc:',Fc)


h = 0.32

Bxy = Beam(L, E=200*10**9)
Bxy.define_I(0,L,I)
Bxy.add_pinned_support('a',0)
Bxy.add_pinned_support('e',a+h)
Bxy.add_pinned_support('d',L)
Bxy.force('b',-1.2*Fb,a)
Bxy.solve_beam()
Bxy.plot_beam_diagrams(save=True)

