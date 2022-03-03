

from fem_main_program import *

U, ReactionForces, VonMises, SXX, SYY, SXY = launch_fem('q8_element_pull_planestress',MagFac=1000)

print('Work done: {}'.format(0.5*np.dot(F.transpose(),U)))
print(U)