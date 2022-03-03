"""

Author of code: Stephan Schmidt

This code is intended for and provided to Mechanical and Aeronautical Engineering Students at the University of Pretoria (UP) as part of the course MOW 323 - Simulation-based Design.

Please contact stephan.schmidt@up.ac.za if you have questions regarding the use of the code.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
__version__ = "20211007"

def newmark_gn22(t,M,C,K,x0,v0,force_function,beta1=0.5,beta2=0.5,calc_energy=False,verbose=False):
    import numpy as np
    import time
    """
    This function uses the Newmark GN22 method for solving a second order differential equation in the form:
    M x'' + C x' + K x = F(t)


    Parameters
    ----------------
    t: numpy.array
        The time vector corresponding to the time of the simulation.

    M: numpy.array
        The mass matrix of the system

    C: numpy.array
        The viscous damping matrix of the system

    K: numpy.array
        The stiffness matrix of the system

    x0: numpy.array
        The initial displacement of the system

    v0: numpy.array
        The initial velocity of the system

    force_function: function(time)
        The force function that takes time as an input and returns the force applied to each degree of freedom.

    beta1: float
        The beta1 parameter of the algorithm.

    beta2: float
        The beta2 parameter of the algorithm.

    calc_energy: bool
        Specifies whether the potential and kinetic energy of the system should be calculated (True) or not (False)

    Returns
    ----------------

    dict
        A dictionary is returned with the following arguments:
            - X: numpy.array
                    The displacement response at each time step.
            - V: numpy.array
                    The velocity response at each time step.
            - A: numpy.array
                    The acceleration response at each time step.

    """
    a0 = np.linalg.solve(M,force_function(t[0]))
    X = np.zeros([t.shape[0],M.shape[0]])
    V = np.zeros(X.shape)
    A = np.zeros(X.shape)
    X[0,:] = np.array(x0)
    V[0,:] = np.array(v0)
    A[0,:] = np.array(a0)

    if calc_energy == True:
        Pe = np.zeros(A.shape[0])
        Ke = np.zeros(A.shape[0])
        Pe[0] = np.dot(X[0,:],np.dot(K,X[0,:]))
        Ke[0] = np.dot(V[0,:],np.dot(M,V[0,:]))

    dt = t[1] - t[0]
    niterprint = int(len(t)/10)
    time_start = time.time()
    for it in range(len(t)-1):

        if (verbose == True) and (it % niterprint)==0:
            time_avg = (time.time() - time_start)/(it+1E-40)
            time_left = time_avg * (len(t) - it)
            print("ITERATIONS:",it,'of',str(len(t)) + str(","),"TIME AVG.:",str(time_avg) + str(","),"TIME LEFT.:",time_left)

        dumx = np.array(X[it,:])
        dumv = np.array(V[it,:])
        duma = np.array(A[it,:])

        uc = dumx + dt * dumv + 0.5 * (1 - beta2) * dt**2 * duma
        vc = dumv + (1 - beta1) * dt * duma

        Aidx = M + beta1 * dt * C + 0.5 * beta2 * dt**2 * K

        F = force_function(t[it])

        ac = np.linalg.solve(-Aidx,-F + np.dot(C,vc) + np.dot(K,uc))

        A[it+1,:] = ac
        V[it+1,:] = vc + beta1 * dt * ac
        X[it+1,:] = uc + 0.5 * beta2 * dt**2 * ac

        if calc_energy == True:
            Pe[it+1] = np.dot(X[it+1,:],np.dot(K,X[it+1,:]))
            Ke[it+1] = np.dot(V[it+1,:],np.dot(M,V[it+1,:]))

    dict_out = {"X": X,"V": V,"A": A}
    if calc_energy == True:
        dict_out["potential_energy"] = Pe
        dict_out["kinetic_energy"]   = Ke
    return dict_out

def state_space_sdof(t,x,m,c,k,force_function):
    """

    """
    import numpy as np
    import time    
    E = np.array([[0,1],[-k/m,-c/m]])
    Q = np.array([0,force_function(t,x)/m])
    return np.dot(E,x) + Q

# if __name__ == "__main__":
#
#     import scipy.linalg as scilin
#     import numpy as np
#     import scipy.optimize as sciopt
#     import matplotlib.pyplot as plt
#     import scipy.integrate as sintegrate
#
#     m1 = 10
#     m2 = 20
#
#     k1 = 10000
#     k2 = 20000
#
#     c1 = 30
#     c2 = 30
#
#     parm_M = np.diag([m1,m2])
#
#     parm_K = np.array([[k1+k2,-k2],[-k2,k2]])
#
#     parm_C = np.array([[c1+c2,-c2],[-c2,c2]])
#
#     parm_F0 = np.array([1,0])
#
#     OMEGA = 10
#
#     t = np.linspace(0,10,100000)
#
#     force_function = lambda t: parm_F0 * np.cos(OMEGA * t)
#
#     beta1 = 0.5
#     beta2 = 0.5
#
#     dict_results = newmark_gn22(t,parm_M,parm_C,parm_K,np.zeros(2),np.zeros(2),force_function,beta1,beta2,calc_energy=True)
#
#     plt.close("all")
#     plt.figure(1)
#     plt.plot(t,y[0,:],'r-')
#     plt.plot(t,X[:,0],'b--')
#
#     plt.figure(2)
#     plt.plot(t,y[1,:])
#     plt.plot(t,X[:,1],'b--')
#
#     plt.figure(3)
#     plt.plot(t,y[2,:],'r-')
#     plt.plot(t,V[:,0],'b--')
#
#     plt.figure(4)
#     plt.plot(t,y[3,:])
#     plt.plot(t,V[:,1],'b--')
#
#     plt.figure(5)
#     plt.plot(dict_results["potential_energy"])
#     plt.plot(+dict_results["kinetic_energy"])
#
#     plt.show(block = False)
