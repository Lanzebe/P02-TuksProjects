import cmath
from math import pi
import matplotlib.pyplot as plt
import numpy as np

class FourierSeries(object):
    def __init__(self,t_0,t_1,N_max_series = 5, Integration_n = 100):
        self.C = []
        self.t_0 = t_0
        self.t_1 = t_1
        self.L_t = t_1-t_0
        self.N_max = N_max_series
        self.num_of_integrations = Integration_n

    def numerical_int(self,fun,idx):
        dx = (self.t_1-self.t_0)/self.num_of_integrations
        sum = 0
        for i in range(self.num_of_integrations):
            t = self.t_0 + i*dx
            #print(t)
            sum = sum + 0.5*(self.e_pow_fun(fun,idx,t) + self.e_pow_fun(fun,idx,t+dx))*dx
        return sum

    def e_pow_fun(self,fun,i,t):
        return cmath.exp(-2*pi*i*1j*((t)/self.L_t))*fun(t)

    def calc_cs(self,fun_to_approx):
        for i in range(self.N_max):
            if i==0:
                self.C.append((1/self.L_t)*self.numerical_int(fun_to_approx,i))
            else:
                self.C.append((2/self.L_t)*self.numerical_int(fun_to_approx,i))
            

    def output(self,t):
        ans = 0
        for i in range(self.N_max):
            ans = ans + self.C[i]*cmath.exp(2*pi*i*1j*(t/self.L_t))
        return ans

    def outputDerivative(self,t):
        ans = 0
        for i in range(self.N_max):
            ans = ans + self.C[i]*(2*pi*i*1j*(1/self.L_t))*cmath.exp(2*pi*i*1j*(t/self.L_t))
        return ans


