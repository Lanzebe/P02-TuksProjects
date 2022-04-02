import numpy as np
import matplotlib.pyplot as plt
from numpy.core.machar import MachAr

class Beam:
    def __init__(self, Length):
        self.lenght_of_beam = Length
        self.shear_forces = 0
        self.bending_moments = 0
        self.x_s = np.linspace(0, Length, 1000)

    def add_distributed_loading(self, q, position_start, position_stop):
        '''
        [q]-> loading per meter and downward loading is positive.
        Position_start is where the loading starts from the left.
        Position_stop is where the loading starts from the left.
        '''
        self.shear_forces = self.shear_forces +((-q*(self.x_s-position_start))*(np.where((self.x_s>position_start) & (self.x_s<=position_stop),1,0)))#+((+q*(self.x_s-position_stop))*(self.on_or_off(position_stop))))
        self.bending_moments = self.bending_moments + ((((-q/2)*(self.x_s-position_start)**2)*(np.where((self.x_s>position_start) & (self.x_s<=position_stop),1,0))))#-(((-q/2)*(self.x_s-position_stop)**2)*(self.on_or_off(position_stop))))
    
    def add_point_load(self, Load, X_position):
        '''
        Add acording to the axis system. A positive point load is upward.
        '''
        if X_position==self.lenght_of_beam:
            self.shear_forces = self.shear_forces + (Load*(np.where((self.x_s>=X_position),1,0)))
        else:
            self.shear_forces = self.shear_forces + (Load*(np.where((self.x_s>X_position),1,0)))#-(Load*(self.on_or_off(self.lenght_of_beam)))
        self.bending_moments = self.bending_moments + ((Load*(self.x_s-X_position))*(np.where((self.x_s>=X_position) & (self.x_s<=self.lenght_of_beam),1,0)))
        #print(self.shear_forces[99])

    def add_moment(self, Moment, X_position):
        '''
        Add acording to the axis system. A positive Moment is clockwise on the whole beam.
        '''
        if X_position==self.lenght_of_beam:
            self.bending_moments = self.bending_moments + (-Moment*(np.where((self.x_s>=X_position),1,0)))
        else:
            self.bending_moments = self.bending_moments + (-Moment*(np.where((self.x_s>X_position),1,0)))
    
    def on_or_off(self, a, b=0):
        '''
        Function working with discontinuity in a loading function.
        x-> set of x axis values.
        a-> point of discontinuity(values of x greather than of equal to a will be "on")
        '''
        if b==0:
            b = self.lenght_of_beam
        if a==self.lenght_of_beam:
            return np.where(self.x_s>=a,1,0)
        else:
            return np.where((self.x_s>a) & (self.x_s<=b),1,0)
        



    def plot_beam_data(self):
        plt.figure(1)
        
        plt.subplot(211)
        plt.plot(self.x_s,self.shear_forces)
        plt.title('Shear forces')
        plt.grid()
        plt.subplot(212)
        plt.plot(self.x_s,self.bending_moments)
        plt.title('Bending moment diagram')
        plt.grid()
        plt.show()
        #add grids to graph and make axis more visable



'''
example form ST2
Ma = 7820/9
Mc = (-7820/9)+680
Beam1 = Beam(3)
Beam1.add_point_load(340,0)
Beam1.add_moment(Ma,0)
Beam1.add_moment(Mc,3)
Beam1.add_distributed_loading(170,1,3)


Beam1.plot_beam_data()
'''


'''
Beam_CT5 = Beam(10)
Beam_CT5.add_distributed_loading(500,0,10)
Beam_CT5.add_point_load(5000,0)
Beam_CT5.add_moment(25000,0)
Beam_CT5.plot_beam_data()
'''
