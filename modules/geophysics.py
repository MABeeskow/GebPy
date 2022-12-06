#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		geophysics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		06.12.2022

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
from scipy import signal

class geophysics:
    #
    def __init__(self, sequences):
        self.sequences = sequences
    #
    def calculatePressure(self):
        data = [0]
        g = 9.81
        for i in range(0, len(self.sequences)):
            data.append(data[i]+self.sequences[i][4]*1000*g*(self.sequences[i][3]-self.sequences[i][2]))
        #
        return data
    #
    def calculateTemperature(self):
        data = [randint(0, 25)]
        a = 30/1000
        for i in range(0, len(self.sequences)):
            data.append(round(data[0]+self.sequences[i][3]*a, 1))
        #
        return data
    #
    def calculateTTI(self):
        dataTTI = []
        for i in range(0, len(self.sequences)):
            dataTTI.append([self.sequences[i][0]])
            dataTTI[i].append([self.sequences[i][2], self.sequences[i][3]])
            dataTTI[i].append([(1-0.01*self.sequences[i][7])*1/self.sequences[i][5][0]+0.01*self.sequences[i][7]*1/self.sequences[i][5][3]])
        return dataTTI
    #
class seismics:
    #
    def __init__(self, input):
        self.input = input
    #
    def calculateImpedance(self):
        # input = sequences
        density = []
        velocity = []
        impedance = []
        for i in range(0, len(self.input)):
            density.append(self.input[i][4]*1000)
            velocity.append(self.input[i][5][0])
            impedance.append(density[i]*velocity[i])
        #
        return impedance
    #
    def calculateReflection(self):
        # input = impedance
        reflection = [0]
        for i in range(1, len(self.input)):
            reflection.append((self.input[i]-self.input[i-1])/(self.input[i]+self.input[i-1]))
        #
        return reflection
    #
    def calculatet0(self):
        # input = sequences
        velocity = []
        thickness = []
        t = []
        t0 = []
        for i in range(0, len(self.input)):
            velocity.append(self.input[i][5][0])
            thickness.append(self.input[i][1])
            t.append((2*thickness[i])/velocity[i])
        for j in range(1, len(t)+1):
            t0.append(sum(t[:j]))
        #
        return  t0
    #
    def calculate_travel_times(self):
        # input = sequences
        velocities = np.array(self.input)
        v_p = velocities[:,0]
        v_s = velocities[:,1]
        dtc = 1/v_p*10**6
        dts = 1/v_s*10**6
        #
        return dtc, dts
#
class Seismology:
    #
    def __init__(self,):
        pass
    #
    def calculate_impedance(self, velocity=None, density=None, data_all=None):
        """Returns an array that contains the impedance values of the previously generated rock units.
        **Arguments**:
            velocity: array, list of velocity values
            density: array, list of density values
        **Outputs**:
            data: array of seismic impedance values
        """
        if data_all == None and velocity != None and density != None:
            data = velocity*density
        else:
            density = []
            velocity = []
            for i in range(len(data_all)):
                for j in range(len(data_all[i])):
                    density.append(data_all[i][j][4][1][0])
                    velocity.append(data_all[i][j][4][3][0])
            data = np.array(velocity)*np.array(density)
        #
        return data
    #
    def calculate_reflection_coefficient(self, impedance):
        """Returns an array that contains the reflection coefficient values of the previously generated rock units.
        **Arguments**:
            impedance: array, list of impedance values
        **Outputs**:
            data: array of reflection coefficient values
        """
        data = [0]
        for i in range(1, len(impedance)):
            data.append((impedance[i]-impedance[i-1])/(impedance[i]+impedance[i-1]))
        #
        return np.array(data)
    #
    def calculate_transmission_coefficient(self, impedance):
        """Returns an array that contains the transmission coefficient values of the previously generated rock units.
        **Arguments**:
            impedance: array, list of impedance values
        **Outputs**:
            data: array of transmission coefficient values
        """
        data = []
        for i in range(1, len(impedance)):
            data.append((2*impedance[i-1])/(impedance[i]+impedance[i-1]))
        #
        return np.array(data)
    #
    def calculate_t0(self, thickness, velocity):
        """Returns an array that contains the two-way travel times of the previously generated rock units.
        **Arguments**:
            thickness: array, list of thickness values
            velocity: array, list of velocity values
        **Outputs**:
            data: array of two-way travel times
        """
        data = (2*thickness)/(velocity)
        #
        return np.array(data)
    #
    def create_seismic_trace(self, reflection=None, data_all=None):
        """Returns an array that contains the seismic trace based on the input data.
        **Arguments**:
            reflection: array, list of reflection coefficient values
        **Outputs**:
            data: array of seismic trace values
        """
        if data_all == None and reflection != None:
            points = len(reflection)
            width = points/100
        else:
            density = []
            velocity = []
            for i in range(len(data_all)):
                for j in range(len(data_all[i])):
                    density.append(data_all[i][j][4][1][0])
                    velocity.append(data_all[i][j][4][3][0])
            impedance = np.array(velocity)*np.array(density)
            reflection = [0]
            for i in range(1, len(impedance)):
                reflection.append((impedance[i]-impedance[i-1])/(impedance[i]+impedance[i-1]))
            points = len(reflection)
            width = points/100
        wavelet = signal.ricker(points, width)
        data = signal.convolve(reflection, wavelet, mode="same")
        #
        return data
#
class Mixing:
    #
    def __init__(self,):
        pass
    #
    def mean_arithmetic(self, weight, value):
        """Returns the arithmetic mean of a property considering a mixture
        **Arguments**:
            weight: array of weights
            value: array of property values
        **Outputs**:
            result: bulk property with respect to the arithmetic mean
        """
        result = np.sum(weight*value)
        return result
    #
    def mean_harmonic(self, weight, value):
        """Returns the harmonic mean of a property considering a mixture
        **Arguments**:
            weight: array of weights
            value: array of property values
        **Outputs**:
            result: bulk property with respect to the harmonic mean
        """
        result = (np.sum(weight/value))**(-1)
        return result
    #
    def mean_geometric(self, weight, value):
        """Returns the geometric mean of a property considering a mixture
        **Arguments**:
            weight: array of weights
            value: array of property values
        **Outputs**:
            result: bulk property with respect to the geometric mean
        """
        result = np.prod(value**weight)
        return result
    #
    def mean_squarerootmean(self, weight, value):
        """Returns the square-root mean of a property considering a mixture
        **Arguments**:
            weight: array of weights
            value: array of property values
        **Outputs**:
            result: bulk property with respect to the square-root mean
        """
        result = (np.sum(weight*np.sqrt(value)))**2
        return result
#
class Elasticity:
    #
    def __init__(self,):
        pass
    #
    def calc_voigt_bound(self, f, m):
        M = 0
        for i in range(len(f)):
            M += f[i]*m[i]
        return M
    #
    def calc_reuss_bound(self, f, m):
        Minv = 0
        for i in range(len(f)):
            Minv += f[i]*(1/m[i])
        M = 1/Minv
        return M
    #
    def calc_vrh(self, mv, mr):
        M = (mv+mr)/2
        return M
    #
    def calc_harmonic_mean(self, f, m):
        n = np.sum(f)
        a = 0
        for i in range(len(f)):
            a += f[i]/m[i]
        M = n/a
        return M
    #
    def calc_geometric_mean(self, f, m):
        a = 1
        for i in range(len(f)):
            if m[i] > 0:
                a *= m[i]**f[i]
        M = a**(1/np.sum(f))
        return M
    #
    def calc_arithmetic_mean(self, f, m):
        a = 0
        for i in range(len(f)):
            a += f[i]*m[i]
        M = a/(np.sum(f))
        return M
#
class BoreholeGeophysics:
    #
    def __init__(self,):
        pass
    #
    def calculate_pe(self, x_list, elements_list):
        contributions = [x_list[i]*(elements_list[i][1]/10)**3.6 for i in range(len(x_list))]
        pe = 0.9115*np.sum(contributions)
        return pe
    #
class WellLog:
    #
    def __init__(self, amounts=None, elements=None, rho_b=None):
        self.amounts = amounts
        self.elements = elements
        self.rho_b = rho_b
    #
    def calculate_gr(self):
        gr = 0
        for i in range(len(self.elements)):
            if self.elements[i][0] == "K":
                gr += Conversions(amount=self.amounts[i][2]).convert_to_percent()*16
            elif self.elements[i][0] == "Th":
                gr += Conversions(amount=self.amounts[i][2]).convert_to_ppm()*4
            elif self.elements[i][0] == "U":
                gr += Conversions(amount=self.amounts[i][2]).convert_to_ppm()*8
            else:
                gr += 0
        return gr
    #
    def calculate_pe(self):
        pe = 0
        for i in range(len(self.elements)):
            pe += self.amounts[i][2]*(self.elements[i][1]/10)**3.6
        return pe
    #
    def calculate_electron_density(self):
        rho_e = 0
        for i in range(len(self.elements)):
            rho_e += self.amounts[i][2]*(2*self.elements[i][1]/self.elements[i][2])*self.rho_b
        return rho_e

class Conversions:
    #
    def __init__(self, amount):
        self.amount = amount
    #
    def convert_to_percent(self):
        result = self.amount*10**2
        return result
    #
    def convert_to_ppm(self):
        result = self.amount*10**6
        return result