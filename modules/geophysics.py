#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		geophysics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		16.02.2020

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
from modules import minerals
from modules.elements import elements
from scipy import signal
import pywt
import matplotlib.pyplot as plt

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
    def calculateRickerWavelet(self):
        # input = reflection
        width = 2.0
        wavelet = signal.ricker(self.input, width)
        #
        return wavelet
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
    def calculateTrace(self):
        # input = reflection
        widths = np.arange(1, len(self.input)+1)
        cwt, freq = pywt.cwt(self.input, widths, "mexh", method="conv")
        data = [cwt, freq]
        #
        #plt.imshow(cwt, extent=[-1, 1, 1, len(self.input)], cmap='viridis', aspect='auto', vmax=abs(cwt).max(), vmin=-abs(cwt).max())
        #plt.show()
        #
        return data
    #
    def seismic_wiggle(self):
        dt=0.004
        ranges = None
        normalize=False
        scale=1.
        color='k'
        input = np.asarray([self.input])
        npts, ntraces = input.shape
        #npts = len(self.input)
        #ntraces = 1
        if ntraces < 1:
            raise IndexError("Nothing to plot")
        if npts < 1:
            raise IndexError("Nothing to plot")
        t = np.linspace(0, dt*npts, npts)
        amp = 1.  # normalization factor
        gmin = 0.  # global minimum
        toffset = 0.  # offset in time to make 0 centered
        if normalize:
            gmax = input.max()
            gmin = input.min()
            amp = (gmax-gmin)
            toffset = 0.5
        plt.ylim(max(t), 0)
        if ranges is None:
            ranges = (0, ntraces)
        x0, x1 = ranges
        # horizontal increment
        dx = float((x1-x0)/ntraces)
        plt.xlim(x0, x1)
        for i, trace in enumerate(input.transpose()):
            tr = (((trace-gmin)/amp)-toffset)*scale*dx
            x = x0+i*dx  # x positon for this trace
            plt.plot(x+tr, t, 'k')
            plt.fill_betweenx(t, x+tr, x, tr > 0, color=color)
            plt.show()
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
        pe = np.sum(contributions)
        return pe