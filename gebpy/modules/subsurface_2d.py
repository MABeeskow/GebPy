#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		subsurface_2d.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.07.2021

#-----------------------------------------------

## MODULES
import numpy as np
import scipy.interpolate as interp
from numpy import round
import random as rd
from itertools import chain
from random import randint
from modules.carbonates import limestone, dolomite
from modules.siliciclastics import sandstone, shale, ore, Soil
from modules.igneous import plutonic, volcanic, Plutonic, Volcanic
from modules.evaporites import evaporites, Evaporites
from modules import minerals, sequences
from modules.elements import elements
from modules import fluids

class structural_geology:
    #
    def __init__(self):
        pass
    #
    #def create_tilted_structures(self, y, dip=30):
    #    x = y*np.tan(dip*np.pi/180)

    #
    def generate_2d_horizontal(self, max_thickness=500, n=10, x=100):
        data_boreholes = []
        data = sequences.SedimentaryBasin()
        data_sedbasin = data.create_sedimentary_basin(maximum_thickness=max_thickness)
        for i in range(n):
            data_boreholes.append([data_sedbasin, x*i])
        #
        return data_boreholes
    #
    #def generate_2d_tilted(self, dip=20, max_thickness=500, n=10, x=50):
    #    data_boreholes = []
    #    data = sequences.SedimentaryBasin()
    #    data_sedbasin = data.create_sedimentary_basin(maximum_thickness=max_thickness)

class Surface:
    #
    def __init__(self, coordinates):
        self.coordinates = coordinates
    #
    def create_linear_surface(self):
        # y = a*x + b
        x1 = self.coordinates[0][0]
        y1 = self.coordinates[0][1]
        x2 = self.coordinates[1][0]
        y2 = self.coordinates[1][1]
        #
        b = (x1*y2 - x2*y1)/(x1 - x2)
        a = (y1 - b)/x1
        #
        return a, b
    #
    def create_quadratic_surface(self):
        # y = a*x**2 + b*x + c
        x1 = self.coordinates[0][0]
        y1 = self.coordinates[0][1]
        x2 = self.coordinates[1][0]
        y2 = self.coordinates[1][1]
        x3 = self.coordinates[2][0]
        y3 = self.coordinates[2][1]
        #
        a = ((x2 - x3)*(y1 - y3) - (x1 - x3)*(y2 - y3))/((x2 - x3)*(x1**2 - x3**2) - (x1 - x3)*(x2**2 - x3**2))
        b = (y1 - y3 - a*(x1**2 - x3**2))/(x1 - x3)
        c = y3 - a*x3**2 - b*x3
        #
        return a, b, c
    #
    def create_interpolated_surface(self):
        self.coordinates = np.array(self.coordinates)
        x = self.coordinates[:, 0]
        y = self.coordinates[:, 1]
        #
        f = interp.CubicSpline(x, y)
        #
        return f

class Units:
    #
    def __init__(self, x_limits, y_limits, number):
        self.x_limits = x_limits
        self.y_limits = y_limits
        self.number = number
    #
    def create_outcrops(self, extrema):
        condition = False
        while condition == False:
            x_outcrops = np.sort(rd.sample(list(range(int(self.x_limits[0]*(1-0.05)),
                                                      int(extrema[-1][0]-5))) + list(range(int(extrema[-1][0]+5),
                                                      int(self.x_limits[1]*(1-0.05)))), self.number))
            x_diff = np.diff(x_outcrops)
            check = any(item in np.arange(0, 5) for item in x_diff)
            if check == False:
                condition = True
            else:
                continue
        #
        return x_outcrops
    #
    def create_units(self, f_surface, kind=None):
        x_units = []
        y_units = []
        f1 = interp.CubicSpline.derivative(f_surface)
        f2 = interp.CubicSpline.derivative(f1)
        x_extreme = np.around(f1.solve(), 4)
        y_extreme = np.around(f_surface(x_extreme), 4)
        extrema_type = np.around(f2(x_extreme), 4)
        data_extrema = []
        for i in range(len(x_extreme)):
            if extrema_type[i] < 0:
                data_extrema.append([round(x_extreme[i], 2), round(y_extreme[i], 2), "Maximum"])
            elif extrema_type[i] > 0:
                data_extrema.append([round(x_extreme[i], 2), round(y_extreme[i], 2), "Minimum"])
            elif extrema_type[i] == 0:
                data_extrema.append([round(x_extreme[i], 2), round(y_extreme[i], 2), "Saddle"])
        if len(x_extreme) == 0:
            data_extrema.append([round(max(self.x_limits), 2), round(max(self.y_limits), 2), "Maximum"])
            data_extrema.append([round(min(self.x_limits), 2), round(min(self.y_limits), 2), "Minimum"])
        data_extrema = np.array(data_extrema, dtype="object")
        data_extrema = data_extrema[data_extrema[:, 1].argsort()]
        print("Global Minimum:", data_extrema[0])
        print("Global Maximum:", data_extrema[-1])

        self.x_outcrops = self.create_outcrops(extrema=data_extrema)
        self.y_outcrops = f_surface(self.x_outcrops)
        n_units = len(self.x_outcrops)

        if kind == None:
            x_delta = rd.randint(-25, 25)
            self.x_outcrops_bottom = self.x_outcrops + x_delta
            y_min = min(self.y_limits)
            for i in range(n_units):
                if self.x_outcrops is not self.x_outcrops_bottom:
                    x1 = self.x_outcrops[i]
                    x2 = self.x_outcrops_bottom[i]
                    y1 = self.y_outcrops[i]
                    y2 = y_min
                    b = (x1*y2 - x2*y1)/(x1 - x2)
                    if x1 != 0:
                        a = (y1 - b)/x1
                    else:
                        a = (y2 - b)/x2
                    x = np.linspace(self.x_outcrops[i], self.x_outcrops_bottom[i], 25+1)
                    upper_limit = np.argwhere(x > self.x_limits[-1])
                    if len(upper_limit) > 0:
                        x = x[:upper_limit[0][0]]
                    x_units.append(x)
                    y_units.append(a*x + b)
                else:
                    x1 = self.x_outcrops[i]
                    x2 = self.x_outcrops_bottom[i]
                    y1 = self.y_outcrops[i]
                    y2 = y_min
                    x_units.append(np.array([x1, x2]))
                    y_units.append(np.array([y1, y2]))
        elif kind == "horizontal":
            for i in range(len(self.x_outcrops)):
                if self.x_outcrops[i] < data_extrema[-1][0]:
                    x_units.append([self.x_outcrops[i], f_surface.solve(y=self.y_outcrops[i])[2]])
                    y_units.append([self.y_outcrops[i], self.y_outcrops[i]])
                elif self.x_outcrops[i] > data_extrema[-1][0]:
                    x_units.append([self.x_outcrops[i], f_surface.solve(y=self.y_outcrops[i])[1]])
                    y_units.append([self.y_outcrops[i], self.y_outcrops[i]])
        #
        return x_units, y_units, self.x_outcrops