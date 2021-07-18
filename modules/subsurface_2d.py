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
    def __init__(self, x_values):
        self.x_values = x_values
    #
    def create_units(self, f_surface):
        n_units = len(self.x_values)
        self.y_values = f_surface(self.x_values)
        print("x:", self.x_values)
        print("y:", self.y_values)
        f_surface1 = f_surface.derivative(nu=1)
        derivatives = f_surface1(self.x_values)
        print("Derivatives", derivatives)
        a_pos = max(derivatives)
        a_neg = min(derivatives)
        if abs(a_pos) >= abs(a_neg):
            print("Steepest point:", a_pos)
        else:
            print("Steepest point:", a_neg)