#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		2d_geology.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.02.2021

#-----------------------------------------------

## MODULES
import numpy as np
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
