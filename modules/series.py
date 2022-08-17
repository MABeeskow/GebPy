#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		series.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		13.08.2020

#-----------------------------------------------

## MODULES
import math
import random as rd
import numpy as np
import scipy
from ore import Ores
from carbonates import limestone
from evaporites import Evaporites

#######################
## SERIES GENERATION ##
#######################
#
## ROTLIEGEND
class Zechstein:
    #
    def __init__(self, thickness=1000, resolution=10, composition=None):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
        #
        self.actual_thickness = 0
    #
    def create_zechstein_z1(self, thickness_z1=100):  # Z1 - Werra Series
        fraction_kupferschiefer_pre = round(rd.uniform(1, 10), 4)
        fraction_kupferschiefer = round(round(fraction_kupferschiefer_pre*2)/2/100, 4)
        fraction_limestone_pre = round(rd.uniform(15, 30), 4)
        fraction_limestone = round(round(fraction_limestone_pre*2)/2/100, 4)
        fraction_anhydrite = round(1 - fraction_kupferschiefer - fraction_limestone, 4)
        #
        thickness_kupferschiefer = round(thickness_z1*fraction_kupferschiefer, 4)
        thickness_limestone = round(thickness_z1*fraction_limestone, 4)
        thickness_anhydrite = round(thickness_z1*fraction_anhydrite, 4)
        #
        self.actual_thickness += thickness_kupferschiefer/self.resolution
        ## Create Kupferschiefer Unit
        container_kupferschiefer = {}
        for i in np.linspace(self.actual_thickness, self.actual_thickness+thickness_kupferschiefer, self.resolution,
                             endpoint=False):
            depth = round(thickness_z1 - i, 2)
            container_kupferschiefer[depth] = Ores(fluid="water", actualThickness=0, porosity=rd.uniform(0.0, 0.05),
                                               data_type=True).create_kupferschiefer()
        self.actual_thickness = thickness_kupferschiefer + thickness_limestone/self.resolution
        ## Create Limestone Unit
        container_limestone = {}
        for i in np.linspace(self.actual_thickness, self.actual_thickness+thickness_limestone, self.resolution,
                             endpoint=False):
            depth = round(thickness_z1 - i, 2)
            container_limestone[depth] = limestone(fluid="water", actualThickness=0).create_simple_limestone(
                dict=True, porosity=rd.uniform(0.1, 0.4))
        self.actual_thickness = thickness_kupferschiefer + thickness_limestone + thickness_anhydrite/self.resolution
        ## Create Anhydrite Unit
        container_anhydrite = {}
        for i in np.linspace(self.actual_thickness, self.actual_thickness+thickness_anhydrite, self.resolution,
                             endpoint=False):
            depth = round(thickness_z1 - i, 2)
            container_anhydrite[depth] = Evaporites(fluid="water", actualThickness=0).create_simple_anhydrite(
                dict=True, porosity=rd.uniform(0.05, 0.2))
        self.actual_thickness = thickness_kupferschiefer + thickness_limestone + thickness_anhydrite
        #
        ## TEST
        # for key, value in container_kupferschiefer.items():
        #     print(key, value)
        # for key, value in container_limestone.items():
        #     print(key, value)
        # for key, value in container_anhydrite.items():
        #     print(key, value)
        #
        return container_kupferschiefer, container_limestone, container_anhydrite