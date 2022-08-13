#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		series.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		13.08.2020

#-----------------------------------------------

## MODULES
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
    def __init__(self, thickness=1000, resolution=1, composition=None):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
    #
    def create_zechstein_z1(self):  # Z1 - Werra Series
        #
        data_kupferschiefer = Ores(
            fluid="water", actualThickness=0, porosity=rd.uniform(0.0, 0.05), data_type=True).create_kupferschiefer()
        data_limestone = limestone(
            fluid="water", actualThickness=0).create_simple_limestone(
            dict=True, porosity=rd.uniform(0.1, 0.4))
        data_anhydrite = Evaporites(
            fluid="water", actualThickness=0).create_simple_anhydrite(
            dict=True, porosity=rd.uniform(0.05, 0.2))
        print("Kupferschiefer:", data_kupferschiefer)
        print("Limestone:", data_limestone)
        print("Anhydrite:", data_anhydrite)