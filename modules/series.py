#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		series.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		09.08.2020

#-----------------------------------------------

## MODULES
import numpy as np
import scipy

#######################
## SERIES GENERATION ##
#######################
#
## ROTLIEGEND
class Rotliegend:
    #
    def __init__(self, thickness, resolution, composition):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
    #
    def create_volcanites(self):
