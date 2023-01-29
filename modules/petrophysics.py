#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		petrophysics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		29.01.2023

#-----------------------------------------------

## MODULES
import numpy as np
import random as rd

## METHODS
class SolidProperties:
    #
    def __init__(self, value_amount, value_property):
        self.value_amount = value_amount
        self.value_property = value_property
    #
    def calculate_solid_density(self):
        pass
