#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_sequences.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		14.01.2021

# -----------------------------------------------

## MODULES
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from modules import sequences

## TESTING
# Test soil generation within SedimentaryBasin class
actual_thickness = 0
data = sequences.SedimentaryBasin(actualThickness=actual_thickness, parts=10)
data_soil = data.create_soil()
for i in range(len(data_soil)):
    print(data_soil[i])

# Test sand generation within SedimentaryBasin class
actual_thickness = 0
data = sequences.SedimentaryBasin(actualThickness=actual_thickness, parts=10)
data_sand = data.create_sand()
for i in range(len(data_sand)):
    print(data_sand[i])