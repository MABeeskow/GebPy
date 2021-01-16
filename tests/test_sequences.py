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
data = sequences.SedimentaryBasin()
data_soil = data.create_soil()
for i in range(len(data_soil)):
    print(data_soil[i])

print("")
# Test sand generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sand = data.create_sand()
for i in range(len(data_sand)):
    print(data_sand[i])

print("")
# Test sandstone generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone()
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
# Test shale generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_shale = data.create_shale()
for i in range(len(data_shale)):
    print(data_shale[i])

print("")
# Test sedimentary basin generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sedbasin = data.create_sedimentary_basin()
for i in range(len(data_sedbasin)):
    print(data_sedbasin[i])