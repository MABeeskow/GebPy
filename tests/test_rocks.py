#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		13.03.2021

# -----------------------------------------------

## MODULES
import sys
import numpy as np
import random as rd
import  matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from modules import sequences, geophysics, siliciclastics

## TESTING
# Test soil generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_soil = data.create_soil(grainsize=True, thickness=10)
for i in range(len(data_soil)):
    print(data_soil[i])

print("")
# Test sand generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sand = data.create_sand(grainsize=True, thickness=10)
for i in range(len(data_sand)):
    print(data_sand[i])

print("")
# Test sandstone generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(thickness=10)
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(fluid="oil", thickness=10)
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(fluid="gas", thickness=10)
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(fluid="air", thickness=10)
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(fluid="air", thickness=10, phi=0.15, Qz_rich=True)
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
# Test limestone generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_limestone = data.create_limestone(thickness=10)
for i in range(len(data_limestone)):
    print(data_limestone[i])

print("")
data = sequences.SedimentaryBasin()
data_limestone = data.create_limestone(fluid="oil", thickness=10)
for i in range(len(data_limestone)):
    print(data_limestone[i])

print("")
data = sequences.SedimentaryBasin()
data_limestone = data.create_limestone(fluid="gas", thickness=10)
for i in range(len(data_limestone)):
    print(data_limestone[i])

print("")
# Test shale generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_shale = data.create_shale(thickness=10)
for i in range(len(data_shale)):
    print(data_shale[i])

#print("")
# Test shale generation within SedimentaryBasin class
#data = sequences.SedimentaryBasin()
#data_shale = data.create_shale_complex(thickness=10)
#for i in range(len(data_shale)):
#    print(data_shale[i])

print("")
# Test rock salt generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_rocksalt = data.create_rocksalt(thickness=10)
for i in range(len(data_rocksalt)):
    print(data_rocksalt[i])

print("")
print(":: IGNEOUS ROCKS ::")
# Test granite generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_granite = data.create_granite(thickness=10)
for i in range(len(data_granite)):
    print(data_granite[i])

print("")
# Test basalt generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_basalt = data.create_basalt(thickness=10)
for i in range(len(data_basalt)):
    print(data_basalt[i])

print("")
print(":: IGNEOUS ROCKS (PLUTONITE) ::")
# Test granite generation within Plutonite class
data = sequences.Plutonite()
data_rock = data.create_granite(thickness=10)
for i in range(len(data_rock)):
    print(data_rock[i])
print("")
# Test gabbro generation within Plutonite class
data = sequences.Plutonite()
data_rock = data.create_gabbro(thickness=10)
for i in range(len(data_rock)):
    print(data_rock[i])
print("")
# Test felsic rock generation within Plutonite class
data = sequences.Plutonite()
data_rock = data.create_felsic(thickness=10)
for i in range(len(data_rock)):
    print(data_rock[i])
print("")
# Test intermediate rock generation within Plutonite class
data = sequences.Plutonite()
data_rock = data.create_intermediate(thickness=10)
for i in range(len(data_rock)):
    print(data_rock[i])

print("")
print(":: NAGRA Benken ::")
# Test rock generation based on the NAGRA Benken dataset
for i in range(10):
    data = siliciclastics.NAGRA()
    data_benken = data.create_benken_rocks("")
    print(data_benken)

print("")
# Test rock generation based on the NAGRA Benken dataset within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_benken = data.create_benken(thickness=10)
for i in range(len(data_benken)):
    print(data_benken[i])