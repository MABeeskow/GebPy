#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_minerals.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		04.04.2021

# -----------------------------------------------

## MODULES
import sys
import numpy as np
import random as rd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from modules import sequences, geophysics, minerals

## TESTING ORGANICS
print("TEST - ORGANICS")
# Test Carbohydrates
data = minerals.Organics.carbohydrates("")
print("Carbohydrates:\n", data)
# Test Lignin
data = minerals.Organics.lignin("")
print("Lignin:\n", data)
# Test Lipid
data = minerals.Organics.lipid("")
print("Lipid:\n", data)
# Test Organic matter
data = minerals.Organics.organic_matter("")
print("Organic matter:\n", data)

print("")
## TESTING PHYLLOSILICATES
print("TEST - PHYLLOSILICATES")
# Test Illite
data = minerals.phyllosilicates.illite("")
print("Illite:\n", data)