#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_minerals.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.05.2021

# -----------------------------------------------

## MODULES
import sys
import numpy as np
import random as rd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from modules import sequences, geophysics, minerals, oxides

class TESTING_MINERALS:
    #
    def __init__(self):
        ## TESTING NATIVES
        print("TEST - NATIVES")
        # Test Carbohydrates
        data = minerals.natives.sulfur("")
        print("Sulfur:\n", data)

        print("")
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
        # Test Vermiculite
        data = minerals.phyllosilicates.vermiculite("")
        print("Vermiculite:\n", data)
        # Test Clays
        data = minerals.phyllosilicates.clay("")
        print("Clays (Ilt,Mnt,Kln):\n", data)

        print("")
        ## TESTING INOSILICATES
        print("TEST - INOSILICATES")
        # Test Riebeckite
        data = minerals.inosilicates.riebeckite("")
        print("Riebeckite:\n", data)
        # Test Arfvedsonite
        data = minerals.inosilicates.arfvedsonite("")
        print("Arfvedsonite:\n", data)
        # Test Sodium Amphiboles
        data = minerals.inosilicates.amphibole_na("")
        print("Sodium Amphiboles (Rbk,Arf):\n", data)
        # Test Calcium Amphiboles
        data = minerals.inosilicates.amphibole_ca("")
        print("Calcium Amphiboles (Tr,Act):\n", data)
        # Test Pyroxenes
        data = minerals.inosilicates.pyroxene("")
        print("Pyroxene (En,Aug):\n", data)
        # Test Calcium Pyroxenes
        data = minerals.inosilicates.pyroxene_ca("")
        print("Calcium Pyroxene (Aug,Di):\n", data)

        print("")
        ## TESTING FELDSPARS
        print("TEST - FELDSPARS")
        # Test Olivine
        data = minerals.feldspars.alkalifeldspar(self, keyword="None")
        print("Alkalifeldspars:\n", data)

        print("")
        ## TESTING NESOSILICATES
        print("TEST - NESOSILICATES")
        # Test Olivine
        data = minerals.nesosilicates.olivine(self, keyword="None")
        print("Olivine:\n", data)

        print("")
        ## TESTING OXIDES
        print("TEST - OXIDES")
        # Test Quartz
        data = minerals.oxides.quartz("")
        print("Quartz:\n", data)
        # Test Quartz (with Al traces)
        data = minerals.oxides.quartz(self, tr_Al=True)
        print("Quartz:\n", data)
        # Test Quartz (with Ti traces)
        data = minerals.oxides.quartz(self, tr_Ti=True)
        print("Quartz:\n", data)
        # Test Quartz (with Li traces)
        data = minerals.oxides.quartz(self, tr_Li=True)
        print("Quartz:\n", data)
        ## TESTING OXIDES
        print("TEST - OXIDES (from oxides.py)")
        # Test Quartz
        data = oxides.Quartz().create_quartz()
        print("Quartz:\n", data)
# RUN
TESTING_MINERALS()