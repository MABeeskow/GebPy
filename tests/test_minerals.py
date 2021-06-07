#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_minerals.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		04.06.2021

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
        print("TEST - OXIDES (from oxides.py)")
        # Test Quartz (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_quartz()
        print("Quartz (pure):\n", data, "\n")
        data = oxides.Oxides(impurity="random").create_quartz()
        print("Quartz (incl. trace elements):\n", data, "\n")
        # Test Uraninite (incl. trace elements)
        data = oxides.Oxides(impurity="random").create_uraninite()
        print("Uraninite (incl. trace elements):\n", data)
        # Test Magnetite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_magnetite()
        print("Magnetite (incl. trace elements):\n", data)
        # Test Hematite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_hematite()
        print("Hematite (incl. trace elements):\n", data)
        # Test Corundum (incl. trace elements)
        data = oxides.Oxides(impurity="random").create_corundum()
        print("Corundum (incl. trace elements):\n", data)
        # Test Wustite (incl. trace elements)
        data = oxides.Oxides(impurity="random").create_wustite()
        print("Wustite (incl. trace elements):\n", data)
        # Test Chromite (incl. trace elements)
        data = oxides.Oxides(impurity="random").create_chromite()
        print("Chromite (incl. trace elements):\n", data)
        # Test Spinel (incl. trace elements)
        data = oxides.Oxides(impurity="random").create_spinel()
        print("Spinel (incl. trace elements):\n", data)



        print("")
        data = oxides.Oxides(impurity="random").create_quartz()
        print("Quartz (incl. trace elements):\n", data)
# RUN
TESTING_MINERALS()