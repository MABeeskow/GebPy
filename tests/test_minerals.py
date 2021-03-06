#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_minerals.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		21.06.2021

# -----------------------------------------------

## MODULES
import sys
import numpy as np
import random as rd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from modules import sequences, geophysics, minerals, oxides, sulfides, sulfates

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
        print("Quartz (pure):\n", data)
        # Test Uraninite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_uraninite()
        print("Uraninite (incl. trace elements):\n", data)
        # Test Magnetite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_magnetite()
        print("Magnetite (incl. trace elements):\n", data)
        # Test Hematite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_hematite()
        print("Hematite (incl. trace elements):\n", data)
        # Test Corundum (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_corundum()
        print("Corundum (incl. trace elements):\n", data)
        # Test Wustite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_wustite()
        print("Wustite (incl. trace elements):\n", data)
        # Test Chromite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_chromite()
        print("Chromite (incl. trace elements):\n", data)
        # Test Spinel (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_spinel()
        print("Spinel (incl. trace elements):\n", data)
        # Test Boehmite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_boehmite()
        print("Boehmite (incl. trace elements):\n", data)
        # Test Diaspore (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_diaspore()
        print("Diaspore (incl. trace elements):\n", data)
        # Test Gibbsite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_gibbsite()
        print("Gibbsite (incl. trace elements):\n", data)
        # Test Cuprite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_cuprite()
        print("Cuprite (incl. trace elements):\n", data)
        # Test Goethite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_goethite()
        print("Goethite (incl. trace elements):\n", data)
        # Test Ilmenite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_ilmenite()
        print("Ilmenite (incl. trace elements):\n", data)
        # Test Rutile (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_rutile()
        print("Rutile (incl. trace elements):\n", data)
        # Test Brookite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_brookite()
        print("Brookite (incl. trace elements):\n", data)
        # Test Anatase (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_anatase()
        print("Anatase (incl. trace elements):\n", data)
        # Test Manganite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_manganite()
        print("Manganite (incl. trace elements):\n", data)
        # Test Groutite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_groutite()
        print("Groutite (incl. trace elements):\n", data)
        # Test Pyrophanite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_pyrophanite()
        print("Pyrophanite (incl. trace elements):\n", data)
        # Test Geikielite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_geikielite()
        print("Geikielite (incl. trace elements):\n", data)
        # Test Eskolaite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_eskolaite()
        print("Eskolaite (incl. trace elements):\n", data)
        # Test Karelianite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_karelianite()
        print("Karelianite (incl. trace elements):\n", data)
        # Test Claudetite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_claudetite()
        print("Claudetite (incl. trace elements):\n", data)
        # Test Arsenolite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_arsenolite()
        print("Arsenolite (incl. trace elements):\n", data)
        # Test Senarmontite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_senarmontite()
        print("Senarmontite (incl. trace elements):\n", data)
        # Test Valentinite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_valentinite()
        print("Valentinite (incl. trace elements):\n", data)
        # Test Bismite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_bismite()
        print("Bismite (incl. trace elements):\n", data)
        # Test Sphaerobismite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_sphaerobismite()
        print("Sphaerobismite (incl. trace elements):\n", data)
        # Test Cassiterite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_cassiterite()
        print("Cassiterite (incl. trace elements):\n", data)
        # Test Pyrolusite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_pyrolusite()
        print("Pyrolusite (incl. trace elements):\n", data)
        # Test Brucite (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_brucite()
        print("Brucite (incl. trace elements):\n", data)
        # Test Ulvöspinel (incl. trace elements)
        data = oxides.Oxides(impurity="pure").create_ulvospinel()
        print("Ulvöspinel (incl. trace elements):\n", data)

        print("")
        ## TESTING SULFIDES
        print("TEST - SULFIDES (from sulfides.py)")
        # Test Cinnabar (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_cinnabar()
        print("Cinnabar:\n", data)
        # Test Pyrite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_pyrite()
        print("Pyrite:\n", data)
        # Test Bornite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_bornite()
        print("Bornite:\n", data)
        # Test Galena (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_galena()
        print("Galena:\n", data)
        # Test Chalcopyrite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_chalcopyrite()
        print("Chalcopyrite:\n", data)
        # Test Molybdenite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_molybdenite()
        print("Molybdenite:\n", data)
        # Test Sphalerite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_sphalerite()
        print("Sphalerite:\n", data)
        # Test Stibnite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_stibnite()
        print("Stibnite:\n", data)
        # Test Arsenopyrite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_arsenopyrite()
        print("Arsenopyrite:\n", data)
        # Test Acanthite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_acanthite()
        print("Acanthite:\n", data)
        # Test Argentite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_argentite()
        print("Argentite:\n", data)
        # Test Alabandite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_alabandite()
        print("Alabandite:\n", data)
        # Test Berthierite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_berthierite()
        print("Berthierite:\n", data)
        # Test Pyrrhotite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_pyrrhotite()
        print("Pyrrhotite:\n", data)
        # Test Cobaltite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_cobaltite()
        print("Cobaltite:\n", data)
        # Test Carrollite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_carrollite()
        print("Carrollite:\n", data)
        # Test Chalcocite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_chalcocite()
        print("Chalcocite:\n", data)
        # Test Digenite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_digenite()
        print("Digenite:\n", data)
        # Test Tenanntite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_tennantite()
        print("Tenanntite:\n", data)
        # Test Tetrahedrite (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_tetrahedrite()
        print("Tetrahedrite:\n", data)
        # Test Tenanntite-Group (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_tennantite_group()
        print("Tenanntite-Group:\n", data)
        # Test Tetrahedrite-Group (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_tetrahedrite_group()
        print("Tetrahedrite-Group:\n", data)
        # Test Fahlore (incl. trace elements)
        data = sulfides.Sulfides(impurity="pure").create_fahlore()
        print("Fahlore:\n", data)

        print("")
        ## TESTING SULFATES
        print("TEST - SULFATES (from sulfates.py)")
        # Test Anhydrite (incl. trace elements)
        data = sulfates.Sulfates(impurity="pure").create_anhydrite()
        print("Anhydrite:\n", data)
        # Test Gpysum (incl. trace elements)
        data = sulfates.Sulfates(impurity="pure").create_gypsum()
        print("Gpysum:\n", data)
        # Test Scheelite (incl. trace elements)
        data = sulfates.Sulfates(impurity="pure").create_scheelite()
        print("Scheelite:\n", data)
        # Test Barite (incl. trace elements)
        data = sulfates.Sulfates(impurity="pure").create_barite()
        print("Barite:\n", data)
# RUN
TESTING_MINERALS()