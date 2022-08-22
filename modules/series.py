#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		series.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		18.08.2020

#-----------------------------------------------

## MODULES
import math
import random as rd
import numpy as np
import scipy
from ore import Ores
from carbonates import limestone, CarbonateRocks
from evaporites import Evaporites
from modules.siliciclastics import shale, Sandstone

#######################
## SERIES GENERATION ##
#######################
#
## ROTLIEGEND
class Zechstein:
    #
    def __init__(self, actual_thickness=0, thickness=1000, resolution=25, composition=None):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
        #
        self.actual_thickness = actual_thickness
    #
    def create_zechstein_z1(self, thickness_z1=100, top_z=0):  # Z1 - Werra Series
        fraction_kupferschiefer_pre = round(rd.uniform(2.5, 15), 4)
        fraction_kupferschiefer = round(round(fraction_kupferschiefer_pre*2)/2/100, 4)
        fraction_limestone_pre = round(rd.uniform(15, 30), 4)
        fraction_limestone = round(round(fraction_limestone_pre*2)/2/100, 4)
        fraction_anhydrite = round(1 - fraction_kupferschiefer - fraction_limestone, 4)
        #
        thickness_kupferschiefer = round(thickness_z1*fraction_kupferschiefer, 4)
        thickness_limestone = round(thickness_z1*fraction_limestone, 4)
        thickness_anhydrite = round(thickness_z1*fraction_anhydrite, 4)
        #
        actual_top = top_z
        actual_bottom = top_z + thickness_anhydrite
        #
        ## Create Anhydrite Unit
        container_anhydrite = {}
        steps_anhydrite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite:
            depth = round(i, 4)
            container_anhydrite[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.2])
        actual_top += thickness_anhydrite
        actual_bottom += thickness_limestone
        #
        ## Create Limestone Unit
        container_limestone = {}
        steps_limestone = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_limestone:
            depth = round(i, 4)
            # container_limestone[depth] = limestone(fluid="water", actualThickness=0).create_simple_limestone(
            #     dict=True, porosity=rd.uniform(0.1, 0.4))
            container_limestone[depth] = CarbonateRocks(fluid="water", actualThickness=0).create_limestone(
                number=1, porosity=[0.1, 0.4])
        actual_top += thickness_limestone
        actual_bottom += thickness_kupferschiefer
        #
        ## Create Kupferschiefer Unit
        container_kupferschiefer = {}
        steps_kupferschiefer = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_kupferschiefer:
            depth = round(i, 4)
            container_kupferschiefer[depth] = Ores(fluid="water", actualThickness=0, porosity=rd.uniform(0.0, 0.05),
                                               data_type=True).create_kupferschiefer()
        #
        ## TEST
        # for key, value in reversed(container_anhydrite.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_anhydrite, container_limestone, container_kupferschiefer
    #
    def create_zechstein_z2(self, thickness_z2=300, top_z=0):  # Z2 - Straßfurt Series
        condition = False
        while condition == False:
            counter = 0
            #
            fraction_dolomite_pre = round(rd.uniform(5, 15), 4)
            fraction_dolomite = round(round(fraction_dolomite_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_lower_pre = round(rd.uniform(20, 30), 4)
            fraction_anhydrite_lower = round(round(fraction_anhydrite_lower_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_lower_pre = round(rd.uniform(50, 70), 4)
            fraction_rocksalt_lower = round(round(fraction_rocksalt_lower_pre * 2) / 2 / 100, 4)
            fraction_potash_pre = round(rd.uniform(10, 20), 4)
            fraction_potash = round(round(fraction_potash_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_upper_pre = round(rd.uniform(1, 3), 4)
            fraction_rocksalt_upper = round(round(fraction_rocksalt_upper_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_upper = round(1 - fraction_dolomite - fraction_anhydrite_lower - fraction_rocksalt_lower
                                             - fraction_potash - fraction_rocksalt_upper, 4)
            #
            thickness_dolomite = round(thickness_z2 * fraction_dolomite, 4)
            thickness_anhydrite_lower = round(thickness_z2 * fraction_anhydrite_lower, 4)
            thickness_rocksalt_lower = round(thickness_z2 * fraction_rocksalt_lower, 4)
            thickness_potash = round(thickness_z2 * fraction_potash, 4)
            thickness_rocksalt_upper = round(thickness_z2 * fraction_rocksalt_upper, 4)
            thickness_anhydrite_upper = round(thickness_z2 * fraction_anhydrite_upper, 4)
            list_units = ["Anhydrite (upper)", "Rock Salt (upper)", "Potash", "Rock Salt (lower)", "Anhydrite (lower)",
                          "Dolomite"]
            list_thickness = [thickness_anhydrite_upper, thickness_rocksalt_upper, thickness_potash,
                              thickness_rocksalt_lower, thickness_anhydrite_lower, thickness_dolomite]
            for thickness in list_thickness:
                if thickness > 0:
                    counter += 1
                else:
                    counter += 0
            if counter == len(list_thickness) and np.sum(list_thickness) == thickness_z2:
                condition = True
        #
        # for index, unit in enumerate(list_units):
        #     print(unit, list_thickness[index])
        # print("Total Thickness:", np.sum(list_thickness))
        self.actual_thickness += thickness_dolomite/self.resolution
        actual_top = top_z
        actual_bottom = top_z + thickness_anhydrite_upper
        # print("Anhydrite (upper):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Anhydrite (upper) Unit
        container_anhydrite_upper = {}
        steps_anhydrite_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite_upper:
            depth = round(i, 4)
            container_anhydrite_upper[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.0, 0.1])
        actual_top += thickness_anhydrite_upper
        actual_bottom += thickness_rocksalt_upper
        # print("Rock Salt (upper):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Rock Salt (upper) Unit
        container_rocksalt_upper = {}
        steps_rocksalt_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rocksalt_upper:
            depth = round(i, 4)
            container_rocksalt_upper[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt_upper
        actual_bottom += thickness_potash
        # print("Potash:", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Potash Unit
        container_potash = {}
        steps_potash = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_potash:
            depth = round(i, 4)
            container_potash[depth] = Evaporites(fluid="water", actualThickness=0).create_potash(
                porosity=[0.0, 0.05])
        actual_top += thickness_potash
        actual_bottom += thickness_rocksalt_lower
        # print("Rock Salt (lower):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Rock Salt (lower) Unit
        container_rocksalt_lower = {}
        steps_rocksalt_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rocksalt_lower:
            depth = round(i, 4)
            container_rocksalt_lower[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt_lower
        actual_bottom += thickness_anhydrite_lower
        # print("Anhydrite (lower):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Anhydrite (lower) Unit
        container_anhydrite_lower = {}
        steps_anhydrite_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite_lower:
            depth = round(i, 4)
            container_anhydrite_lower[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.0, 0.1])
        actual_top += thickness_anhydrite_lower
        actual_bottom += thickness_dolomite
        # print("Dolomite:", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Dolomite Unit
        container_dolomite = {}
        steps_dolomite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_dolomite:
            depth = round(i, 4)
            container_dolomite[depth] = CarbonateRocks().create_dolomite(number=1, porosity=[0.1, 0.4])
        actual_top += thickness_dolomite
        actual_bottom += thickness_dolomite
        #
        self.actual_thickness = thickness_dolomite + thickness_anhydrite_lower + thickness_rocksalt_lower \
                                + thickness_potash + thickness_rocksalt_upper + thickness_anhydrite_upper
        #
        ## TEST
        # for key, value in reversed(container_anhydrite_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_rocksalt_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_potash.items()):
        #     print(key, value)
        # for key, value in reversed(container_rocksalt_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_anhydrite_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_dolomite.items()):
        #     print(key, value)
        #
        return container_anhydrite_upper, container_rocksalt_upper, container_potash, container_rocksalt_lower, \
               container_anhydrite_lower, container_dolomite
    #
    def create_zechstein_z3(self, thickness_z3=300, top_z=0):  # Z3 - Leine Series
        condition = False
        while condition == False:
            counter = 0
            #
            fraction_shale_pre = round(rd.uniform(5, 10), 4)
            fraction_shale = round(round(fraction_shale_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_pre = round(rd.uniform(15, 25), 4)
            fraction_anhydrite = round(round(fraction_anhydrite_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_lower_pre = round(rd.uniform(10, 20), 4)
            fraction_rocksalt_lower = round(round(fraction_rocksalt_lower_pre * 2) / 2 / 100, 4)
            fraction_potash_lower_pre = round(rd.uniform(5, 10), 4)
            fraction_potash_lower = round(round(fraction_potash_lower_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_medium_pre = round(rd.uniform(25, 35), 4)
            fraction_rocksalt_medium = round(round(fraction_rocksalt_medium_pre * 2) / 2 / 100, 4)
            fraction_potash_upper_pre = round(rd.uniform(5, 10), 4)
            fraction_potash_upper = round(round(fraction_potash_upper_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_upper = round(
                1 - fraction_shale - fraction_anhydrite - fraction_rocksalt_lower - fraction_potash_lower - fraction_rocksalt_medium - fraction_potash_upper, 4)
            #
            thickness_shale = round(thickness_z3 * fraction_shale, 4)
            thickness_anhydrite = round(thickness_z3 * fraction_anhydrite, 4)
            thickness_rocksalt_lower = round(thickness_z3 * fraction_rocksalt_lower, 4)
            thickness_potash_lower = round(thickness_z3 * fraction_potash_lower, 4)
            thickness_rocksalt_medium = round(thickness_z3 * fraction_rocksalt_medium, 4)
            thickness_potash_upper = round(thickness_z3 * fraction_potash_upper, 4)
            thickness_rocksalt_upper = round(thickness_z3 * fraction_rocksalt_upper, 4)
            list_units = ["Rock Salt (upper)", "Potash (upper)", "Rock Salt (medium)", "Potash (lower)",
                          "Rock Salt (lower)", "Anhydrite", "Shale"]
            list_thickness = [thickness_rocksalt_upper, thickness_potash_upper, thickness_rocksalt_medium,
                              thickness_potash_lower, thickness_rocksalt_lower, thickness_anhydrite, thickness_shale]
            for thickness in list_thickness:
                if thickness > 0:
                    counter += 1
                else:
                    counter += 0
            if counter == len(list_thickness) and np.sum(list_thickness) == thickness_z3:
                condition = True
        #
        # for index, unit in enumerate(list_units):
        #     print(unit, list_thickness[index])
        # print("Total Thickness:", np.sum(list_thickness))
        # self.actual_thickness += thickness_dolomite / self.resolution
        actual_top = top_z
        actual_bottom = top_z + thickness_rocksalt_upper
        # print("Anhydrite (upper):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Rock Salt (upper) Unit
        container_rocksalt_upper = {}
        steps_rocksalt_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rocksalt_upper:
            depth = round(i, 4)
            container_rocksalt_upper[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt_upper
        actual_bottom += thickness_potash_upper
        # print("Rock Salt (upper):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Potash (upper) Unit
        container_potash_upper = {}
        steps_potash_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_potash_upper:
            depth = round(i, 4)
            container_potash_upper[depth] = Evaporites(fluid="water", actualThickness=0).create_potash(
                porosity=[0.0, 0.05])
        actual_top += thickness_potash_upper
        actual_bottom += thickness_rocksalt_medium
        # print("Potash:", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Rock Salt (medium Unit
        container_rocksalt_medium = {}
        steps_rocksalt_medium = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rocksalt_medium:
            depth = round(i, 4)
            container_rocksalt_medium[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt_medium
        actual_bottom += thickness_potash_lower
        # print("Rock Salt (lower):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Potash (lower) Unit
        container_potash_lower = {}
        steps_potash_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_potash_lower:
            depth = round(i, 4)
            container_potash_lower[depth] = Evaporites(fluid="water", actualThickness=0).create_potash(
                porosity=[0.0, 0.05])
        actual_top += thickness_potash_lower
        actual_bottom += thickness_rocksalt_lower
        # print("Anhydrite (lower):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Rock Salt (lower) Unit
        container_rocksalt_lower = {}
        steps_rocksalt_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rocksalt_lower:
            depth = round(i, 4)
            container_rocksalt_lower[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.1])
        actual_top += thickness_rocksalt_lower
        actual_bottom += thickness_anhydrite
        # print("Dolomite:", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Anhydrite Unit
        container_anhydrite = {}
        steps_anhydrite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite:
            depth = round(i, 4)
            container_anhydrite[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite
        actual_bottom += thickness_shale
        #
        ## Create Shale Unit
        container_shale = {}
        steps_shale = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_shale:
            depth = round(i, 4)
            container_shale[depth] = Sandstone().create_mudstone(number=1, porosity=[0.0, 0.1])
            # container_shale[depth] = shale(fluid="water").create_simple_shale(dict_output=True,
            #                                                                   porosity=rd.uniform(0, 0.1))
        # actual_top += thickness_shale
        # actual_bottom += thickness_shale
        #
        # self.actual_thickness = thickness_dolomite + thickness_anhydrite_lower + thickness_rocksalt_lower \
        #                         + thickness_potash + thickness_rocksalt_upper + thickness_anhydrite_upper
        #
        ## TEST
        # for key, value in reversed(container_anhydrite_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_rocksalt_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_potash.items()):
        #     print(key, value)
        # for key, value in reversed(container_rocksalt_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_anhydrite_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_dolomite.items()):
        #     print(key, value)
        #
        return container_rocksalt_upper, container_potash_upper, container_rocksalt_medium, container_potash_lower, \
               container_rocksalt_lower, container_anhydrite, container_shale
    #
    def create_zechstein_z4(self, thickness_z4=200, top_z=0):  # Z4 - Aller Series
        condition = False
        while condition == False:
            counter = 0
            #
            fraction_shale_pre = round(rd.uniform(7, 15), 4)
            fraction_shale = round(round(fraction_shale_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_lower_pre = round(rd.uniform(7, 15), 4)
            fraction_anhydrite_lower = round(round(fraction_anhydrite_lower_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_pre = round(rd.uniform(63, 71), 4)
            fraction_rocksalt = round(round(fraction_rocksalt_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_upper = round(1 - fraction_shale - fraction_anhydrite_lower - fraction_rocksalt, 4)
            #
            thickness_shale = round(thickness_z4 * fraction_shale, 4)
            thickness_anhydrite_lower = round(thickness_z4 * fraction_anhydrite_lower, 4)
            thickness_rocksalt = round(thickness_z4 * fraction_rocksalt, 4)
            thickness_anhydrite_upper = round(thickness_z4 * fraction_anhydrite_upper, 4)
            list_units = ["Anhydrite (upper)", "Rock Salt", "Anhydrite (lower)", "Shale"]
            list_thickness = [thickness_anhydrite_upper, thickness_rocksalt, thickness_anhydrite_lower, thickness_shale]
            for thickness in list_thickness:
                if thickness > 0:
                    counter += 1
                else:
                    counter += 0
            if counter == len(list_thickness) and np.sum(list_thickness) == thickness_z4:
                condition = True
        #
        # for index, unit in enumerate(list_units):
        #     print(unit, list_thickness[index])
        # print("Total Thickness:", np.sum(list_thickness))
        # self.actual_thickness += thickness_dolomite / self.resolution
        actual_top = top_z
        actual_bottom = top_z + thickness_anhydrite_upper
        # print("Anhydrite (upper):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Anhydrite (upper) Unit
        container_anhydrite_upper = {}
        steps_anhydrite_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite_upper:
            depth = round(i, 4)
            container_anhydrite_upper[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_upper
        actual_bottom += thickness_rocksalt
        # print("Rock Salt (upper):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Rock Salt Unit
        container_rocksalt = {}
        steps_rocksalt_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rocksalt_upper:
            depth = round(i, 4)
            container_rocksalt[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt
        actual_bottom += thickness_anhydrite_lower
        # print("Potash:", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Anhydrite (lower) Unit
        container_anhydrite_lower = {}
        steps_anhydrite_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite_lower:
            depth = round(i, 4)
            container_anhydrite_lower[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_lower
        actual_bottom += thickness_shale
        # print("Rock Salt (lower):", "Top =", actual_top, "Bottom =", actual_bottom, "Thickness =",
        #       actual_bottom - actual_top)
        #
        ## Create Shale Unit
        container_shale = {}
        steps_shale = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_shale:
            depth = round(i, 4)
            container_shale[depth] = Sandstone().create_mudstone(number=1, porosity=[0.0, 0.1])
            # container_shale[depth] = shale(fluid="water").create_simple_shale(dict_output=True,
            #                                                                   porosity=rd.uniform(0, 0.1))
        # actual_top += thickness_shale
        # actual_bottom += thickness_shale
        #
        # self.actual_thickness = thickness_dolomite + thickness_anhydrite_lower + thickness_rocksalt_lower \
        #                         + thickness_potash + thickness_rocksalt_upper + thickness_anhydrite_upper
        #
        ## TEST
        # for key, value in reversed(container_anhydrite_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_rocksalt_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_potash.items()):
        #     print(key, value)
        # for key, value in reversed(container_rocksalt_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_anhydrite_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_dolomite.items()):
        #     print(key, value)
        #
        return container_anhydrite_upper, container_rocksalt, container_anhydrite_lower, container_shale