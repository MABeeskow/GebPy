#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		series.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		08.08.2023

#-----------------------------------------------

## MODULES
import random as rd
import numpy as np
from modules.ore import Ores
from modules.carbonates import CarbonateRocks
from modules.evaporites import Evaporites
from modules.siliciclastics import SiliciclasticRocks
from modules.sedimentary_rocks import SedimentaryRocks

#######################
## SERIES GENERATION ##
#######################
#
## ROTLIEGEND
#
## ZECHSTEIN
class Zechstein:
    #
    def __init__(self, actual_thickness=0, thickness=1000, resolution=25, composition=None):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
        #
        self.actual_thickness = actual_thickness
    #
    def export_lithological_keys(self):
        list_keys = ["Kupferschiefer", "Limestone", "Anhydrite", "Dolostone", "Rock Salt", "Potash", "Mudstone"]
        list_keys.sort()
        #
        return list_keys
    #
    def create_zechstein_z1(self, thickness_z1=100, top_z=0):  # Z1 - Werra Series
        fraction_kupferschiefer_pre = round(rd.uniform(11, 18), 4)
        fraction_kupferschiefer = round(round(fraction_kupferschiefer_pre*2)/2/100, 4)
        fraction_limestone_pre = round(rd.uniform(21, 36), 4)
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
            container_limestone[depth] = CarbonateRocks(fluid="water", actualThickness=0).create_limestone(
                number=1, porosity=[0.0, 0.4])
        actual_top += thickness_limestone
        actual_bottom += thickness_kupferschiefer
        #
        ## Create Kupferschiefer Unit
        container_kupferschiefer = {}
        steps_kupferschiefer = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_kupferschiefer:
            depth = round(i, 4)
            container_kupferschiefer[depth] = Ores(
                fluid="water", actualThickness=0, porosity=rd.uniform(0.0, 0.05),
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
            fraction_dolomite_pre = round(rd.uniform(6, 10), 4)
            fraction_dolomite = round(round(fraction_dolomite_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_lower_pre = round(rd.uniform(14, 24), 4)
            fraction_anhydrite_lower = round(round(fraction_anhydrite_lower_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_lower_pre = round(rd.uniform(42, 71), 4)
            fraction_rocksalt_lower = round(round(fraction_rocksalt_lower_pre * 2) / 2 / 100, 4)
            fraction_potash_pre = round(rd.uniform(8, 13), 4)
            fraction_potash = round(round(fraction_potash_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_upper_pre = round(rd.uniform(3, 5), 4)
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
                          "Dolostone"]
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
        self.actual_thickness += thickness_dolomite/self.resolution
        actual_top = top_z
        actual_bottom = top_z + thickness_anhydrite_upper
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
        #
        ## Create Dolostone Unit
        container_dolomite = {}
        steps_dolomite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_dolomite:
            depth = round(i, 4)
            container_dolomite[depth] = CarbonateRocks().create_dolostone(number=1, porosity=[0.1, 0.4])
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
            fraction_mudstone_pre = round(rd.uniform(5, 8), 4)
            fraction_mudstone = round(round(fraction_mudstone_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_pre = round(rd.uniform(14, 24), 4)
            fraction_anhydrite = round(round(fraction_anhydrite_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_lower_pre = round(rd.uniform(13, 22), 4)
            fraction_rocksalt_lower = round(round(fraction_rocksalt_lower_pre * 2) / 2 / 100, 4)
            fraction_potash_lower_pre = round(rd.uniform(4, 7), 4)
            fraction_potash_lower = round(round(fraction_potash_lower_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_medium_pre = round(rd.uniform(25, 42), 4)
            fraction_rocksalt_medium = round(round(fraction_rocksalt_medium_pre * 2) / 2 / 100, 4)
            fraction_potash_upper_pre = round(rd.uniform(4, 6), 4)
            fraction_potash_upper = round(round(fraction_potash_upper_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_upper = round(
                1 - fraction_mudstone - fraction_anhydrite - fraction_rocksalt_lower - fraction_potash_lower -
                fraction_rocksalt_medium - fraction_potash_upper, 4)
            #
            thickness_mudstone = round(thickness_z3 * fraction_mudstone, 4)
            thickness_anhydrite = round(thickness_z3 * fraction_anhydrite, 4)
            thickness_rocksalt_lower = round(thickness_z3 * fraction_rocksalt_lower, 4)
            thickness_potash_lower = round(thickness_z3 * fraction_potash_lower, 4)
            thickness_rocksalt_medium = round(thickness_z3 * fraction_rocksalt_medium, 4)
            thickness_potash_upper = round(thickness_z3 * fraction_potash_upper, 4)
            thickness_rocksalt_upper = round(thickness_z3 * fraction_rocksalt_upper, 4)
            list_units = ["Rock Salt (upper)", "Potash (upper)", "Rock Salt (medium)", "Potash (lower)",
                          "Rock Salt (lower)", "Anhydrite", "Mudstone"]
            list_thickness = [thickness_rocksalt_upper, thickness_potash_upper, thickness_rocksalt_medium,
                              thickness_potash_lower, thickness_rocksalt_lower, thickness_anhydrite, thickness_mudstone]
            for thickness in list_thickness:
                if thickness > 0:
                    counter += 1
                else:
                    counter += 0
            if counter == len(list_thickness) and np.sum(list_thickness) == thickness_z3:
                condition = True
        #
        actual_top = top_z
        actual_bottom = top_z + thickness_rocksalt_upper
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
        #
        ## Create Anhydrite Unit
        container_anhydrite = {}
        steps_anhydrite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite:
            depth = round(i, 4)
            container_anhydrite[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite
        actual_bottom += thickness_mudstone
        #
        ## Create Mudstone Unit
        container_mudstone = {}
        steps_mudstone = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_mudstone:
            depth = round(i, 4)
            container_mudstone[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
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
               container_rocksalt_lower, container_anhydrite, container_mudstone
    #
    def create_zechstein_z4(self, thickness_z4=200, top_z=0):  # Z4 - Aller Series
        condition = False
        while condition == False:
            counter = 0
            #
            fraction_mudstone_pre = round(rd.uniform(10, 16), 4)
            fraction_mudstone = round(round(fraction_mudstone_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_lower_pre = round(rd.uniform(5, 9), 4)
            fraction_anhydrite_lower = round(round(fraction_anhydrite_lower_pre * 2) / 2 / 100, 4)
            fraction_rocksalt_pre = round(rd.uniform(55, 92), 4)
            fraction_rocksalt = round(round(fraction_rocksalt_pre * 2) / 2 / 100, 4)
            fraction_anhydrite_upper = round(1 - fraction_mudstone - fraction_anhydrite_lower - fraction_rocksalt, 4)
            #
            thickness_mudstone = round(thickness_z4 * fraction_mudstone, 4)
            thickness_anhydrite_lower = round(thickness_z4 * fraction_anhydrite_lower, 4)
            thickness_rocksalt = round(thickness_z4 * fraction_rocksalt, 4)
            thickness_anhydrite_upper = round(thickness_z4 * fraction_anhydrite_upper, 4)
            list_units = ["Anhydrite (upper)", "Rock Salt", "Anhydrite (lower)", "Shale"]
            list_thickness = [thickness_anhydrite_upper, thickness_rocksalt, thickness_anhydrite_lower, thickness_mudstone]
            for thickness in list_thickness:
                if thickness > 0:
                    counter += 1
                else:
                    counter += 0
            if counter == len(list_thickness) and np.sum(list_thickness) == thickness_z4:
                condition = True
        #
        actual_top = top_z
        actual_bottom = top_z + thickness_anhydrite_upper
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
        #
        ## Create Anhydrite (lower) Unit
        container_anhydrite_lower = {}
        steps_anhydrite_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite_lower:
            depth = round(i, 4)
            container_anhydrite_lower[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_lower
        actual_bottom += thickness_mudstone
        #
        ## Create Mudstone Unit
        container_mudstone = {}
        steps_mudstone = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_mudstone:
            depth = round(i, 4)
            container_mudstone[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
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
        return container_anhydrite_upper, container_rocksalt, container_anhydrite_lower, container_mudstone
    #
    def create_zechstein_z5(self, thickness_z5=100, top_z=0):  # 5 - Ohre Series
        fraction_mudstone_lower_pre = round(rd.uniform(22, 36), 4)
        fraction_mudstone_lower = round(round(fraction_mudstone_lower_pre*2)/2/100, 4)
        fraction_anhydrite_pre = round(rd.uniform(11, 19), 4)
        fraction_anhydrite = round(round(fraction_anhydrite_pre*2)/2/100, 4)
        fraction_mudstone_upper = round(1 - fraction_mudstone_lower - fraction_anhydrite, 4)
        #
        thickness_mudstone_lower = round(thickness_z5*fraction_mudstone_lower, 4)
        thickness_anhydrite = round(thickness_z5*fraction_anhydrite, 4)
        thickness_mudstone_upper = round(thickness_z5*fraction_mudstone_upper, 4)
        #
        actual_top = top_z
        actual_bottom = top_z + thickness_mudstone_upper
        #
        ## Create Mudstone Upper Unit
        container_mudstone_upper = {}
        steps_mudstone_upper = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_mudstone_upper:
            depth = round(i, 4)
            container_mudstone_upper[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_upper
        actual_bottom += thickness_anhydrite
        #
        ## Create Anhydrite Unit
        container_anhydrite = {}
        steps_anhydrite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_anhydrite:
            depth = round(i, 4)
            container_anhydrite[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite
        actual_bottom += thickness_mudstone_lower
        #
        ## Create Mudstone Lower Unit
        container_mudstone_lower = {}
        steps_mudstone_lower = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_mudstone_lower:
            depth = round(i, 4)
            container_mudstone_lower[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_upper.items()):
        #     print(key, value)
        # for key, value in reversed(container_anhydrite.items()):
        #     print(key, value)
        # for key, value in reversed(container_mudstone_lower.items()):
        #     print(key, value)
        #
        return container_mudstone_upper, container_anhydrite, container_mudstone_lower
    #
class Muschelkalk:
    #
    def __init__(self, actual_thickness=0, thickness=1000, resolution=25, composition=None):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
        #
        self.actual_thickness = actual_thickness
    #
    def export_lithological_keys(self):
        list_keys = ["Marl", "Dolostone", "Limestone", "Anhydrite", "Mudstone"]
        list_keys.sort()
        #
        return list_keys
    #
    def create_muschelkalk_unterer(self, thickness_unit=100, top_unit=0):  # Unterer Muschelkalk
        fraction_marl_pre = round(rd.uniform(10, 15), 4)
        fraction_marl = round(round(fraction_marl_pre * 2) / 2 / 100, 4)
        fraction_dolomite_pre = round(rd.uniform(30, 40), 4)
        fraction_dolomite = round(round(fraction_dolomite_pre * 2) / 2 / 100, 4)
        fraction_limestone = round(1 - fraction_marl - fraction_dolomite, 4)
        #
        thickness_marl = round(thickness_unit * fraction_marl, 4)
        thickness_dolomite = round(thickness_unit * fraction_dolomite, 4)
        thickness_limestone = round(thickness_unit * fraction_limestone, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_limestone
        #
        ## Create Limestone Unit
        container_limestone = {}
        steps_limestone = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_limestone:
            depth = round(i, 4)
            # container_limestone[depth] = limestone(fluid="water", actualThickness=0).create_simple_limestone(
            #     dict=True, porosity=rd.uniform(0.15, 0.4))
            container_limestone[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.0, 0.4])
        actual_top += thickness_limestone
        actual_bottom += thickness_dolomite
        #
        ## Create Dolostone Unit
        container_dolomite = {}
        steps_dolomite = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_dolomite:
            depth = round(i, 4)
            container_dolomite[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_dolostone(number=1, porosity=[0.1, 0.2])
        actual_top += thickness_dolomite
        actual_bottom += thickness_marl
        #
        ## Create Marl Unit
        container_marl = {}
        steps_marl = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_marl:
            depth = round(i, 4)
            container_marl[depth] = SedimentaryRocks(
                fluid="water", actualThickness=0).create_marl_alt(
                number=1, porosity=rd.uniform(0.1, 0.3))
        #
        ## TEST
        # for key, value in reversed(container_anhydrite.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_limestone, container_dolomite, container_marl
    #
    def create_muschelkalk_mittlerer(self, thickness_unit=100, top_unit=0):  # Mittlerer Muschelkalk
        fraction_dolomite_lower_pre = round(rd.uniform(10, 15), 4)
        fraction_dolomite_lower = round(round(fraction_dolomite_lower_pre * 2) / 2 / 100, 4)
        fraction_anhydrite_pre = round(rd.uniform(15, 25), 4)
        fraction_anhydrite = round(round(fraction_anhydrite_pre * 2) / 2 / 100, 4)
        fraction_dolomite_medium_pre = round(rd.uniform(30, 40), 4)
        fraction_dolomite_medium = round(round(fraction_dolomite_medium_pre * 2) / 2 / 100, 4)
        fraction_marl_pre = round(rd.uniform(5, 10), 4)
        fraction_marl = round(round(fraction_marl_pre * 2) / 2 / 100, 4)
        fraction_dolomite_upper = round(
            1 - fraction_dolomite_lower - fraction_anhydrite - fraction_dolomite_medium - fraction_marl, 4)
        #
        thickness_dolomite_upper = round(thickness_unit * fraction_dolomite_upper, 4)
        thickness_marl = round(thickness_unit * fraction_marl, 4)
        thickness_dolomite_medium = round(thickness_unit * fraction_dolomite_medium, 4)
        thickness_anhydrite = round(thickness_unit * fraction_anhydrite, 4)
        thickness_dolomite_lower = round(thickness_unit * fraction_dolomite_lower, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_dolomite_upper
        #
        ## Create Dolostone upper Unit
        container_dolomite_upper = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_dolomite_upper[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_dolostone(number=1, porosity=[0.1, 0.2])
        actual_top += thickness_dolomite_upper
        actual_bottom += thickness_marl
        #
        ## Create Marl Unit
        container_marl = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_marl[depth] = SedimentaryRocks(
                fluid="water", actualThickness=0).create_marl_alt(
                number=1, porosity=rd.uniform(0.1, 0.3))
        actual_top += thickness_marl
        actual_bottom += thickness_dolomite_medium
        #
        ## Create Dolostone medium Unit
        container_dolomite_medium = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_dolomite_medium[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_dolostone(number=1, porosity=[0.1, 0.2])
        actual_top += thickness_dolomite_medium
        actual_bottom += thickness_anhydrite
        #
        ## Create Anhydrite Unit
        container_anhydrite = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_anhydrite[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite
        actual_bottom += thickness_dolomite_lower
        #
        ## Create Dolostone lower Unit
        container_dolomite_lower = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_dolomite_lower[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_dolostone(number=1, porosity=[0.1, 0.2])
        #
        ## TEST
        # for key, value in reversed(container_anhydrite.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_dolomite_upper, container_marl, container_dolomite_medium, container_anhydrite, \
               container_dolomite_lower
    #
    def create_muschelkalk_oberer(self, thickness_unit=100, top_unit=0):  # Oberer Muschelkalk
        fraction_limestone_upper_pre = round(rd.uniform(10, 15), 4)
        fraction_limestone_upper = round(round(fraction_limestone_upper_pre * 2) / 2 / 100, 4)
        fraction_mudstone_upper_pre = round(rd.uniform(5, 10), 4)
        fraction_mudstone_upper = round(round(fraction_mudstone_upper_pre * 2) / 2 / 100, 4)
        fraction_limestone_medium_upper_pre = round(rd.uniform(10, 15), 4)
        fraction_limestone_medium_upper = round(round(fraction_limestone_medium_upper_pre * 2) / 2 / 100, 4)
        fraction_marl_pre = round(rd.uniform(5, 10), 4)
        fraction_marl = round(round(fraction_marl_pre * 2) / 2 / 100, 4)
        fraction_limestone_medium_medium_pre = round(rd.uniform(10, 15), 4)
        fraction_limestone_medium_medium = round(round(fraction_limestone_medium_medium_pre * 2) / 2 / 100, 4)
        fraction_mudstone_medium_pre = round(rd.uniform(5, 10), 4)
        fraction_mudstone_medium = round(round(fraction_mudstone_medium_pre * 2) / 2 / 100, 4)
        fraction_limestone_medium_lower_pre = round(rd.uniform(10, 15), 4)
        fraction_limestone_medium_lower = round(round(fraction_limestone_medium_lower_pre * 2) / 2 / 100, 4)
        fraction_mudstone_lower_pre = round(rd.uniform(5, 10), 4)
        fraction_mudstone_lower = round(round(fraction_mudstone_lower_pre * 2) / 2 / 100, 4)
        fraction_limestone_lower = round(
            1 - fraction_limestone_upper - fraction_mudstone_upper - fraction_limestone_medium_upper - fraction_marl
            - fraction_limestone_medium_medium - fraction_mudstone_medium - fraction_limestone_medium_lower
            - fraction_mudstone_lower, 4)
        #
        thickness_limestone_upper = round(thickness_unit * fraction_limestone_upper, 4)
        thickness_mudstone_upper = round(thickness_unit * fraction_mudstone_upper, 4)
        thickness_limestone_medium_upper = round(thickness_unit * fraction_limestone_medium_upper, 4)
        thickness_marl = round(thickness_unit * fraction_marl, 4)
        thickness_limestone_medium_medium = round(thickness_unit * fraction_limestone_medium_medium, 4)
        thickness_mudstone_medium = round(thickness_unit * fraction_mudstone_medium, 4)
        thickness_limestone_medium_lower = round(thickness_unit * fraction_limestone_medium_lower, 4)
        thickness_mudstone_lower = round(thickness_unit * fraction_mudstone_lower, 4)
        thickness_limestone_lower = round(thickness_unit * fraction_limestone_lower, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_limestone_upper
        #
        ## Create Limestone upper Unit
        container_limestone_upper = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_limestone_upper[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.0, 0.4])
        actual_top += thickness_limestone_upper
        actual_bottom += thickness_mudstone_upper
        #
        ## Create Mudstone upper Unit
        container_mudstone_upper = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_mudstone_upper[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_upper
        actual_bottom += thickness_limestone_medium_upper
        #
        ## Create Limestone medium upper Unit
        container_limestone_medium_upper = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_limestone_medium_upper[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.0, 0.4])
        actual_top += thickness_limestone_medium_upper
        actual_bottom += thickness_marl
        #
        ## Create Marl Unit
        container_marl = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_marl[depth] = SedimentaryRocks(
                fluid="water", actualThickness=0).create_marl_alt(
                number=1, porosity=rd.uniform(0.1, 0.3))
        actual_top += thickness_marl
        actual_bottom += thickness_limestone_medium_medium
        #
        ## Create Limestone medium medium Unit
        container_limestone_medium_medium = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_limestone_medium_medium[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.0, 0.4])
        actual_top += thickness_limestone_medium_medium
        actual_bottom += thickness_mudstone_medium
        #
        ## Create Mudstone medium Unit
        container_mudstone_medium = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_mudstone_medium[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_medium
        actual_bottom += thickness_limestone_medium_lower
        #
        ## Create Limestone medium lower Unit
        container_limestone_medium_lower = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_limestone_medium_lower[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.0, 0.4])
        actual_top += thickness_limestone_medium_lower
        actual_bottom += thickness_mudstone_lower
        #
        ## Create Mudstone lower Unit
        container_mudstone_lower = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_mudstone_lower[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_lower
        actual_bottom += thickness_limestone_lower
        #
        ## Create Limestone lower Unit
        container_limestone_lower = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_limestone_lower[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.0, 0.4])
        #
        ## TEST
        # for key, value in reversed(container_limestone_lower.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_limestone_upper, container_mudstone_upper, container_limestone_medium_upper, container_marl, \
               container_limestone_medium_medium, container_mudstone_medium, container_limestone_medium_lower, \
               container_mudstone_lower, container_limestone_lower
    #
class Buntsandstein:
    #
    def __init__(self, actual_thickness=0, thickness=1000, resolution=25, composition=None):
        self.thickness = thickness
        self.resolution = resolution
        self.composition = composition
        #
        self.actual_thickness = actual_thickness
    #
    def export_lithological_keys(self):
        list_keys = ["Sandstone", "Mudstone", "Limestone", "Rock Salt", "Anhydrite"]
        list_keys.sort()
        #
        return list_keys
    #
    def create_buntsandstein_lower(self, thickness_unit=100, top_unit=0):
        thickness_bemburg = int(rd.uniform(0.42, 0.46)*thickness_unit)
        thickness_nordhausen = int(thickness_unit - thickness_bemburg)
        #
        data_bemburg = self.create_bemburg_series(
            top_unit=top_unit, thickness_unit=thickness_bemburg)
        data_nordhausen = self.create_nordhausen_series(
            top_unit=top_unit + thickness_bemburg, thickness_unit=thickness_nordhausen)
        #
        data_units = data_bemburg + data_nordhausen
        #
        return data_units
    #
    def create_buntsandstein_medium(self, thickness_unit=100, top_unit=0):
        thickness_solling = int(rd.uniform(0.05, 0.09)*thickness_unit)
        thickness_hardegsen = int(rd.uniform(0.19, 0.23)*thickness_unit)
        thickness_detfurth = int(rd.uniform(0.31, 0.35)*thickness_unit)
        thickness_volpriehausen = int(thickness_unit - thickness_solling - thickness_hardegsen - thickness_detfurth)
        #
        data_solling = self.create_solling_series(
            top_unit=top_unit, thickness_unit=thickness_solling)
        data_hardegsen = self.create_hardegsen_series(
            top_unit=top_unit + thickness_solling, thickness_unit=thickness_hardegsen)
        data_detfurth = self.create_detfurth_series(
            top_unit=top_unit + thickness_solling + thickness_hardegsen, thickness_unit=thickness_detfurth)
        data_volpriehausen = self.create_detfurth_series(
            top_unit=top_unit + thickness_solling + thickness_hardegsen + thickness_detfurth,
            thickness_unit=thickness_volpriehausen)
        #
        data_units = data_solling + data_hardegsen + data_detfurth + data_volpriehausen
        #
        return data_units
    #
    def create_buntsandstein_upper(self, thickness_unit=100, top_unit=0):
        thickness_myophorien = int(rd.uniform(0.10, 0.14)*thickness_unit)
        thickness_pelitroet = int(rd.uniform(0.44, 0.48)*thickness_unit)
        thickness_salinarroet = int(thickness_unit - thickness_myophorien - thickness_pelitroet)
        #
        data_myophorien = self.create_myophorien_series(
            top_unit=top_unit, thickness_unit=thickness_myophorien)
        data_pelitroet = self.create_pelitroet_series(
            top_unit=top_unit + thickness_myophorien, thickness_unit=thickness_pelitroet)
        data_salinarroet = self.create_salinarroet_series(
            top_unit=top_unit + thickness_myophorien + thickness_pelitroet, thickness_unit=thickness_salinarroet)
        #
        data_units = data_myophorien + data_pelitroet + data_salinarroet
        #
        return data_units
    #
    def create_nordhausen_series(self, thickness_unit=100, top_unit=0):
        fraction_sandstone_05 = round(rd.uniform(0.08, 0.10), 4)
        fraction_mudstone_05 = round(rd.uniform(0.11, 0.13), 4)
        fraction_sandstone_04 = round(rd.uniform(0.06, 0.08), 4)
        fraction_mudstone_04 = round(rd.uniform(0.06, 0.08), 4)
        fraction_sandstone_03 = round(rd.uniform(0.18, 0.20), 4)
        fraction_mudstone_03 = round(rd.uniform(0.06, 0.08), 4)
        fraction_sandstone_02 = round(rd.uniform(0.04, 0.06), 4)
        fraction_mudstone_02 = round(rd.uniform(0.04, 0.06), 4)
        fraction_sandstone_01 = round(rd.uniform(0.04, 0.06), 4)
        fraction_mudstone_01 = round(rd.uniform(0.04, 0.06), 4)
        #
        fraction_mudstone_06 = round(
            1 - fraction_mudstone_01 - fraction_mudstone_02 - fraction_mudstone_03 - fraction_mudstone_04 -
            fraction_mudstone_05 - fraction_sandstone_01 - fraction_sandstone_02 - fraction_sandstone_03 -
            fraction_sandstone_04 - fraction_sandstone_05, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_mudstone_02 = round(thickness_unit*fraction_mudstone_02, 4)
        thickness_mudstone_03 = round(thickness_unit*fraction_mudstone_03, 4)
        thickness_mudstone_04 = round(thickness_unit*fraction_mudstone_04, 4)
        thickness_mudstone_05 = round(thickness_unit*fraction_mudstone_05, 4)
        thickness_mudstone_06 = round(thickness_unit*fraction_mudstone_06, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        thickness_sandstone_02 = round(thickness_unit*fraction_sandstone_02, 4)
        thickness_sandstone_03 = round(thickness_unit*fraction_sandstone_03, 4)
        thickness_sandstone_04 = round(thickness_unit*fraction_sandstone_04, 4)
        thickness_sandstone_05 = round(thickness_unit*fraction_sandstone_05, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_mudstone_06
        #
        thickness_list = [thickness_mudstone_06, thickness_sandstone_05, thickness_mudstone_05, thickness_sandstone_04,
                          thickness_mudstone_04, thickness_sandstone_03, thickness_mudstone_03, thickness_sandstone_02,
                          thickness_mudstone_02, thickness_sandstone_01, thickness_mudstone_01]
        #
        ## Create Mudstone 06
        container_mudstone_06 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_06[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_06
        actual_bottom += thickness_sandstone_05
        #
        ## Create Sandstone 05
        container_sandstone_05 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_05[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.15, 0.25])
        actual_top += thickness_sandstone_05
        actual_bottom += thickness_mudstone_05
        #
        ## Create Mudstone 05
        container_mudstone_05 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_05[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_05
        actual_bottom += thickness_sandstone_04
        #
        ## Create Sandstone 04
        container_sandstone_04 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_04[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.15, 0.25])
        actual_top += thickness_sandstone_04
        actual_bottom += thickness_mudstone_04
        #
        ## Create Mudstone 04
        container_mudstone_04 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_04[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_04
        actual_bottom += thickness_sandstone_03
        #
        ## Create Sandstone 03
        container_sandstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_03[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.15, 0.25])
        actual_top += thickness_sandstone_03
        actual_bottom += thickness_mudstone_03
        #
        ## Create Mudstone 03
        container_mudstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_03[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_03
        actual_bottom += thickness_sandstone_02
        #
        ## Create Sandstone 02
        container_sandstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_02[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.15, 0.25])
        actual_top += thickness_sandstone_02
        actual_bottom += thickness_mudstone_02
        #
        ## Create Mudstone 02
        container_mudstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_02[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_02
        actual_bottom += thickness_sandstone_01
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.2, 0.3])
        actual_top += thickness_sandstone_01
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_mudstone_06, container_sandstone_05, container_mudstone_05, container_sandstone_04, \
            container_mudstone_04, container_sandstone_03, container_mudstone_03, container_sandstone_02, \
            container_mudstone_02, container_sandstone_01, container_mudstone_01
    #
    def create_bemburg_series(self, thickness_unit=100, top_unit=0):
        fraction_sandstone_06 = round(rd.uniform(0.09, 0.11), 4)
        fraction_mudstone_03 = round(rd.uniform(0.10, 0.12), 4)
        fraction_sandstone_05 = round(rd.uniform(0.09, 0.11), 4)
        fraction_sandstone_04 = round(rd.uniform(0.06, 0.08), 4)
        fraction_sandstone_03 = round(rd.uniform(0.09, 0.11), 4)
        fraction_mudstone_02 = round(rd.uniform(0.09, 0.11), 4)
        fraction_sandstone_02 = round(rd.uniform(0.08, 0.10), 4)
        fraction_mudstone_01 = round(rd.uniform(0.05, 0.07), 4)
        fraction_sandstone_01 = round(rd.uniform(0.11, 0.13), 4)
        #
        fraction_mudstone_04 = round(
            1 - fraction_mudstone_01 - fraction_mudstone_02 - fraction_mudstone_03 - fraction_sandstone_01 -
            fraction_sandstone_02 - fraction_sandstone_03 - fraction_sandstone_04 - fraction_sandstone_05 -
            fraction_sandstone_06, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_mudstone_02 = round(thickness_unit*fraction_mudstone_02, 4)
        thickness_mudstone_03 = round(thickness_unit*fraction_mudstone_03, 4)
        thickness_mudstone_04 = round(thickness_unit*fraction_mudstone_04, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        thickness_sandstone_02 = round(thickness_unit*fraction_sandstone_02, 4)
        thickness_sandstone_03 = round(thickness_unit*fraction_sandstone_03, 4)
        thickness_sandstone_04 = round(thickness_unit*fraction_sandstone_04, 4)
        thickness_sandstone_05 = round(thickness_unit*fraction_sandstone_05, 4)
        thickness_sandstone_06 = round(thickness_unit*fraction_sandstone_06, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_mudstone_04
        #
        thickness_list = [thickness_mudstone_04, thickness_sandstone_06, thickness_mudstone_03, thickness_sandstone_05,
                          thickness_sandstone_04, thickness_sandstone_03, thickness_mudstone_02, thickness_sandstone_02,
                          thickness_mudstone_01, thickness_sandstone_01]
        #
        ## Create Mudstone 04
        container_mudstone_04 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_04[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_04
        actual_bottom += thickness_sandstone_06
        #
        ## Create Sandstone 06
        container_sandstone_06 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_06[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_06
        actual_bottom += thickness_mudstone_03
        #
        ## Create Mudstone 03
        container_mudstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_03[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_03
        actual_bottom += thickness_sandstone_05
        #
        ## Create Sandstone 05
        container_sandstone_05 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_05[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_05
        actual_bottom += thickness_sandstone_04
        #
        ## Create Sandstone 04
        container_sandstone_04 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_04[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.2, 0.3])
        actual_top += thickness_sandstone_04
        actual_bottom += thickness_sandstone_03
        #
        ## Create Sandstone 03
        container_sandstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_03[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_03
        actual_bottom += thickness_mudstone_02
        #
        ## Create Mudstone 02
        container_mudstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_02[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_02
        actual_bottom += thickness_sandstone_02
        #
        ## Create Sandstone 02
        container_sandstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_02[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_02
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_01
        actual_bottom += thickness_sandstone_01
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_mudstone_04, container_sandstone_06, container_mudstone_03, container_sandstone_05, \
            container_sandstone_04, container_sandstone_03, container_mudstone_02, container_sandstone_02, \
            container_mudstone_01, container_sandstone_01
    #
    def create_volpriehausen_series(self, thickness_unit=100, top_unit=0):
        fraction_sandstone_03 = round(rd.uniform(0.13, 0.15), 4)
        fraction_mudstone_02 = round(rd.uniform(0.22, 0.24), 4)
        fraction_sandstone_02 = round(rd.uniform(0.13, 0.15), 4)
        fraction_mudstone_01 = round(rd.uniform(0.09, 0.11), 4)
        fraction_sandstone_01 = round(rd.uniform(0.28, 0.30), 4)
        #
        fraction_mudstone_03 = round(
            1 - fraction_mudstone_01 - fraction_mudstone_02 - fraction_sandstone_01 - fraction_sandstone_02 -
            fraction_sandstone_03, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_mudstone_02 = round(thickness_unit*fraction_mudstone_02, 4)
        thickness_mudstone_03 = round(thickness_unit*fraction_mudstone_03, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        thickness_sandstone_02 = round(thickness_unit*fraction_sandstone_02, 4)
        thickness_sandstone_03 = round(thickness_unit*fraction_sandstone_03, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_mudstone_03
        #
        thickness_list = [thickness_mudstone_03, thickness_sandstone_03, thickness_mudstone_02, thickness_sandstone_02,
                          thickness_mudstone_01, thickness_sandstone_01]
        #
        ## Create Mudstone 03
        container_mudstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_03[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_03
        actual_bottom += thickness_sandstone_03
        #
        ## Create Sandstone 03
        container_sandstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_03[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_03
        actual_bottom += thickness_mudstone_02
        #
        ## Create Mudstone 02
        container_mudstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_02[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_02
        actual_bottom += thickness_sandstone_02
        #
        ## Create Sandstone 02
        container_sandstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_02[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_02
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_01
        actual_bottom += thickness_sandstone_01
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.25, 0.3])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_mudstone_03, container_sandstone_03, container_mudstone_02, container_sandstone_02, \
            container_mudstone_01, container_sandstone_01
    #
    def create_detfurth_series(self, thickness_unit=100, top_unit=0):
        fraction_sandstone_03 = round(rd.uniform(0.16, 0.18), 4)
        fraction_mudstone_02 = round(rd.uniform(0.16, 0.18), 4)
        fraction_sandstone_02 = round(rd.uniform(0.16, 0.18), 4)
        fraction_mudstone_01 = round(rd.uniform(0.18, 0.20), 4)
        fraction_sandstone_01 = round(rd.uniform(0.16, 0.18), 4)
        #
        fraction_mudstone_03 = round(
            1 - fraction_mudstone_01 - fraction_mudstone_02 - fraction_sandstone_01 - fraction_sandstone_02 -
            fraction_sandstone_03, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_mudstone_02 = round(thickness_unit*fraction_mudstone_02, 4)
        thickness_mudstone_03 = round(thickness_unit*fraction_mudstone_03, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        thickness_sandstone_02 = round(thickness_unit*fraction_sandstone_02, 4)
        thickness_sandstone_03 = round(thickness_unit*fraction_sandstone_03, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_mudstone_03
        #
        thickness_list = [thickness_mudstone_03, thickness_sandstone_03, thickness_mudstone_02, thickness_sandstone_02,
                          thickness_mudstone_01, thickness_sandstone_01]
        #
        ## Create Mudstone 03
        container_mudstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_03[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_03
        actual_bottom += thickness_sandstone_03
        #
        ## Create Sandstone 03
        container_sandstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_03[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_03
        actual_bottom += thickness_mudstone_02
        #
        ## Create Mudstone 02
        container_mudstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_02[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_02
        actual_bottom += thickness_sandstone_02
        #
        ## Create Sandstone 02
        container_sandstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_02[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_02
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_01
        actual_bottom += thickness_sandstone_01
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.25, 0.3])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_mudstone_03, container_sandstone_03, container_mudstone_02, container_sandstone_02, \
            container_mudstone_01, container_sandstone_01
    #
    def create_hardegsen_series(self, thickness_unit=100, top_unit=0):
        fraction_mudstone_01 = round(rd.uniform(0.37, 0.39), 4)
        fraction_sandstone_01 = round(rd.uniform(0.30, 0.32), 4)
        #
        fraction_sandstone_02 = round(1 - fraction_mudstone_01 - fraction_sandstone_01, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        thickness_sandstone_02 = round(thickness_unit*fraction_sandstone_02, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_sandstone_02
        #
        thickness_list = [thickness_sandstone_02, thickness_mudstone_01, thickness_sandstone_01]
        #
        ## Create Sandstone 02
        container_sandstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_02[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_02
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_01
        actual_bottom += thickness_sandstone_01
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.2, 0.3])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_sandstone_02, container_mudstone_01, container_sandstone_01
    #
    def create_solling_series(self, thickness_unit=100, top_unit=0):
        fraction_mudstone_01 = round(rd.uniform(0.32, 0.34), 4)
        #
        fraction_sandstone_01 = round(1 - fraction_mudstone_01, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_sandstone_01
        #
        thickness_list = [thickness_sandstone_01, thickness_mudstone_01]
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.25, 0.3])
        actual_top += thickness_sandstone_01
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_sandstone_01, container_mudstone_01
    #
    def create_salinarroet_series(self, thickness_unit=100, top_unit=0):
        fraction_rocksalt_02 = round(rd.uniform(0.15, 0.17), 4)
        fraction_anhydrite_02 = round(rd.uniform(0.09, 0.11), 4)
        fraction_rocksalt_01 = round(rd.uniform(0.35, 0.37), 4)
        fraction_anhydrite_01 = round(rd.uniform(0.15, 0.17), 4)
        #
        fraction_anhydrite_03 = round(1 - fraction_anhydrite_02 - fraction_anhydrite_01 - fraction_rocksalt_02 -
                                      fraction_rocksalt_01, 4)
        #
        thickness_anhydrite_01 = round(thickness_unit*fraction_anhydrite_01, 4)
        thickness_anhydrite_02 = round(thickness_unit*fraction_anhydrite_02, 4)
        thickness_anhydrite_03 = round(thickness_unit*fraction_anhydrite_03, 4)
        thickness_rocksalt_01 = round(thickness_unit*fraction_rocksalt_01, 4)
        thickness_rocksalt_02 = round(thickness_unit*fraction_rocksalt_02, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_anhydrite_03
        #
        ## Create Anhydrite 03
        container_anhydrite_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_anhydrite_03[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_03
        actual_bottom += thickness_rocksalt_02
        #
        ## Create Rock Salt 02
        container_rocksalt_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_rocksalt_02[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt_02
        actual_bottom += thickness_anhydrite_02
        #
        ## Create Anhydrite 02
        container_anhydrite_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_anhydrite_02[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_02
        actual_bottom += thickness_rocksalt_01
        #
        ## Create Rock Salt 01
        container_rocksalt_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_rocksalt_01[depth] = Evaporites(fluid="water", actualThickness=0).create_rocksalt(
                porosity=[0.0, 0.05])
        actual_top += thickness_rocksalt_01
        actual_bottom += thickness_anhydrite_01
        #
        ## Create Anhydrite 01
        container_anhydrite_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_anhydrite_01[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_anhydrite_03, container_rocksalt_02, container_anhydrite_02, container_rocksalt_01, \
            container_anhydrite_01
    #
    def create_pelitroet_series(self, thickness_unit=100, top_unit=0):
        fraction_anhydrite_03 = round(rd.uniform(0.05, 0.07), 4)
        fraction_mudstone_05 = round(rd.uniform(0.07, 0.09), 4)
        fraction_anhydrite_02 = round(rd.uniform(0.05, 0.07), 4)
        fraction_mudstone_04 = round(rd.uniform(0.14, 0.16), 4)
        fraction_sandstone_02 = round(rd.uniform(0.11, 0.13), 4)
        fraction_mudstone_03 = round(rd.uniform(0.09, 0.11), 4)
        fraction_anhydrite_01 = round(rd.uniform(0.05, 0.07), 4)
        fraction_mudstone_02 = round(rd.uniform(0.09, 0.11), 4)
        fraction_sandstone_01 = round(rd.uniform(0.11, 0.13), 4)
        fraction_mudstone_01 = round(rd.uniform(0.07, 0.09), 4)
        #
        fraction_mudstone_06 = round(
            1 - fraction_anhydrite_03 - fraction_anhydrite_02 - fraction_anhydrite_01 - fraction_mudstone_05 -
            fraction_mudstone_04 - fraction_mudstone_03 - fraction_mudstone_02 - fraction_mudstone_01 -
            fraction_sandstone_02 - fraction_sandstone_01, 4)
        #
        thickness_anhydrite_01 = round(thickness_unit*fraction_anhydrite_01, 4)
        thickness_anhydrite_02 = round(thickness_unit*fraction_anhydrite_02, 4)
        thickness_anhydrite_03 = round(thickness_unit*fraction_anhydrite_03, 4)
        thickness_sandstone_01 = round(thickness_unit*fraction_sandstone_01, 4)
        thickness_sandstone_02 = round(thickness_unit*fraction_sandstone_02, 4)
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_mudstone_02 = round(thickness_unit*fraction_mudstone_02, 4)
        thickness_mudstone_03 = round(thickness_unit*fraction_mudstone_03, 4)
        thickness_mudstone_04 = round(thickness_unit*fraction_mudstone_04, 4)
        thickness_mudstone_05 = round(thickness_unit*fraction_mudstone_05, 4)
        thickness_mudstone_06 = round(thickness_unit*fraction_mudstone_06, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_mudstone_06
        #
        ## Create Mudstone 06
        container_mudstone_06 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_06[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_06
        actual_bottom += thickness_anhydrite_03
        #
        ## Create Anhydrite 03
        container_anhydrite_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_anhydrite_03[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_03
        actual_bottom += thickness_mudstone_05
        #
        ## Create Mudstone 05
        container_mudstone_05 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_05[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_05
        actual_bottom += thickness_anhydrite_02
        #
        ## Create Anhydrite 02
        container_anhydrite_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_anhydrite_02[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_02
        actual_bottom += thickness_mudstone_04
        #
        ## Create Mudstone 04
        container_mudstone_04 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_04[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_04
        actual_bottom += thickness_sandstone_02
        #
        ## Create Sandstone 02
        container_sandstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_02[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_02
        actual_bottom += thickness_mudstone_03
        #
        ## Create Mudstone 03
        container_mudstone_03 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_03[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_03
        actual_bottom += thickness_anhydrite_01
        #
        ## Create Anhydrite 01
        container_anhydrite_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_anhydrite_01[depth] = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                porosity=[0.05, 0.1])
        actual_top += thickness_anhydrite_01
        actual_bottom += thickness_mudstone_02
        #
        ## Create Mudstone 02
        container_mudstone_02 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_02[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_02
        actual_bottom += thickness_sandstone_01
        #
        ## Create Sandstone 01
        container_sandstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_sandstone_01[depth] = SiliciclasticRocks().create_sandstone(number=1, porosity=[0.1, 0.25])
        actual_top += thickness_sandstone_01
        actual_bottom += thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_mudstone_06, container_anhydrite_03, container_mudstone_05, container_anhydrite_02, \
            container_mudstone_04, container_sandstone_02, container_mudstone_03, container_anhydrite_01, \
            container_mudstone_02, container_sandstone_01, container_mudstone_01
    #
    def create_myophorien_series(self, thickness_unit=100, top_unit=0):
        fraction_limestone_01 = round(rd.uniform(0.35, 0.37), 4)
        #
        fraction_mudstone_01 = round(1 - fraction_limestone_01, 4)
        #
        thickness_mudstone_01 = round(thickness_unit*fraction_mudstone_01, 4)
        thickness_limestone_01 = round(thickness_unit*fraction_limestone_01, 4)
        #
        actual_top = top_unit
        actual_bottom = top_unit + thickness_mudstone_01
        #
        ## Create Mudstone 01
        container_mudstone_01 = {}
        steps_rock = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_rock:
            depth = round(i, 4)
            container_mudstone_01[depth] = SiliciclasticRocks().create_mudstone_alt(number=1, porosity=[0.0, 0.1])
        actual_top += thickness_mudstone_01
        actual_bottom += thickness_limestone_01
        #
        ## Create Limestone 01
        container_limestone_01 = {}
        steps_unit = np.linspace(actual_bottom, actual_top, self.resolution, endpoint=False)[::-1]
        for i in steps_unit:
            depth = round(i, 4)
            container_limestone_01[depth] = CarbonateRocks(
                fluid="water", actualThickness=0).create_limestone(number=1, porosity=[0.2, 0.4])
        #
        ## TEST
        # for key, value in reversed(container_mudstone_01.items()):
        #     print(key, value)
        # for key, value in reversed(container_limestone.items()):
        #     print(key, value)
        # for key, value in reversed(container_kupferschiefer.items()):
        #     print(key, value)
        #
        return container_mudstone_01, container_limestone_01
    #