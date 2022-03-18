#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		metamorphics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.03.2022

# -----------------------------------------------

## MODULES
import numpy as np
import random as rd
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Nesosilicates
from modules.fluids import Water

#######################
## METAMORPHIC ROCKS ##
#######################
class MetamorphicRocks:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_kyanite = Nesosilicates(impurity="pure", data_type=True).create_kyanite()
        self.data_sillimanite = Nesosilicates(impurity="pure", data_type=True).create_sillimanite()
        self.data_water = Water.water("")
    #
    def create_granulite(self, number):
        #
        self.data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        self.data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        self.data_garnet_al = Nesosilicates(impurity="pure", data_type=True).create_aluminium_garnet()
        #
        assemblage = [self.data_quartz, self.data_alkalifeldspar, self.data_plagioclase, self.data_garnet_al,
                      self.data_kyanite, self.data_sillimanite]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        mineral_list = []
        elements = []
        for mineral in assemblage:
            amounts_mineralogy[mineral["mineral"]] = []
            mineral_list.append(mineral["mineral"])
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        elements.sort()
        #
        n = 0
        amounts_helper = []
        while n < number:
            w_total = 0
            n_minerals = 0
            for mineral in mineral_list:
                if mineral == "Qz":
                    value = round(rd.uniform(0.35, 0.65), 4)
                    if value >= 0.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    value = round(rd.uniform(0.15, 0.60), 4)
                    if value >= 0.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    value = round(rd.uniform(0.05, 0.15), 4)
                    if value >= 0.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Grt":
                    value = round(rd.uniform(0.0, 0.05), 4)
                    if value >= 0.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ky":
                    value = round(rd.uniform(0.0, 0.03), 4)
                    if value >= 0.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Sil":
                    value = round(1-w_total, 4)
                    if value >= 0.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
            #
            if np.sum(amounts_helper) == 1.0 and n_minerals == len(mineral_list):
                for index, mineral in enumerate(mineral_list):
                    amounts_mineralogy[mineral].append(amounts_helper[index])
                n += 1
                amounts_helper.clear()
            else:
                n += 0
                amounts_helper.clear()
        #
        for key, value in amounts_mineralogy.items():
            print(key, value)
        for key, value in amounts_chemistry.items():
            print(key, value)

## TEST
MetamorphicRocks().create_granulite(number=10)