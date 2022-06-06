#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sedimentary_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		06.06.2022

#-----------------------------------------------

## MODULES
import datetime
import sys

import numpy as np
from numpy import round
import random as rd
from modules.oxides import Oxides
from modules.carbonates import Carbonates
from modules.silicates import Phyllosilicates
from modules.silicates import Tectosilicates
from modules.fluids import Water
import time

#####################
## SANDSTONE ROCKS ##
#####################
class Sandstone:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        self.data_water = Water.water("")
    #
    def create_sandstone(self, number, porosity=None):
        #
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()   # variable
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()         # variable
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()              # variable
        #
        assemblage = [self.data_quartz, data_alkalifeldspar, data_plagioclase,  self.data_calcite, self.data_hematite,
                      data_chlorite]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
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
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        amounts_helper = []
        while n < number:
            w_total = 0
            n_minerals = 0
            for mineral in mineral_list:
                if mineral == "Qz":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(1 - w_total, 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.85 <= value <= 1.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Cal":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
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
        n = 0
        amounts_helper = {}
        while n < number:
            w_total = 0
            n_elements = 0
            rho_s_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.05), 4)
            else:
                phi_helper = porosity
            #
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Qz"][n] * self.data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_alkalifeldspar["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Kfs"][n] * data_alkalifeldspar["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_plagioclase["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Pl"][n] * data_plagioclase["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_chlorite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Chl"][n] * data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_hematite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Hem"][n] * self.data_hematite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                #
                n_elements += 1
            #
            shear_factor = 1.0
            #
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_alkalifeldspar["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["PE"], 3)
                    elif mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["PE"], 3)
                #
                rho_helper = round((1 - phi_helper) * rho_s_helper + phi_helper * self.data_water[2] / 1000, 3)
                youngsmod_helper = round(
                    (9 * bulkmod_helper * shearmod_helper) / (3 * bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round(
                    (3 * bulkmod_helper - 2 * shearmod_helper) / (6 * bulkmod_helper + 2 * shearmod_helper), 3)
                vP_helper = round(
                    ((bulkmod_helper * 10 ** 9 + 4 / 3 * shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vS_helper = round(((shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vPvS_helper_helper = round(vP_helper / vS_helper, 3)
                #
                bulk_properties["rho_s"].append(round(rho_s_helper, 3))
                bulk_properties["rho"].append(rho_helper)
                bulk_properties["K"].append(round(bulkmod_helper, 3))
                bulk_properties["G"].append(round(shearmod_helper, 3))
                bulk_properties["E"].append(youngsmod_helper)
                bulk_properties["nu"].append(poisson_helper)
                bulk_properties["vP"].append(vP_helper)
                bulk_properties["vS"].append(vS_helper)
                bulk_properties["vPvS"].append(vPvS_helper_helper)
                bulk_properties["GR"].append(round(gr_helper, 3))
                bulk_properties["PE"].append(round(pe_helper, 3))
                bulk_properties["phi"].append(round(phi_helper, 3))
                n += 1
        #
        results = {}
        results["rock"] = "Sandstone"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = bulk_properties["phi"][0]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"][0]
            results["rho"] = bulk_properties["rho"][0]
            results["vP"] = bulk_properties["vP"][0]
            results["vS"] = bulk_properties["vS"][0]
            results["vP/vS"] = bulk_properties["vPvS"][0]
            results["K"] = bulk_properties["K"][0]
            results["G"] = bulk_properties["G"][0]
            results["E"] = bulk_properties["E"][0]
            results["nu"] = bulk_properties["nu"][0]
            results["GR"] = bulk_properties["GR"][0]
            results["PE"] = bulk_properties["PE"][0]
        #
        return results
    #
#############################
## OTHER SEDIMENTARY ROCKS ##
#############################
class SedimentaryRocks:
    #
    def __init__(self, fluid="water", actualThickness=100):
        # Settings
        self.fluid = fluid
        self.actualThickness = actualThickness
        # Oxide Minerals
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        # Carbonate Minerals
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_dolomite = Carbonates(impurity="pure", data_type=True).create_dolomite()
        # Phyllosilicate Minerals
        self.data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        # Fluids
        self.data_water = Water.water("")
    #
    def create_marl(self, number, porosity=None):
        #
        data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()                     # variable
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()            # variable
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()                 # variable
        data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()   # variable
        #
        assemblage = [self.data_calcite, self.data_dolomite, data_illite, data_chlorite, self.data_kaolinite,
                      data_montmorillonite, self.data_quartz, data_plagioclase]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        amounts_helper_minerals = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage:
            amounts_mineralogy[mineral["mineral"]] = []
            mineral_list.append(mineral["mineral"])
            elements_mineral = list(mineral["chemistry"].keys())
            amounts_helper_minerals[mineral["mineral"]] = 0
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        amounts_helper = []
        amounts_fractions_helper = {}
        fraction_list = ["carbonate", "clay", "accessory"]
        fraction_list.sort()
        for fraction in fraction_list:
            amounts_fractions_helper[fraction] = 0
        carbonate_list = ["Cal"]
        clay_list = ["Ilt", "Chl", "Kln", "Mnt"]
        accessory_list = ["Dol", "Pl", "Qz"]
        while n < number:
            w_total_acc = 0
            w_total_clay = 0
            w_total_fractions = 0
            n_minerals = 0
            n_fractions = 0
            for fraction in fraction_list:
                if fraction == "carbonate":
                    if n_fractions < len(fraction_list) - 1:
                        value_fract = round(rd.uniform(0.35, 0.65), 4)
                    else:
                        value_fract = round(1 - w_total_fractions, 4)
                    if 0.35 <= value_fract <= 0.65:
                        amounts_helper_minerals["Cal"] = value_fract
                        #
                        amounts_fractions_helper[fraction] = value_fract
                        w_total_fractions += value_fract
                        n_fractions += 1
                        n_minerals += 1
                elif fraction == "clay":
                    if n_fractions < len(fraction_list) - 1:
                        value_fract = round(1 - w_total_fractions, 4)
                    else:
                        value_fract = round(1 - w_total_fractions, 4)
                    if 0.35 <= value_fract <= 0.65:
                        for clay_mineral in clay_list:
                            if clay_mineral == "Chl":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, 1.0), 4)
                                else:
                                    value = round(value_fract(1 - w_total_clay), 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
                                # if 0.0 <= value <= value_fract:
                                #     amounts_helper_minerals[clay_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_clay += value
                                #     n_minerals += 1
                            elif clay_mineral == "Ilt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, 1.0), 4)
                                else:
                                    value = round(value_fract*(1 - w_total_clay), 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
                                # if 0.50 <= value <= value_fract:
                                #     amounts_helper_minerals[clay_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_clay += value
                                #     n_minerals += 1
                            elif clay_mineral == "Kln":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, 1.0), 4)
                                else:
                                    value = round(value_fract*(1 - w_total_clay), 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
                                # if 0.0 <= value <= value_fract:
                                #     amounts_helper_minerals[clay_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_clay += value
                                #     n_minerals += 1
                            elif clay_mineral == "Mnt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract - w_total_clay, 4)
                                else:
                                    value = round(value_fract - w_total_clay, 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
                                # if 0.0 <= value <= value_fract:
                                #     amounts_helper_minerals[clay_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_clay += value
                                #     n_minerals += 1
                        #
                        amounts_fractions_helper[fraction] = value_fract
                        w_total_fractions += value_fract
                        n_fractions += 1
                elif fraction == "accessory":
                    if n_fractions < len(fraction_list) - 1:
                        value_fract = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value_fract = round(1 - w_total_fractions, 4)
                    if 0.0 <= value_fract <= 0.05:
                        for accessory_mineral in accessory_list:
                            if accessory_mineral == "Dol":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, value_fract), 4)
                                else:
                                    value = round(value_fract - w_total_acc, 4)
                                amounts_helper_minerals[accessory_mineral] = value
                                amounts_helper.append(value)
                                w_total_acc += value
                                n_minerals += 1
                                # if 0.0 <= value <= value_fract:
                                #     amounts_helper_minerals[accessory_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_acc += value
                                #     n_minerals += 1
                            elif accessory_mineral == "Pl":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, value_fract), 4)
                                else:
                                    value = round(value_fract - w_total_acc, 4)
                                amounts_helper_minerals[accessory_mineral] = value
                                amounts_helper.append(value)
                                w_total_acc += value
                                n_minerals += 1
                                # if 0.0 <= value <= value_fract:
                                #     amounts_helper_minerals[accessory_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_acc += value
                                #     n_minerals += 1
                            elif accessory_mineral == "Qz":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract - w_total_acc, 4)
                                else:
                                    value = round(value_fract - w_total_acc, 4)
                                amounts_helper_minerals[accessory_mineral] = value
                                amounts_helper.append(value)
                                w_total_acc += value
                                n_minerals += 1
                                # if 0.0 <= value <= value_fract:
                                #     amounts_helper_minerals[accessory_mineral] = value
                                #     amounts_helper.append(value)
                                #     w_total_acc += value
                                #     n_minerals += 1
                        #
                        amounts_fractions_helper[fraction] = value_fract
                        w_total_fractions += value_fract
                        n_fractions += 1
            #
            # print("Fractions:", amounts_fractions_helper)
            # print(np.sum([*amounts_fractions_helper.values()]))
            # print("Mineralogy:", amounts_helper_minerals)
            # print(np.sum([*amounts_helper_minerals.values()]))
            #
            if np.sum([*amounts_helper_minerals.values()]) == 1.0 and n_minerals == len(mineral_list):
                for index, mineral in enumerate(mineral_list):
                    amounts_mineralogy[mineral].append(amounts_helper_minerals[mineral])
                n += 1
                amounts_helper_minerals.clear()
            else:
                n += 0
                for mineral in amounts_helper_minerals:
                    amounts_helper_minerals[mineral] = 0
                for fraction in amounts_fractions_helper:
                    amounts_fractions_helper[fraction] = 0
        #
        n = 0
        amounts_helper = {}
        while n < number:
            w_total = 0
            n_elements = 0
            rho_s_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.25), 4)
            else:
                phi_helper = porosity
            #
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
            data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_chlorite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Chl"][n] * data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_dolomite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Dol"][n] * self.data_dolomite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_illite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Ilt"][n] * data_illite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_kaolinite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Kln"][n] * self.data_kaolinite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_montmorillonite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Mnt"][n] * data_montmorillonite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_plagioclase["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Pl"][n] * data_plagioclase["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Qz"][n] * self.data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                #
                n_elements += 1
            #
            bulk_factor = 0.15
            shear_factor = 0.05
            #
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * data_chlorite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["PE"], 3)
                    elif mineral == "Dol":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * self.data_dolomite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_dolomite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["PE"], 3)
                    elif mineral == "Ilt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_illite["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * data_illite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_illite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_illite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_illite["PE"], 3)
                    elif mineral == "Kln":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * self.data_kaolinite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_kaolinite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["PE"], 3)
                    elif mineral == "Mnt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * data_montmorillonite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_montmorillonite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * data_plagioclase["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["PE"], 3)
                    elif mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(bulk_factor*amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                #
                rho_helper = round((1 - phi_helper) * rho_s_helper + phi_helper * self.data_water[2] / 1000, 3)
                youngsmod_helper = round(
                    (9 * bulkmod_helper * shearmod_helper) / (3 * bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round(
                    (3 * bulkmod_helper - 2 * shearmod_helper) / (6 * bulkmod_helper + 2 * shearmod_helper), 3)
                vP_helper = round(
                    ((bulkmod_helper * 10 ** 9 + 4 / 3 * shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vS_helper = round(((shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vPvS_helper_helper = round(vP_helper / vS_helper, 3)
                #
                bulk_properties["rho_s"].append(round(rho_s_helper, 3))
                bulk_properties["rho"].append(rho_helper)
                bulk_properties["K"].append(round(bulkmod_helper, 3))
                bulk_properties["G"].append(round(shearmod_helper, 3))
                bulk_properties["E"].append(youngsmod_helper)
                bulk_properties["nu"].append(poisson_helper)
                bulk_properties["vP"].append(vP_helper)
                bulk_properties["vS"].append(vS_helper)
                bulk_properties["vPvS"].append(vPvS_helper_helper)
                bulk_properties["GR"].append(round(gr_helper, 3))
                bulk_properties["PE"].append(round(pe_helper, 3))
                bulk_properties["phi"].append(round(phi_helper, 3))
                n += 1
        #
        results = {}
        results["rock"] = "Marl"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = bulk_properties["phi"][0]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"][0]
            results["rho"] = bulk_properties["rho"][0]
            results["vP"] = bulk_properties["vP"][0]
            results["vS"] = bulk_properties["vS"][0]
            results["vP/vS"] = bulk_properties["vPvS"][0]
            results["K"] = bulk_properties["K"][0]
            results["G"] = bulk_properties["G"][0]
            results["E"] = bulk_properties["E"][0]
            results["nu"] = bulk_properties["nu"][0]
            results["GR"] = bulk_properties["GR"][0]
            results["PE"] = bulk_properties["PE"][0]
        #
        return results
    #
## TEST
# print("Test: Marl")
# start = time.process_time()
# print(SedimentaryRocks(fluid="water", actualThickness=0).create_marl(number=1))
# end = time.process_time()
# print("Time:", round(end - start, 3))
# start = time.process_time()
# print(SedimentaryRocks(fluid="water", actualThickness=0).create_marl(number=10))
# end = time.process_time()
# print("Time:", round(end - start, 3))
# start = time.process_time()
# print(SedimentaryRocks(fluid="water", actualThickness=0).create_marl(number=100))
# end = time.process_time()
# print("Time:", round(end - start, 3))