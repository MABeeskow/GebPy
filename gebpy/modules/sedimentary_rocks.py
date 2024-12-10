#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sedimentary_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		16.05.2023

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
    def create_marl_alt(self, rock="Marl", number=1, composition=None, classification="Marl", porosity=None):
        #
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["phi"] = []
        results_container["fluid"] = self.fluid
        results_container["rho_s"] = []
        results_container["rho"] = []
        results_container["vP"] = []
        results_container["vS"] = []
        results_container["vP/vS"] = []
        results_container["K"] = []
        results_container["G"] = []
        results_container["E"] = []
        results_container["nu"] = []
        results_container["GR"] = []
        results_container["PE"] = []
        #
        n = 0
        while n < number:
            data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            #
            mineralogy = {"Cal": self.data_calcite, "Ilt": data_illite, "Pl": data_plagioclase}
            #
            minerals_list = list(mineralogy.keys())
            #
            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []
            #
            condition = False
            #
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}
                #
                if composition != None:
                    phi_cal = composition["Cal"]
                    phi_ilt = composition["Ilt"]
                    phi_pl = composition["Pl"]
                    #
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Pl"] = phi_pl
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        if classification == "Marl":
                            cal_limits = [0.35, 0.65]
                            ilt_limits = [0.35, 0.65]
                            pl_limits = [0.0, 0.3]
                        #
                        phi_cal = round(rd.uniform(cal_limits[0], cal_limits[1]), 4)
                        phi_ilt = round(rd.uniform(ilt_limits[0], (1 - phi_cal)), 4)
                        phi_pl = round(1 - phi_cal - phi_ilt, 4)
                        #
                        phi_total = phi_cal + phi_ilt + phi_pl
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if cal_limits[0] <= phi_cal <= cal_limits[1] \
                                    and ilt_limits[0] <= phi_ilt <= ilt_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1]:
                                condition_2 = True
                        #
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Pl"] = phi_pl
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                    #
                    for element, value in mineralogy[key]["chemistry"].items():
                        if element not in elements_list:
                            elements_list.append(element)
                            w_elements[element] = 0.0
                #
                if elements_list[0] not in results_container["chemistry"]:
                    for element in elements_list:
                        results_container["chemistry"][element] = []
                #
                rho_s = round(rho_s, 3)
                for key, value in phi_minerals.items():
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                elif type(porosity) == float:
                    var_porosity = porosity
                else:
                    var_porosity = round(rd.uniform(porosity[0], porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                #
                old_index = elements_list.index("O")
                elements_list += [elements_list.pop(old_index)]
                #
                w_elements_total = 0.0
                for element in elements_list:
                    if element != "O":
                        for mineral, w_mineral in w_minerals.items():
                            if element in mineralogy[mineral]["chemistry"]:
                                if element == "U":
                                    n_digits = 4
                                else:
                                    n_digits = 4
                                #
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        #
                        w_elements[element] = round(w_elements[element], 4)
                #
                total_w_minerals = round(sum(w_minerals.values()), 4)
                total_w_elements = round(sum(w_elements.values()), 4)
                if total_w_minerals == 1.0 and total_w_elements == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    #
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    #
                    condition = True
            #
            rho_solid = 0
            velocity_solid = {"vP": 0, "vS": 0}
            gamma_ray = 0.0
            photoelectricity = 0.0
            for key, value in phi_minerals.items():
                rho_solid += phi_minerals[key]*mineralogy[key]["rho"]
                velocity_solid["vP"] += phi_minerals[key]*mineralogy[key]["vP"]
                velocity_solid["vS"] += phi_minerals[key]*mineralogy[key]["vS"]
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
            #
            ## Density
            rho_solid = round(rho_solid, 3)
            rho = round((1 - var_porosity)*rho_solid + var_porosity*data_fluid[2]/1000, 3)
            ## Seismic Velocities
            vP_factor = 4/3
            vP = round(velocity_solid["vP"]*(1 - vP_factor*var_porosity), 3)
            vS_factor = 3/2
            vS = round(velocity_solid["vS"]*(1 - vS_factor*var_porosity), 3)
            vPvS = round(vP/vS, 6)
            ## Elastic Parameters
            bulk_modulus = round(rho*(vP**2 - 4/3*vS**2)*10**(-9), 3)
            shear_modulus = round((rho*vS**2)*10**(-9), 3)
            youngs_modulus = round((9*bulk_modulus*shear_modulus)/(3*bulk_modulus + shear_modulus), 3)
            poisson_ratio = round((3*bulk_modulus - 2*shear_modulus)/(6*bulk_modulus + 2*shear_modulus), 4)
            ## Gamma Ray
            gamma_ray = round(gamma_ray, 3)
            ## Photoelectricity
            photoelectricity = round(photoelectricity, 3)
            #
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(float(value))
            #
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(float(value))
            #
            results_container["phi"].append(float(var_porosity))
            results_container["rho_s"].append(float(rho_s))
            results_container["rho"].append(float(rho))
            results_container["vP"].append(float(vP))
            results_container["vS"].append(float(vS))
            results_container["vP/vS"].append(float(vPvS))
            results_container["K"].append(float(bulk_modulus))
            results_container["G"].append(float(shear_modulus))
            results_container["E"].append(float(youngs_modulus))
            results_container["nu"].append(float(poisson_ratio))
            results_container["GR"].append(float(gamma_ray))
            results_container["PE"].append(float(photoelectricity))
            #
            n += 1
        #
        return results_container
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
                            elif clay_mineral == "Ilt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, 1.0), 4)
                                else:
                                    value = round(value_fract*(1 - w_total_clay), 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
                            elif clay_mineral == "Kln":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, 1.0), 4)
                                else:
                                    value = round(value_fract*(1 - w_total_clay), 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
                            elif clay_mineral == "Mnt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract - w_total_clay, 4)
                                else:
                                    value = round(value_fract - w_total_clay, 4)
                                amounts_helper_minerals[clay_mineral] = value
                                amounts_helper.append(value)
                                w_total_clay += value
                                n_minerals += 1
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
                            elif accessory_mineral == "Pl":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract*rd.uniform(0.0, value_fract), 4)
                                else:
                                    value = round(value_fract - w_total_acc, 4)
                                amounts_helper_minerals[accessory_mineral] = value
                                amounts_helper.append(value)
                                w_total_acc += value
                                n_minerals += 1
                            elif accessory_mineral == "Qz":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(value_fract - w_total_acc, 4)
                                else:
                                    value = round(value_fract - w_total_acc, 4)
                                amounts_helper_minerals[accessory_mineral] = value
                                amounts_helper.append(value)
                                w_total_acc += value
                                n_minerals += 1
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