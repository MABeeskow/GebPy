#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		metamorphics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		19.03.2022

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
    def create_granulite(self, number, porosity=None):
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
            #x = round(rd.uniform(0.0, 1.0), 4)
            for mineral in mineral_list:
                if mineral == "Qz":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.25, 0.65), 4)
                        #value = round(0.65 - 0.35*x, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.65:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.60), 4)
                        #value = round(0.10 + 0.50*x, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.60:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.60), 4)
                        #value = round(0.15 - 0.5*x, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.60:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Grt":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                        #value = round(0.05*(1-x), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ky":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.03), 4)
                        #value = round(0.025*(1-x), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.025:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Sil":
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                        #value = round(0.025*(1-x), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.025:
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
            rho_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            youngsmod_helper = 0
            poisson_helper = 0
            vP_helper = 0
            vS_helper = 0
            vPvS_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.05), 4)
            else:
                phi_helper = porosity
            self.data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            self.data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            self.data_garnet_al = Nesosilicates(impurity="pure", data_type=True).create_aluminium_garnet()
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Qz"][n]*self.data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_alkalifeldspar["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Kfs"][n]*self.data_alkalifeldspar["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_plagioclase["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Pl"][n]*self.data_plagioclase["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_garnet_al["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Grt"][n]*self.data_garnet_al["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_kyanite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ky"][n]*self.data_kyanite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_sillimanite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Sil"][n]*self.data_sillimanite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                n_elements += 1
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_alkalifeldspar["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["PE"], 3)
                    elif mineral == "Grt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_garnet_al["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_garnet_al["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_garnet_al["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_garnet_al["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_garnet_al["PE"], 3)
                    elif mineral == "Ky":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_kyanite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_kyanite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_kyanite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_kyanite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_kyanite["PE"], 3)
                    elif mineral == "Sil":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_sillimanite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_sillimanite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_sillimanite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_sillimanite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_sillimanite["PE"], 3)
                #
                rho_helper = round((1-phi_helper)*rho_s_helper + phi_helper*self.data_water[2]/1000, 3)
                youngsmod_helper = round((9*bulkmod_helper*shearmod_helper)/(3*bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round((3*bulkmod_helper - 2*shearmod_helper)/(6*bulkmod_helper + 2*shearmod_helper), 3)
                vP_helper = round(((bulkmod_helper*10**9 + 4/3*shearmod_helper*10**9)/(rho_helper))**0.5, 3)
                vS_helper = round(((shearmod_helper*10**9)/(rho_helper))**0.5, 3)
                vPvS_helper_helper = round(vP_helper/vS_helper, 3)
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
        # for key, value in amounts_mineralogy.items():
        #     print("Mineral:", key, "Mean:", round(np.mean(value)*100, 4), "STD:", round(np.std(value, ddof=1)*100, 4))
        # print("")
        # for key, value in amounts_chemistry.items():
        #     print("Element:", key, "Mean:", round(np.mean(value)*100, 4), "STD:", round(np.std(value, ddof=1)*100, 4))
        # print("")
        # for key, value in bulk_properties.items():
        #     print("Property:", key, "Mean:", round(np.mean(value), 4), "STD:", round(np.std(value, ddof=1), 4))
        #
        results = {}
        results["rock"] = "Granulite"
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

## TEST
# MetamorphicRocks().create_granulite(number=1)
# MetamorphicRocks().create_granulite(number=10)