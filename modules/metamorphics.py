#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		metamorphics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		24.05.2022

# -----------------------------------------------

## MODULES
import numpy as np
import random as rd
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Nesosilicates, Sorosilicates, Phyllosilicates, Inosilicates
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
        self.data_tremolite = Inosilicates(impurity="pure", data_type=True).create_tremolite()
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
    #
    def create_greenschist(self, number, porosity=None):
        #
        self.data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_epidote = Sorosilicates(impurity="pure", data_type=True).create_epidote()
        self.data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
        self.data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
        self.data_antigorite = Phyllosilicates(impurity="pure", data_type=True).create_antigorite()
        #
        assemblage = [self.data_quartz, self.data_plagioclase, self.data_chlorite, self.data_antigorite,
                      self.data_epidote, self.data_muscovite]
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
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.25:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.20, 0.50), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.20 <= value <= 0.50:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.15, 0.30), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.25:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ant":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.05, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ep":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.05, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ms":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1-w_total, 4)
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
            self.data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            self.data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            self.data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Qz"][n]*self.data_quartz["chemistry"][element], 4)
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
                if element in self.data_chlorite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Chl"][n]*self.data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_antigorite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ant"][n]*self.data_antigorite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_muscovite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ms"][n]*self.data_muscovite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_epidote["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ep"][n]*self.data_epidote["chemistry"][element], 4)
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
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["PE"], 3)
                    elif mineral == "Ant":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_antigorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_antigorite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_antigorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_antigorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_antigorite["PE"], 3)
                    elif mineral == "Ep":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_epidote["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["PE"], 3)
                    elif mineral == "Ms":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_muscovite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_muscovite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_muscovite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_muscovite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_muscovite["PE"], 3)
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
        results = {}
        results["rock"] = "Greenschist"
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
    def create_greenschist_basaltic(self, number, porosity=None):
        #
        self.data_actinolite = Inosilicates(impurity="pure", data_type=True).create_actinolite()
        self.data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
        self.data_epidote = Sorosilicates(impurity="pure", data_type=True).create_epidote()
        self.data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        assemblage = [self.data_chlorite, self.data_actinolite, self.data_plagioclase, self.data_epidote]
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
                if mineral == "Act":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.20, 0.30), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.20 <= value <= 0.30:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.25, 0.50), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.50:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ep":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.20 <= value <= 0.55:
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
            x = rd.uniform(0.9, 1.0)
            self.data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(x_value=x)
            self.data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            self.data_actinolite = Inosilicates(impurity="pure", data_type=True).create_actinolite()
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_actinolite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Act"][n]*self.data_actinolite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_chlorite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Chl"][n]*self.data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_epidote["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ep"][n]*self.data_epidote["chemistry"][element], 4)
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
                #
                n_elements += 1
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Act":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_actinolite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_actinolite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_actinolite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_actinolite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_actinolite["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_plagioclase["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["PE"], 3)
                    elif mineral == "Ep":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_epidote["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_epidote["PE"], 3)
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
        results = {}
        results["rock"] = "Greenschist"
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
    def create_greenschist_ultramafic(self, number, porosity=None):
        #
        self.data_brucite = Oxides(impurity="pure", data_type=True).create_brucite()
        self.data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
        self.data_chrysotile = Phyllosilicates(impurity="pure", data_type=True).create_chrysotile()
        self.data_diopside = Inosilicates(impurity="pure", data_type=True).create_diopside()
        self.data_talc = Phyllosilicates(impurity="pure", data_type=True).create_talc()
        self.data_tremolite = Inosilicates(impurity="pure", data_type=True).create_tremolite()
        #
        assemblage = [self.data_chlorite, self.data_chrysotile, self.data_talc, self.data_tremolite,
                      self.data_diopside, self.data_brucite]
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
                if mineral == "Bru":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.25, 0.50), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.50:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ctl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.25, 0.50), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.50:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Tlc":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.15, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.15 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Tr":
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Di":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1-w_total, 4)
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
            self.data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_brucite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Bru"][n]*self.data_brucite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_chlorite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Chl"][n]*self.data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_chrysotile["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ctl"][n]*self.data_chrysotile["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_diopside["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Di"][n]*self.data_diopside["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_talc["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Tlc"][n]*self.data_talc["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_tremolite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Tr"][n]*self.data_tremolite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                #
                n_elements += 1
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Bru":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_brucite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_brucite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_brucite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_brucite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_brucite["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_chlorite["PE"], 3)
                    elif mineral == "Ctl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_chrysotile["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_chrysotile["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_chrysotile["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_chrysotile["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_chrysotile["PE"], 3)
                    elif mineral == "Di":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_diopside["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_diopside["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_diopside["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_diopside["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_diopside["PE"], 3)
                    elif mineral == "Tlc":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_talc["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_talc["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_talc["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_talc["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_talc["PE"], 3)
                    elif mineral == "Tr":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*self.data_tremolite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["PE"], 3)
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
        results = {}
        results["rock"] = "Greenschist"
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
    def create_greenschist_pelitic(self, number, porosity=None):
        #
        data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()                          # fixed
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()  # variable
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()        # variable
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()             # variable
        data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()           # variable
        data_garnet_al = Nesosilicates(impurity="pure", data_type=True).create_aluminium_garnet()      # variable
        data_pyrophyllite = Phyllosilicates(impurity="pure", data_type=True).create_pyrophyllite()     # fixed
        #
        assemblage = [data_quartz, data_alkalifeldspar, data_plagioclase, data_chlorite, data_muscovite,
                      data_garnet_al, data_pyrophyllite]
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
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.70:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.05, 0.15), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.05, 0.15), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.15, 0.30), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.15 <= value <= 0.30:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Ms":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.10), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.10:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Grt":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Prl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1-w_total, 4)
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
            x_kfs = rd.uniform(0.9, 1.0)
            x_pl = rd.uniform(0.9, 1.0)
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar(x_value=x_kfs)
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(x_value=x_pl)
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
            data_garnet_al = Nesosilicates(impurity="pure", data_type=True).create_aluminium_garnet()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in data_quartz["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Qz"][n]*data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_alkalifeldspar["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Kfs"][n]*data_alkalifeldspar["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_plagioclase["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Pl"][n]*data_plagioclase["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_chlorite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Chl"][n]*data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_muscovite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Ms"][n]*data_muscovite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_garnet_al["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Grt"][n]*data_garnet_al["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_pyrophyllite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Prl"][n]*data_pyrophyllite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                #
                n_elements += 1
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_quartz["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_quartz["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_alkalifeldspar["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["PE"], 3)
                    elif mineral == "Ms":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_muscovite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_muscovite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_muscovite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_muscovite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_muscovite["PE"], 3)
                    elif mineral == "Grt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_garnet_al["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_garnet_al["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_garnet_al["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_garnet_al["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_garnet_al["PE"], 3)
                    elif mineral == "Prl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_pyrophyllite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_pyrophyllite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_pyrophyllite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_pyrophyllite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_pyrophyllite["PE"], 3)
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
        results = {}
        results["rock"] = "Greenschist"
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
    def create_amphibolite_ortho(self, number, porosity=None):
        #
        data_actinolite = Inosilicates(impurity="pure", data_type=True).create_actinolite()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
        #
        assemblage = [self.data_tremolite, data_actinolite, data_plagioclase, data_biotite, self.data_quartz,
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
            value_act = 0
            for mineral in mineral_list:
                if mineral == "Qz":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if 0.0 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.15, 0.40), 4)
                    else:
                        value = round(1-w_total, 4)
                    if 0.15 <= value <= 0.40:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Tr":
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                    else:
                        value = round(1-w_total, 4)
                    if 0.0 <= value <= 0.50 - value_act:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.10), 4)
                    else:
                        value = round(1-w_total, 4)
                    if 0.0 <= value <= 0.10:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Act":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.50), 4)
                    else:
                        value = round(1-w_total, 4)
                    if 0.0 <= value <= 0.50:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                        value_act += value
                elif mineral == "Bt":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if 0.0 <= value <= 0.20:
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
            x_pl = rd.uniform(0.9, 1.0)
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(x_value=x_pl)
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            data_actinolite = Inosilicates(impurity="pure", data_type=True).create_actinolite()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Qz"][n]*self.data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_plagioclase["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Pl"][n]*data_plagioclase["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_chlorite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Chl"][n]*data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_actinolite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Act"][n]*data_actinolite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_tremolite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Tr"][n]*self.data_tremolite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_biotite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = round(amounts_mineralogy["Bt"][n]*data_biotite["chemistry"][element], 4)
                    else:
                        value = round(1-w_total, 4)
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
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["K"], 3)
                        shearmod_helper += round(shear_factor*amounts_mineralogy[mineral][n]*self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_quartz["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["K"], 3)
                        shearmod_helper += round(shear_factor*amounts_mineralogy[mineral][n]*data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_plagioclase["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["K"], 3)
                        shearmod_helper += round(shear_factor*amounts_mineralogy[mineral][n]*data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_chlorite["PE"], 3)
                    elif mineral == "Bt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_biotite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_biotite["K"], 3)
                        shearmod_helper += round(shear_factor*amounts_mineralogy[mineral][n]*data_biotite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_biotite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_biotite["PE"], 3)
                    elif mineral == "Act":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_actinolite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_actinolite["K"], 3)
                        shearmod_helper += round(shear_factor*amounts_mineralogy[mineral][n]*data_actinolite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_actinolite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_actinolite["PE"], 3)
                    elif mineral == "Tr":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["K"], 3)
                        shearmod_helper += round(shear_factor*amounts_mineralogy[mineral][n]*self.data_tremolite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*self.data_tremolite["PE"], 3)
                #
                rho_helper = round((1-phi_helper)*rho_s_helper + phi_helper*self.data_water[2]/1000, 3)
                youngsmod_helper = round((9*bulkmod_helper*shearmod_helper)/(3*bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round((3*bulkmod_helper - 2*shearmod_helper)/(6*bulkmod_helper + 2*shearmod_helper), 4)
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
        results = {}
        results["rock"] = "Ortho-Amphibolite"
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
# print("Test: Granulite")
# print(MetamorphicRocks().create_granulite(number=1))
# print(MetamorphicRocks().create_granulite(number=10)
# print("Test: Greenschist")
# print(MetamorphicRocks().create_greenschist(number=1))
# print(MetamorphicRocks().create_greenschist(number=10))
# print("Test: Amphibolite")
# print(MetamorphicRocks().create_amphibolite_ortho(number=1))