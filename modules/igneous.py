#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		igneous.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		10.12.2024

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules.chemistry import PeriodicSystem, OxideCompounds
from modules import minerals, geochemistry, oxides
from modules import fluids
from modules.geophysics import Elasticity as elast
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Phyllosilicates, Inosilicates, Nesosilicates
from modules.fluids import Water, Hydrocarbons
from modules.petrophysics import SeismicVelocities

## CLASSES
class UltraMafic:
    #
    def __init__(self, fluid, actualThickness, dict_output=False, porosity=None):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.dict_output = dict_output
        self.porosity = porosity
        #
        self.data_qz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_water = Water.water("")
        self.data_oil = Hydrocarbons.oil("")
        self.data_gas = Hydrocarbons.natural_gas("")
    #
    def create_ultramafic_rock(self, rock="Granite", number=1, composition=None):
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
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
            data_clinopyroxene = Inosilicates(
                mineral="Clinopyroxene", data_type=True, traces_list=[]).generate_dataset(
                number=1)
            data_orthopyroxene = Inosilicates(
                mineral="Orthopyroxene", data_type=True, traces_list=[]).generate_dataset(
                number=1)
            data_olivine = Nesosilicates(
                mineral="Olivine", data_type=True, traces_list=[]).generate_dataset(
                number=1)
            #
            mineralogy = {"Opx": data_orthopyroxene, "Cpx": data_clinopyroxene, "Ol": data_olivine}
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
                    phi_opx = composition["Opx"]
                    phi_cpx = composition["Cpx"]
                    phi_ol = composition["Ol"]
                    #
                    phi_minerals["Opx"] = phi_opx
                    phi_minerals["Cpx"] = phi_cpx
                    phi_minerals["Ol"] = phi_ol
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        if rock == "Orthopyroxenite":
                            opx_limits = [0.9, 1.0]
                            cpx_limits = [0.0, 0.1]
                            ol_limits = [0.0, 0.1]
                        elif rock == "Clinopyroxenite":
                            opx_limits = [0.0, 0.1]
                            cpx_limits = [0.9, 1.0]
                            ol_limits = [0.0, 0.1]
                        elif rock == "Dunite":
                            opx_limits = [0.0, 0.1]
                            cpx_limits = [0.0, 0.1]
                            ol_limits = [0.9, 1.0]
                        elif rock == "Harzburgite":
                            opx_limits = [0.05, 0.6]
                            cpx_limits = [0.0, 0.1]
                            ol_limits = [0.4, 0.9]
                        elif rock == "Wehrlite":
                            opx_limits = [0.0, 0.1]
                            cpx_limits = [0.05, 0.6]
                            ol_limits = [0.4, 0.9]
                        elif rock == "Websterite":
                            opx_limits = [0.05, 0.9]
                            cpx_limits = [0.05, 0.9]
                            ol_limits = [0.0, 0.05]
                        elif rock == "Lherzolite":
                            opx_limits = [0.05, 0.55]
                            cpx_limits = [0.05, 0.55]
                            ol_limits = [0.4, 0.9]
                        elif rock == "Olivine-Websterite":
                            opx_limits = [0.05, 0.9]
                            cpx_limits = [0.05, 0.9]
                            ol_limits = [0.05, 0.4]
                        elif rock == "Olivine-Orthopyroxenite":
                            opx_limits = [0.55, 0.9]
                            cpx_limits = [0.0, 0.05]
                            ol_limits = [0.05, 0.4]
                        elif rock == "Olivine-Clinopyroxenite":
                            opx_limits = [0.0, 0.05]
                            cpx_limits = [0.55, 0.9]
                            ol_limits = [0.05, 0.4]
                        elif rock == "Peridotite":
                            opx_limits = [0.0, 0.6]
                            cpx_limits = [0.0, 0.6]
                            ol_limits = [0.4, 1.0]
                        elif rock == "Pyroxenite":
                            opx_limits = [0.0, 1.0]
                            cpx_limits = [0.0, 1.0]
                            ol_limits = [0.0, 0.4]
                        #
                        phi_opx = round(rd.uniform(opx_limits[0], opx_limits[1]), 4)
                        phi_cpx = round(rd.uniform(cpx_limits[0], (1.0 - phi_opx)), 4)
                        phi_ol = round(1 - phi_opx - phi_cpx, 4)
                        phi_total = phi_opx + phi_cpx + phi_ol
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if opx_limits[0] <= phi_opx <= opx_limits[1] \
                                    and cpx_limits[0] <= phi_cpx <= cpx_limits[1] \
                                    and ol_limits[0] <= phi_ol <= ol_limits[1]:
                                condition_2 = True
                        #
                    phi_minerals["Opx"] = phi_opx
                    phi_minerals["Cpx"] = phi_cpx
                    phi_minerals["Ol"] = phi_ol
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"][0]
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"][0])/rho_s, n_digits)
                #
                if self.porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                else:
                    var_porosity = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
                #
                rho = round((1 - var_porosity)*rho_s + var_porosity*data_fluid[2]/1000, 3)
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
                                    n_digits = 6
                                else:
                                    n_digits = 6
                                #
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element][0], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 6)
                        #
                        w_elements[element] = round(w_elements[element], 6)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    #
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    #
                    condition = True
            #
            bulk_mod = 0.0
            shear_mod = 0.0
            gamma_ray = 0.0
            photoelectricity = 0.0
            #
            K_list = []
            G_list = []
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"][0]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"][0]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"][0], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"][0], 3))
                phi_list.append(phi_minerals[key])
            #
            K_geo = elast.calc_geometric_mean(self, phi_list, K_list)
            G_geo = elast.calc_geometric_mean(self, phi_list, G_list)
            #
            anisotropic_factor = 1.0
            #
            bulk_mod = K_geo/anisotropic_factor
            shear_mod = G_geo/anisotropic_factor
            #
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 6)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 6)
            # Composition data
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            amounts = []
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_elements = list(w_elements.keys())
            list_oxides = OxideCompounds(var_list_elements=list_elements).find_oxides()
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)

            if list_oxides[0] not in results_container["compounds"]:
                for oxide in list_oxides:
                    results_container["compounds"][oxide] = []

            for key, value in composition_oxides.items():
                results_container["compounds"][key].append(value)

            results_container["mineralogy"] = dict(sorted(
                results_container["mineralogy"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
            results_container["chemistry"] = dict(sorted(
                results_container["chemistry"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
            results_container["compounds"] = dict(sorted(
                results_container["compounds"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))

            # Results
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_s)
            results_container["rho"].append(rho)
            results_container["vP"].append(vP)
            results_container["vS"].append(vS)
            results_container["vP/vS"].append(vPvS)
            results_container["K"].append(bulk_mod)
            results_container["G"].append(shear_mod)
            results_container["E"].append(youngs_mod)
            results_container["nu"].append(poisson_rat)
            results_container["GR"].append(gamma_ray)
            results_container["PE"].append(photoelectricity)
            #
            n += 1
        #
        return results_container
    #
class Plutonic:
    #
    def __init__(self, fluid, actualThickness, dict_output=False, porosity=None):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.dict_output = dict_output
        self.porosity = porosity
        #
        self.data_qz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_water = Water.water("")
        self.data_oil = Hydrocarbons.oil("")
        self.data_gas = Hydrocarbons.natural_gas("")
    #
    def create_granite_streckeisen_alt(self, number=1, composition=None):
        results_container = {}
        results_container["rock"] = "Granite"
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["phi"] = []
        results_container["fluid"] = "water"
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
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_bt = composition["Bt"]
                    #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        phi_qz = round(rd.uniform(0.2, 0.6), 4)
                        phi_kfs = round(rd.uniform(0.15, (1.0 - phi_qz)), 4)
                        phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                        phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                        phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.2 <= phi_qz <= 0.6 and 0.15 <= phi_kfs <= 0.8 and 0.0 <= phi_pl <= 0.52 \
                                    and 0.0 <= phi_bt <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if self.porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                else:
                    var_porosity = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
                #
                rho = round((1 - var_porosity)*rho_s + var_porosity*data_fluid[2]/1000, 3)
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
                                    n_digits = 6
                                else:
                                    n_digits = 6
                                #
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 6)
                        #
                        w_elements[element] = round(w_elements[element], 6)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    condition = True
            #
            bulk_mod = 0.0
            shear_mod = 0.0
            gamma_ray = 0.0
            photoelectricity = 0.0
            #
            K_list = []
            G_list = []
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            K_geo = elast.calc_geometric_mean(self, phi_list, K_list)
            G_geo = elast.calc_geometric_mean(self, phi_list, G_list)
            #
            # anisotropic_factor = round(rd.uniform(1.0, 2.0), 4)
            anisotropic_factor = 1.0
            #
            bulk_mod = K_geo/anisotropic_factor
            shear_mod = G_geo/anisotropic_factor
            #
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 6)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 6)
            #
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)
            #
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
            #
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_s)
            results_container["rho"].append(rho)
            results_container["vP"].append(vP)
            results_container["vS"].append(vS)
            results_container["vP/vS"].append(vPvS)
            results_container["K"].append(bulk_mod)
            results_container["G"].append(shear_mod)
            results_container["E"].append(youngs_mod)
            results_container["nu"].append(poisson_rat)
            results_container["GR"].append(gamma_ray)
            results_container["PE"].append(photoelectricity)
            #
            n += 1
        #
        return results_container
    #
    def create_granodiorite_streckeisen_alt(self, number=1, composition=None):
        results_container = {}
        results_container["rock"] = "Granodiorite"
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["phi"] = []
        results_container["fluid"] = "water"
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
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_bt = composition["Bt"]
                    #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        phi_qz = round(rd.uniform(0.2, 0.6), 4)
                        phi_kfs = round(rd.uniform(0.05, (1.0 - phi_qz)), 4)
                        phi_pl = round(rd.uniform(0.25 , (1.0 - phi_qz - phi_kfs)), 4)
                        phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                        phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.2 <= phi_qz <= 0.6 and 0.05 <= phi_kfs <= 0.28 and 0.25 <= phi_pl <= 0.72 \
                                    and 0.0 <= phi_bt <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if self.porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                else:
                    var_porosity = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
                #
                rho = round((1 - var_porosity)*rho_s + var_porosity*data_fluid[2]/1000, 3)
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
                                    n_digits = 6
                                else:
                                    n_digits = 6
                                #
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 6)
                        #
                        w_elements[element] = round(w_elements[element], 6)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    condition = True
            #
            bulk_mod = 0.0
            shear_mod = 0.0
            gamma_ray = 0.0
            photoelectricity = 0.0
            #
            K_list = []
            G_list = []
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            K_geo = elast.calc_geometric_mean(self, phi_list, K_list)
            G_geo = elast.calc_geometric_mean(self, phi_list, G_list)
            #
            # anisotropic_factor = round(rd.uniform(1.0, 2.0), 4)
            anisotropic_factor = 1.0
            #
            bulk_mod = K_geo/anisotropic_factor
            shear_mod = G_geo/anisotropic_factor
            #
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 6)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 6)
            #
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)
            #
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
            #
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_s)
            results_container["rho"].append(rho)
            results_container["vP"].append(vP)
            results_container["vS"].append(vS)
            results_container["vP/vS"].append(vPvS)
            results_container["K"].append(bulk_mod)
            results_container["G"].append(shear_mod)
            results_container["E"].append(youngs_mod)
            results_container["nu"].append(poisson_rat)
            results_container["GR"].append(gamma_ray)
            results_container["PE"].append(photoelectricity)
            #
            n += 1
        #
        return results_container
    #
    def create_plutonic_rock_streckeisen(self, rock="Granite", number=1, composition=None, enrichment_kfs=None,
                                         enrichment_pl=None, upper_streckeisen=True, porosity=None):
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
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
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar(
                enrichment=enrichment_kfs)
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(
                enrichment=enrichment_pl)
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            if upper_streckeisen == True:
                mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                              "Bt": data_biotite}
            elif upper_streckeisen == False:
                data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
                mineralogy = {"Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                              "Bt": data_biotite}
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
                    if upper_streckeisen == True:
                        phi_qz = composition["Qz"]
                        phi_minerals["Qz"] = phi_qz
                    elif upper_streckeisen == False:
                        phi_nph = composition["Nph"]
                        phi_minerals["Nph"] = phi_nph
                    #
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_bt = composition["Bt"]
                    #
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        #
                        ## UPPER STRECKEISEN DIAGRAM
                        #
                        if rock == "Granite":
                            qz_limits = [0.2, 0.6]
                            kfs_limits = [0.15, 0.8]
                            pl_limits = [0.0, 0.52]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Granodiorite":
                            qz_limits = [0.2, 0.6]
                            kfs_limits = [0.05, 0.28]
                            pl_limits = [0.25, 0.72]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Tonalite":
                            qz_limits = [0.2, 0.6]
                            kfs_limits = [0.0, 0.08]
                            pl_limits = [0.35, 0.8]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Monzodiorite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.52, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Monzogabbro":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.08, 0.35]
                            pl_limits = [0.52, 0.9]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Gabbro":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.1]
                            pl_limits = [0.72, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Diorite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.1]
                            pl_limits = [0.72, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Monzonite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.28, 0.65]
                            pl_limits = [0.28, 0.65]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Syenite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.52, 1.0]
                            pl_limits = [0.0, 0.35]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Granitoid":
                            qz_limits = [0.6, 0.9]
                            kfs_limits = [0.0, 0.4]
                            pl_limits = [0.0, 0.4]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Quarzolite":
                            qz_limits = [0.9, 1.0]
                            kfs_limits = [0.0, 0.1]
                            pl_limits = [0.0, 0.1]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Norite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.1]
                            pl_limits = [0.72, 1.0]
                            bt_limits = [0.0, 0.05]
                        #
                        ## LOWER STRECKEISEN DIAGRAM
                        #
                        elif rock == "Foid-bearing Syenite":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.58, 1.0]
                            pl_limits = [0.0, 0.42]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid-bearing Monzonite":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.32, 0.65]
                            pl_limits = [0.32, 0.65]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid-bearing Monzodiorite":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.58, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid-bearing Monzogabbro":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.58, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid Monzosyenite":
                            upper_streckeisen = False
                            nph_limits = [0.1, 0.6]
                            kfs_limits = [0.2, 0.9]
                            pl_limits = [0.0, 0.45]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid Monzodiorite":
                            upper_streckeisen = False
                            nph_limits = [0.1, 0.6]
                            kfs_limits = [0.0, 0.45]
                            pl_limits = [0.2, 0.9]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid Monzogabbro":
                            upper_streckeisen = False
                            nph_limits = [0.1, 0.6]
                            kfs_limits = [0.0, 0.45]
                            pl_limits = [0.2, 0.9]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foidolite":
                            upper_streckeisen = False
                            nph_limits = [0.6, 1.0]
                            kfs_limits = [0.0, 0.4]
                            pl_limits = [0.0, 0.4]
                            bt_limits = [0.0, 0.05]

                        if upper_streckeisen == True:
                            phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                            phi_kfs = round(rd.uniform(kfs_limits[0], (1.0 - phi_qz)), 4)
                            phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                            phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                            phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                            #
                            if np.isclose(phi_total, 1.0000) == True:
                                if qz_limits[0] <= phi_qz <= qz_limits[1] \
                                        and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                        and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                        and bt_limits[0] <= phi_bt <= bt_limits[1]:
                                    condition_2 = True
                            #
                        elif upper_streckeisen == False:
                            phi_nph = round(rd.uniform(nph_limits[0], nph_limits[1]), 4)
                            phi_kfs = round(rd.uniform(kfs_limits[0], (1.0 - phi_nph)), 4)
                            phi_pl = round(rd.uniform(0.0, (1.0 - phi_nph - phi_kfs)), 4)
                            phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                            phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                            #
                            if np.isclose(phi_total, 1.0000) == True:
                                if nph_limits[0] <= phi_nph <= nph_limits[1] \
                                        and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                        and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                        and bt_limits[0] <= phi_bt <= bt_limits[1]:
                                    condition_2 = True
                        #
                    if upper_streckeisen == True:
                        phi_minerals["Qz"] = phi_qz
                    elif upper_streckeisen == False:
                        phi_minerals["Nph"] = phi_nph
                    #
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if self.porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                else:
                    var_porosity = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
                #
                rho = round((1 - var_porosity)*rho_s + var_porosity*data_fluid[2]/1000, 3)
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
                                    n_digits = 6
                                else:
                                    n_digits = 6
                                #
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 6)
                        #
                        w_elements[element] = round(w_elements[element], 6)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    #
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    #
                    condition = True
            #
            gamma_ray = 0.0
            photoelectricity = 0.0
            #
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
            #
            ## Bulk Density, Porosity, Seismic Velocities
            rho_solid = round(rho_s, 3)
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                rho_limits=[2500, 3000], vP_limits=[6000, 7500], vS_limits=[3500, 4000], delta=0.05,
                porosity=porosity)
            #
            ## Elastic Parameters
            bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                rho=rho, vP=vP, vS=vS)
            #
            ## Gamma Ray
            gamma_ray = round(gamma_ray, 3)
            ## Photoelectricity
            photoelectricity = round(photoelectricity, 3)
            # Composition data
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            amounts = []
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_elements = list(w_elements.keys())
            list_oxides = OxideCompounds(var_list_elements=list_elements).find_oxides()
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)

            if list_oxides[0] not in results_container["compounds"]:
                for oxide in list_oxides:
                    results_container["compounds"][oxide] = []

            for key, value in composition_oxides.items():
                results_container["compounds"][key].append(value)

            results_container["mineralogy"] = dict(sorted(
                results_container["mineralogy"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
            results_container["chemistry"] = dict(sorted(
                results_container["chemistry"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
            results_container["compounds"] = dict(sorted(
                results_container["compounds"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))

            # Results
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_s)
            results_container["rho"].append(rho)
            results_container["vP"].append(vP)
            results_container["vS"].append(vS)
            results_container["vP/vS"].append(vPvS)
            results_container["K"].append(bulk_modulus)
            results_container["G"].append(shear_modulus)
            results_container["E"].append(youngs_modulus)
            results_container["nu"].append(poisson_ratio)
            results_container["GR"].append(gamma_ray)
            results_container["PE"].append(photoelectricity)
            #
            n += 1
        #
        return results_container

    def create_felsic(self, w_Na=None, w_Mg=None, w_K=None, w_Ca=None, w_Fe = None, amounts=None):
        #
        self.w_Na = w_Na
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # Mineralogy + Fluids
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        amphibole = minerals.inosilicates.amphibole_na("")
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, amphibole]
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Na == None and self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                x = rd.uniform(0.0, 1.0)
                phi_qz = round(abs(-24.95*x**3 + 50.32*x**2 - 21.44*x + 96 - (574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82))/100, 4)
                phi_kfs = round(abs(574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82 - (45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4))/100, 4)
                phi_pl = round(abs(45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4 - (-23.27*x**3 + 54.32*x**2 + 0.37*x))/100, 4)
                phi_bt = round(abs(-23.27*x**3 + 54.32*x**2 + 0.37*x - (1.9*x**3 + 17.13*x**2 + 2.24*x))/100, 4)
                phi_ms = round(abs(100 - (-24.95*x**3 + 50.32*x**2 - 21.44*x + 96))/100, 4)
                phi_amph = round(abs(1.9*x**3 + 17.13*x**2 + 2.24*x)/100, 4)
                m_total = phi_qz*quartz[2] + phi_kfs*alkalifeldspar[2] + phi_pl*plagioclase[2] + phi_bt*biotite[2] + phi_ms*muscovite[2] + phi_amph*amphibole[2]
                w_qz = round(phi_qz*quartz[2]/m_total, 4)
                w_kfs = round(phi_kfs*alkalifeldspar[2]/m_total, 4)
                w_pl = round(phi_pl*plagioclase[2]/m_total, 4)
                w_bt = round(phi_bt*biotite[2]/m_total, 4)
                w_ms = round(phi_ms*muscovite[2]/m_total, 4)
                w_amph = round(phi_amph*amphibole[2]/m_total, 4)
            elif self.w_Na != None:
                w_kfs = round(abs(574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82 - (45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4))/100, 4)
                w_pl = round(abs(45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4 - (-23.27*x**3 + 54.32*x**2 + 0.37*x))/100, 4)
                w_amph = round((self.w_Na - w_kfs*alkalifeldspar[6][1] - w_pl*plagioclase[6][1])/(amphibole[6][2]), 4)
                w_qz = round(abs(-24.95*x**3 + 50.32*x**2 - 21.44*x + 96 - (574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82))/100, 4)
                w_bt = round(abs(-23.27*x**3 + 54.32*x**2 + 0.37*x - (1.9*x**3 + 17.13*x**2 + 2.24*x))/100, 4)
                w_ms = round(1-w_qz-w_kfs-w_pl-w_bt-w_amph, 4)
            elif self.w_Mg != None:
                w_bt = round((self.w_Mg)/(biotite[6][3]), 4)
                w_qz = round(abs(-24.95*x**3 + 50.32*x**2 - 21.44*x + 96 - (574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82))/100, 4)
                w_kfs = round(abs(574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82 - (45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4))/100, 4)
                w_pl = round(abs(45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4 - (-23.27*x**3 + 54.32*x**2 + 0.37*x))/100, 4)
                w_ms = round(abs(100 - (-24.95*x**3 + 50.32*x**2 - 21.44*x + 96))/100, 4)
                w_amph = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms, 4)
            elif self.w_K != None:
                w_bt = round(abs(-23.27*x**3 + 54.32*x**2 + 0.37*x - (1.9*x**3 + 17.13*x**2 + 2.24*x))/100, 4)
                w_ms = round(abs(100 - (-24.95*x**3 + 50.32*x**2 - 21.44*x + 96))/100, 4)
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_qz = round(abs(-24.95*x**3 + 50.32*x**2 - 21.44*x + 96 - (574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82))/100, 4)
                w_pl = round(abs(45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4 - (-23.27*x**3 + 54.32*x**2 + 0.37*x))/100, 4)
                w_amph = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms, 4)
            elif self.w_Ca != None:
                w_pl = round((self.w_Ca)/(plagioclase[6][4]), 4)
                w_qz = round(abs(-24.95*x**3 + 50.32*x**2 - 21.44*x + 96 - (574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82))/100, 4)
                w_kfs = round(abs(574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82 - (45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4))/100, 4)
                w_bt = round(abs(-23.27*x**3 + 54.32*x**2 + 0.37*x - (1.9*x**3 + 17.13*x**2 + 2.24*x))/100, 4)
                w_ms = round(abs(100 - (-24.95*x**3 + 50.32*x**2 - 21.44*x + 96))/100, 4)
                w_amph = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms, 4)
            elif self.w_Fe != None:
                w_bt = round(abs(-23.27*x**3 + 54.32*x**2 + 0.37*x - (1.9*x**3 + 17.13*x**2 + 2.24*x))/100, 4)
                w_amph = round((self.w_Fe - w_bt*biotite[6][7])/(amphibole[6][4]), 4)
                w_qz = round(abs(-24.95*x**3 + 50.32*x**2 - 21.44*x + 96 - (574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82))/100, 4)
                w_kfs = round(abs(574.08*x**5 - 1466.01*x**4 + 1220.19*x**3 - 285.57*x**2 - 33.17*x + 82 - (45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4))/100, 4)
                w_pl = round(abs(45.95*x**3 - 0.2*x**2 + 24.45*x + 13.4 - (-23.27*x**3 + 54.32*x**2 + 0.37*x))/100, 4)
                w_ms = round(1-w_qz-w_kfs-w_pl-w_bt-w_amph, 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.0025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.0025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.0025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.0025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.00025)), 4)
                w_amph = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_amph >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_amph, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_amph*amphibole[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1] + w_amph*amphibole[6][2], 4)
            w_Mg = round(w_bt*biotite[6][3], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_amph*amphibole[6][3], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_amph*amphibole[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Amph"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_amph]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                if 0.20 <= phi_V[0] <= 0.6 and 0.15 <= phi_V[1] <= 0.8 and 0.0 <= phi_V[2] <= 0.52:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_amph*amphibole[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_amph]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = 0.85*np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = 0.85*np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if self.porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 4000:
                phi = rd.uniform(0.0, 0.025)
        else:
            phi = self.porosity
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_amph*amphibole[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_amph*amphibole[5][1]
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_amph*amphibole[3][3]
        #
        if self.dict_output == False:
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("water")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
            #
            return data
        else:
            results = {}
            results["rock"] = "Felsic Rock"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Amph"]
            data.append(composition)
            results["chemistry"] = {}
            results["mineralogy"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = concentrations[index]
            for index, mineral in enumerate(mineral_list, start=0):
                results["mineralogy"][mineral] = amounts[index]
            #
            results["phi"] = round(phi, 4)
            results["fluid"] = self.fluid
            #
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            return results
    #
    def create_intermediate(self, amounts=None):
        #
        self.amounts = amounts
        #
        results = {}
        results["rock"] = "Intermediate Rock"
        #
        # Mineralogy + Fluids
        data_kfs = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_pl = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_bt = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        data_amph = Inosilicates(impurity="pure", data_type=True).create_calcium_amphibole()
        data_px = Inosilicates(impurity="pure", data_type=True).create_calium_pyroxene()
        #
        minerals_list = [self.data_qz, data_kfs, data_pl, data_bt, data_amph, data_px]
        water = fluids.Water.water("")
        #
        elements_list = []
        for mineral in minerals_list:
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements_list:
                    elements_list.append(element)
        elements_list.sort()
        #
        data = []
        #
        cond = False
        while cond == False:
            if self.amounts == None:
                x = rd.uniform(0.0, 1.0)
                phi_qz = round(abs(100 - (17.97*x**3 - 32.67*x**2 + 26.71*x + 91 - 3*x))/100, 4)
                phi_kfs = round(abs(17.97*x**3 - 32.67*x**2 + 26.71*x + 91 - 3*x - (121.55*x**5 - 358.74*x**4 + 402.63*x**3 - 224.97*x**2 + 77.15*x + 81.5))/100, 4)
                phi_pl = round(abs(121.55*x**5 - 358.74*x**4 + 402.63*x**3 - 224.97*x**2 + 77.15*x + 81.5 - (9.48*x**3 - 14.12*x**2 + 30.33*x + 31))/100, 4)
                phi_bt = round(abs(9.48*x**3 - 14.12*x**2 + 30.33*x + 31 - (-274.33*x**5 + 771.98*x**4 - 755.27*x**3 + 275.7*x**2 + 17.43*x + 20.9))/100, 4)
                phi_amph = round(abs(-274.33*x**5 + 771.98*x**4 - 755.27*x**3 + 275.7*x**2 + 17.43*x + 20.9 - (-2.07*x**3 + 30.77*x**2 + 0.5*x))/100, 4)
                phi_pyx = round(abs(-2.07*x**3 + 30.77*x**2 + 0.5*x)/100, 4)
                m_total = phi_qz*self.data_qz["rho"] + phi_kfs*data_kfs["rho"] + phi_pl*data_pl["rho"] + phi_bt*data_bt["rho"] + phi_amph*data_amph["rho"] + phi_pyx*data_px["rho"]
                w_qz = round(phi_qz*self.data_qz["rho"]/m_total, 4)
                w_kfs = round(phi_kfs*data_kfs["rho"]/m_total, 4)
                w_pl = round(phi_pl*data_pl["rho"]/m_total, 4)
                w_bt = round(phi_bt*data_bt["rho"]/m_total, 4)
                w_amph = round(phi_amph*data_amph["rho"]/m_total, 4)
                w_pyx = round(phi_pyx*data_px["rho"]/m_total, 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.00025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.0025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.0025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.0025)), 4)
                w_amph = round(abs(np.random.normal(self.amounts[4], 0.0025)), 4)
                w_pyx = round(1-w_qz-w_kfs-w_pl-w_bt-w_amph, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_amph >= 0.0 and w_pyx >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_amph + w_pyx, 4)
            else:
                sumMin = 0
            #
            mineral_amounts = {}
            mineral_amounts["Qz"] = w_qz
            mineral_amounts["Kfs"] = w_kfs
            mineral_amounts["Pl"] = w_pl
            mineral_amounts["Bt"] = w_bt
            mineral_amounts["Amph"] = w_amph
            mineral_amounts["Px"] = w_pyx
            #
            element_amounts = {}
            w_O = 1
            w_sum = 0
            for element in elements_list:
                if element != "O":
                    element_amounts[element] = 0
                    for mineral in minerals_list:
                        if element in mineral["chemistry"]:
                            element_amounts[element] += round(mineral["chemistry"][element]*mineral_amounts[mineral["mineral"]], 4)
                    element_amounts[element] = round(element_amounts[element], 4)
                    w_sum += element_amounts[element]
                    w_O -= element_amounts[element]
            element_amounts["O"] = round(w_O, 4)
            w_sum += element_amounts["O"]
            if sumMin == 1 and w_sum == 1:
                cond = True
            else:
                cond = False
        #
        results["mineralogy"] = mineral_amounts
        results["chemistry"] = element_amounts
        #
        rhoSolid = 0
        K_list = []
        G_list = []
        for mineral in minerals_list:
            rhoSolid += mineral["rho"]*mineral_amounts[mineral["mineral"]]/1000
            K_list.append(round(mineral["K"]*mineral_amounts[mineral["mineral"]], 3))
            G_list.append(round(mineral["G"]*mineral_amounts[mineral["mineral"]], 3))
        X = [w_qz, w_kfs, w_pl, w_bt, w_amph, w_pyx]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        #
        if self.porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 4000:
                phi = rd.uniform(0.0, 0.025)
        else:
            phi = self.porosity
        phi = round(phi, 4)
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        rho = (1 - phi)*rhoSolid + phi*water[2]/1000
        vP = (1 - phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi)*vS_solid
        #
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        #
        GR = 0
        PE = 0
        poisson_mineralogical = 0
        for mineral in minerals_list:
            GR += mineral["GR"]*mineral_amounts[mineral["mineral"]]
            PE += mineral["PE"]*mineral_amounts[mineral["mineral"]]
            poisson_mineralogical += mineral["nu"]*mineral_amounts[mineral["mineral"]]
        #
        if self.dict_output == False:
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(self.fluid)
            data.append([round(GR, 3), round(PE, 3)])
            data.append(amounts)
            #
            return data
        else:
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            return results

class Volcanic:
    #
    def __init__(self, fluid, actualThickness, dict_output=False, porosity=None):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.dict_output = dict_output
        self.porosity = porosity
        #
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        #
        self.data_water = Water.water("")
        self.data_oil = Hydrocarbons.oil("")
        self.data_gas = Hydrocarbons.natural_gas("")
    #
    def create_volcanic_rock_streckeisen(self, rock="Rhyolite", number=1, composition=None, enrichment_kfs=None,
                                         enrichment_pl=None, upper_streckeisen=True):
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
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
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar(
                enrichment=enrichment_kfs)
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(
                enrichment=enrichment_pl)
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            if upper_streckeisen == True:
                mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                              "Bt": data_biotite}
            elif upper_streckeisen == False:
                data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
                mineralogy = {"Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                              "Bt": data_biotite}
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
                    if upper_streckeisen == True:
                        phi_qz = composition["Qz"]
                        phi_minerals["Qz"] = phi_qz
                    elif upper_streckeisen == False:
                        phi_nph = composition["Nph"]
                        phi_minerals["Nph"] = phi_nph
                    #
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_bt = composition["Bt"]
                    #
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        #
                        ## UPPER STRECKEISEN DIAGRAM
                        #
                        if rock == "Rhyolite":
                            qz_limits = [0.2, 0.6]
                            kfs_limits = [0.15, 0.8]
                            pl_limits = [0.0, 0.52]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Dacite":
                            qz_limits = [0.2, 0.6]
                            kfs_limits = [0.0, 0.28]
                            pl_limits = [0.25, 0.8]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Trachyte":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.52, 1.0]
                            pl_limits = [0.0, 0.35]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Latite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.28, 0.65]
                            pl_limits = [0.28, 0.65]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Andesite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.52, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Basalt":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.52, 1.0]
                            bt_limits = [0.0, 0.05]
                        #
                        ## LOWER STRECKEISEN DIAGRAM
                        #
                        elif rock == "Foid-bearing Trachyte":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.58, 1.0]
                            pl_limits = [0.0, 0.35]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid-bearing Latite":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.32, 0.65]
                            pl_limits = [0.35, 0.68]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid-bearing Andesite":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.58, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foid-bearing Basalt":
                            upper_streckeisen = False
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.58, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Phonolite":
                            upper_streckeisen = False
                            nph_limits = [0.1, 0.6]
                            kfs_limits = [0.2, 0.9]
                            pl_limits = [0.0, 0.45]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Tephrite":
                            upper_streckeisen = False
                            nph_limits = [0.1, 0.6]
                            kfs_limits = [0.0, 0.45]
                            pl_limits = [0.2, 0.9]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Foidite":
                            upper_streckeisen = False
                            nph_limits = [0.6, 1.0]
                            kfs_limits = [0.0, 0.4]
                            pl_limits = [0.0, 0.4]
                            bt_limits = [0.0, 0.05]

                        if upper_streckeisen == True:
                            phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                            phi_kfs = round(rd.uniform(kfs_limits[0], (1.0 - phi_qz)), 4)
                            phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                            phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                            phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                            #
                            if np.isclose(phi_total, 1.0000) == True:
                                if qz_limits[0] <= phi_qz <= qz_limits[1] \
                                        and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                        and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                        and bt_limits[0] <= phi_bt <= bt_limits[1]:
                                    condition_2 = True
                            #
                        elif upper_streckeisen == False:
                            phi_nph = round(rd.uniform(nph_limits[0], nph_limits[1]), 4)
                            phi_kfs = round(rd.uniform(kfs_limits[0], (1.0 - phi_nph)), 4)
                            phi_pl = round(rd.uniform(0.0, (1.0 - phi_nph - phi_kfs)), 4)
                            phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                            phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                            #
                            if np.isclose(phi_total, 1.0000) == True:
                                if nph_limits[0] <= phi_nph <= nph_limits[1] \
                                        and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                        and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                        and bt_limits[0] <= phi_bt <= bt_limits[1]:
                                    condition_2 = True
                        #
                    if upper_streckeisen == True:
                        phi_minerals["Qz"] = phi_qz
                    elif upper_streckeisen == False:
                        phi_minerals["Nph"] = phi_nph
                    #
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if self.porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                else:
                    var_porosity = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
                #
                rho = round((1 - var_porosity)*rho_s + var_porosity*data_fluid[2]/1000, 3)
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
                                    n_digits = 6
                                else:
                                    n_digits = 6
                                #
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 6)
                        #
                        w_elements[element] = round(w_elements[element], 6)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    #
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    #
                    condition = True
            #
            bulk_mod = 0.0
            shear_mod = 0.0
            gamma_ray = 0.0
            photoelectricity = 0.0
            #
            K_list = []
            G_list = []
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            K_geo = elast.calc_geometric_mean(self, phi_list, K_list)
            G_geo = elast.calc_geometric_mean(self, phi_list, G_list)
            #
            # anisotropic_factor = round(rd.uniform(1.0, 2.0), 4)
            anisotropic_factor = 1.0
            #
            bulk_mod = K_geo/anisotropic_factor
            shear_mod = G_geo/anisotropic_factor
            #
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 6)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 6)
            # Composition data
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            amounts = []
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_elements = list(w_elements.keys())
            list_oxides = OxideCompounds(var_list_elements=list_elements).find_oxides()
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)

            if list_oxides[0] not in results_container["compounds"]:
                for oxide in list_oxides:
                    results_container["compounds"][oxide] = []

            for key, value in composition_oxides.items():
                results_container["compounds"][key].append(value)

            results_container["mineralogy"] = dict(sorted(
                results_container["mineralogy"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
            results_container["chemistry"] = dict(sorted(
                results_container["chemistry"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
            results_container["compounds"] = dict(sorted(
                results_container["compounds"].items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))

            # Results
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_s)
            results_container["rho"].append(rho)
            results_container["vP"].append(vP)
            results_container["vS"].append(vS)
            results_container["vP/vS"].append(vPvS)
            results_container["K"].append(bulk_mod)
            results_container["G"].append(shear_mod)
            results_container["E"].append(youngs_mod)
            results_container["nu"].append(poisson_rat)
            results_container["GR"].append(gamma_ray)
            results_container["PE"].append(photoelectricity)
            #
            n += 1
        #
        return results_container
    #
    def create_basaltic_rock_yoder_tilley(self, rock="Basalt", number=1, composition=None, enrichment_pl=None):
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
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
            data_olivine = Nesosilicates(impurity="pure", data_type=True).create_olivine()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(
                enrichment=enrichment_pl)
            data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
            data_orthopyroxene = Inosilicates(mineral="Orthopyroxene", data_type=True,
                                              traces_list=[]).generate_dataset(number=1)
            data_augite = Inosilicates(impurity="pure", data_type=True).create_augite()
            #
            mineralogy = {"Qz": self.data_quartz, "Pl": data_plagioclase, "Nph": data_nepheline,
                          "Opx": data_orthopyroxene, "Aug": data_augite, "Ol": data_olivine}
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
                    phi_qz = composition["Qz"]
                    phi_nph = composition["Nph"]
                    #
                    phi_opx = composition["Opx"]
                    phi_pl = composition["Pl"]
                    phi_ol = composition["Ol"]
                    phi_aug = composition["Aug"]
                    #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Nph"] = phi_nph
                    phi_minerals["Opx"] = phi_opx
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Ol"] = phi_ol
                    phi_minerals["Aug"] = phi_aug
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        magicnumber = rd.randint(0, 2)
                        if magicnumber == 0:        # Quartz-Tholeiite
                            qz_limits = [0.0, 1.0]
                            nph_limits = [0.0, 0.0]
                            opx_limits = [0.0, 1.0]
                            pl_limits = [0.0, 1.0]
                            ol_limits = [0.0, 0.0]
                            aug_limits = [0.0, 1.0]
                        elif magicnumber == 1:      # Olivine-Tholeiite
                            qz_limits = [0.0, 0.0]
                            nph_limits = [0.0, 0.0]
                            opx_limits = [0.0, 1.0]
                            pl_limits = [0.0, 1.0]
                            ol_limits = [0.0, 1.0]
                            aug_limits = [0.0, 1.0]
                        elif magicnumber == 2:      # Alkaline Basalt
                            qz_limits = [0.0, 0.0]
                            nph_limits = [0.0, 1.0]
                            opx_limits = [0.0, 0.0]
                            pl_limits = [0.0, 1.0]
                            ol_limits = [0.0, 1.0]
                            aug_limits = [0.0, 1.0]
                        #
                        phi_aug = round(rd.uniform(aug_limits[0], aug_limits[1]), 4)
                        phi_pl = round(rd.uniform(pl_limits[0], (1.0 - phi_aug)), 4)
                        phi_qz = round(rd.uniform(qz_limits[0], (1.0 - phi_aug - phi_pl)), 4)
                        phi_nph = round(rd.uniform(nph_limits[0], (1.0 - phi_aug - phi_pl - phi_qz)), 4)
                        phi_opx = round(rd.uniform(opx_limits[0], (1.0 - phi_aug - phi_pl - phi_qz - phi_nph)), 4)
                        phi_ol = round(1.0 - phi_aug - phi_pl - phi_qz - phi_nph - phi_opx, 4)
                        #
                        phi_total = phi_aug + phi_pl + phi_qz + phi_nph + phi_opx + phi_ol
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if aug_limits[0] <= phi_aug <= aug_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                    and qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and nph_limits[0] <= phi_nph <= nph_limits[1] \
                                    and opx_limits[0] <= phi_opx <= opx_limits[1] \
                                    and ol_limits[0] <= phi_ol <= ol_limits[1]:
                                condition_2 = True
                        #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Nph"] = phi_nph
                    #
                    phi_minerals["Aug"] = phi_aug
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Ol"] = phi_ol
                    phi_minerals["Opx"] = phi_opx
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    try:
                        rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                    except:
                        rho_s += phi_minerals[key]*mineralogy[key]["rho"][0]
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    try:
                        w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                    except:
                        w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"][0])/rho_s, n_digits)
                #
                if self.porosity == None:
                    var_porosity = round(rd.uniform(0.0, 0.1), 4)
                else:
                    var_porosity = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
                #
                rho = round((1 - var_porosity)*rho_s + var_porosity*data_fluid[2]/1000, 3)
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
                                    n_digits = 6
                                else:
                                    n_digits = 6
                                #
                                try:
                                    value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                except:
                                    value = round(w_mineral*mineralogy[mineral]["chemistry"][element][0], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 6)
                        #
                        w_elements[element] = round(w_elements[element], 6)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    #
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    #
                    condition = True
            #
            bulk_mod = 0.0
            shear_mod = 0.0
            gamma_ray = 0.0
            photoelectricity = 0.0
            #
            K_list = []
            G_list = []
            phi_list = []
            for key, value in phi_minerals.items():
                try:
                    gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                    photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                    #
                    gamma_ray = round(gamma_ray, 3)
                    photoelectricity = round(photoelectricity, 3)
                    #
                    K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                    G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                    phi_list.append(phi_minerals[key])
                except:
                    gamma_ray += phi_minerals[key]*mineralogy[key]["GR"][0]
                    photoelectricity += phi_minerals[key]*mineralogy[key]["PE"][0]
                    #
                    gamma_ray = round(gamma_ray, 3)
                    photoelectricity = round(photoelectricity, 3)
                    #
                    K_list.append(round(phi_minerals[key]*mineralogy[key]["K"][0], 3))
                    G_list.append(round(phi_minerals[key]*mineralogy[key]["G"][0], 3))
                    phi_list.append(phi_minerals[key])
            #
            K_geo = elast.calc_geometric_mean(self, phi_list, K_list)
            G_geo = elast.calc_geometric_mean(self, phi_list, G_list)
            #
            # anisotropic_factor = round(rd.uniform(1.0, 2.0), 4)
            anisotropic_factor = 1.0
            #
            bulk_mod = K_geo/anisotropic_factor
            shear_mod = G_geo/anisotropic_factor
            #
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 6)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 6)
            # Composition data
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            amounts = []
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_elements = list(w_elements.keys())
            list_oxides = OxideCompounds(var_list_elements=list_elements).find_oxides()
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)

            if list_oxides[0] not in results_container["compounds"]:
                for oxide in list_oxides:
                    results_container["compounds"][oxide] = []

            for key, value in composition_oxides.items():
                results_container["compounds"][key].append(value)

            # Results
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_s)
            results_container["rho"].append(rho)
            results_container["vP"].append(vP)
            results_container["vS"].append(vS)
            results_container["vP/vS"].append(vPvS)
            results_container["K"].append(bulk_mod)
            results_container["G"].append(shear_mod)
            results_container["E"].append(youngs_mod)
            results_container["nu"].append(poisson_rat)
            results_container["GR"].append(gamma_ray)
            results_container["PE"].append(photoelectricity)
            #
            n += 1
        #
        return results_container

class Pyroclastic:
    #
    def __init__(self, fluid, actualThickness):
        self.fluid = fluid
        self.actualThickness = actualThickness
        #
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_water = Water.water("")
    #
    def create_pyroclastic_rock(self, number=1, composition=None, porosity=None):
        data_olivine = Nesosilicates(impurity="pure", data_type=True).create_olivine()
        data_orthopyroxene = Inosilicates(impurity="pure", data_type=True).create_orthopyroxene()
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        #
        mineralogy = {"Qz": self.data_quartz, "Ol": data_olivine, "Opx": data_orthopyroxene, "Nph": data_nepheline}
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
                phi_qz = composition["Qz"]
                phi_Ol = composition["Ol"]
                phi_Opx = composition["Opx"]
                phi_Nph = composition["Nph"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Ol"] = phi_Ol
                phi_minerals["Opx"] = phi_Opx
                phi_minerals["Nph"] = phi_Nph
            else:
                condition_2 = False
                while condition_2 == False:
                    x = rd.randint(0, 1)
                    phi_qz = round(x*rd.uniform(0.2, 0.6), 4)
                    phi_nph = round((1 - x)*rd.uniform(0.2, 0.6), 4)
                    phi_ol = round(rd.uniform(0.25, (1.0 - phi_qz)), 4)
                    phi_opx = round(1 - phi_qz - phi_ol - phi_nph, 4)
                    phi_total = phi_qz + phi_ol + phi_opx + phi_nph
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.6 and 0.25 <= phi_ol <= 0.75 and 0.05 <= phi_opx <= 0.15 and 0.0 <= phi_nph <= 0.6:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Ol"] = phi_ol
                phi_minerals["Opx"] = phi_opx
                phi_minerals["Nph"] = phi_nph
            #
            rho_s = 0
            for key, value in phi_minerals.items():
                rho_s += phi_minerals[key] * mineralogy[key]["rho"]
                for element, value in mineralogy[key]["chemistry"].items():
                    if element not in elements_list:
                        elements_list.append(element)
                        w_elements[element] = 0.0
            rho_s = round(rho_s, 3)
            for key, value in phi_minerals.items():
                w_minerals[key] = round((phi_minerals[key] * mineralogy[key]["rho"]) / rho_s, 4)
            #
            if porosity == None:
                porosity = round(rd.uniform(0.1, 0.5), 4)
            rho = round((1 - porosity) * rho_s + porosity * self.data_water[2] / 1000, 3)
            #
            old_index = elements_list.index("O")
            elements_list += [elements_list.pop(old_index)]
            #
            w_elements_total = 0.0
            for element in elements_list:
                if element != "O":
                    for mineral, w_mineral in w_minerals.items():
                        if element in mineralogy[mineral]["chemistry"]:
                            value = round(w_mineral * mineralogy[mineral]["chemistry"][element], 4)
                            w_elements[element] += value
                            w_elements_total += value
                            #
                            w_elements[element] = round(w_elements[element], 4)
                elif element == "O":
                    w_elements[element] += round(1 - w_elements_total, 4)
                    #
                    w_elements[element] = round(w_elements[element], 4)
            #
            if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                condition = True
        #
        bulk_mod = 0.0
        shear_mod = 0.0
        gamma_ray = 0.0
        photoelectricity = 0.0
        for key, value in phi_minerals.items():
            bulk_mod += phi_minerals[key] * mineralogy[key]["K"]
            shear_mod += phi_minerals[key] * mineralogy[key]["G"]
            gamma_ray += phi_minerals[key] * mineralogy[key]["GR"]
            photoelectricity += phi_minerals[key] * mineralogy[key]["PE"]
            #
            bulk_mod = round(bulk_mod, 3)
            shear_mod = round(shear_mod, 3)
            gamma_ray = round(gamma_ray, 3)
            photoelectricity = round(photoelectricity, 3)
        #
        youngs_mod = round((9 * bulk_mod * shear_mod) / (3 * bulk_mod + shear_mod), 3)
        poisson_rat = round((3 * bulk_mod - 2 * shear_mod) / (6 * bulk_mod + 2 * shear_mod), 4)
        vP = round(((bulk_mod * 10 ** 9 + 4 / 3 * shear_mod * 10 ** 9) / (rho)) ** 0.5, 3)
        vS = round(((shear_mod * 10 ** 9) / (rho)) ** 0.5, 3)
        vPvS = round(vP / vS, 3)
        #
        results = {}
        results["rock"] = "Pyroclastic Rock"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = porosity
            results["fluid"] = "water"
            results["rho_s"] = rho_s
            results["rho"] = rho
            results["vP"] = vP
            results["vS"] = vS
            results["vP/vS"] = vPvS
            results["K"] = bulk_mod
            results["G"] = shear_mod
            results["E"] = youngs_mod
            results["nu"] = poisson_rat
            results["GR"] = gamma_ray
            results["PE"] = photoelectricity
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in w_minerals.items():
                single_amounts_mineralogy[mineral] = value
            for element, value in w_elements.items():
                single_amounts_chemistry[element] = value
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = porosity
            results["fluid"] = "water"
            results["rho_s"] = rho_s
            results["rho"] = rho
            results["vP"] = vP
            results["vS"] = vS
            results["vP/vS"] = vPvS
            results["K"] = bulk_mod
            results["G"] = shear_mod
            results["E"] = youngs_mod
            results["nu"] = poisson_rat
            results["GR"] = gamma_ray
            results["PE"] = photoelectricity
        #
        return results