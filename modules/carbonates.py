#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		carbonates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		10.12.2024

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals
from modules import fluids
from modules.geophysics import Elasticity as elast
from modules.chemistry import PeriodicSystem, OxideCompounds
from modules.minerals import CrystalPhysics
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry
from modules.minerals import Organics
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Phyllosilicates
from modules.sulfides import Sulfides
from modules.pyllosilicates import Pyllosilicates
from modules.sulfates import Sulfates
from modules.fluids import Water
from modules.petrophysics import SeismicVelocities

## CLASSES
class CarbonateRocks:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        #
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_aragonite = Carbonates(impurity="pure", data_type=True).create_aragonite()
        self.data_dolomite = Carbonates(impurity="pure", data_type=True).create_dolomite()
        self.data_siderite = Carbonates(impurity="pure", data_type=True).create_siderite()
        self.data_magnesite = Carbonates(impurity="pure", data_type=True).create_magnesite()
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        self.data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        self.data_pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()
        #
        self.data_water = Water.water("")

    def create_limestone(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
        #
        mineralogy = {"Cal": self.data_calcite, "Dol": self.data_dolomite, "Qz": self.data_quartz,
                      "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Kln": self.data_kaolinite,
                      "Ilt": data_illite, "Py": self.data_pyrite}
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        amounts_compounds = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral, dataset in mineralogy.items():
            amounts_mineralogy[dataset["mineral"]] = []
            mineral_list.append(dataset["mineral"])
            elements_mineral = list(dataset["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        while n < number:
            condition = False
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}
                #
                if composition != None:
                    phi_cal = composition["Cal"]
                    phi_dol = composition["Dol"]
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_kln = composition["Kln"]
                    phi_ilt = composition["Ilt"]
                    phi_py = composition["Py"]
                    #
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Dol"] = phi_dol
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Py"] = phi_py
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        magicnumber = rd.randint(0, 12)
                        if 0 <= magicnumber <= 8:   # Carbonate-dominated
                            w_carb = round(rd.uniform(0.95, 1.0), 4)
                            w_clast = round(rd.uniform(0.0, (1.0 - w_carb)), 4)
                            w_clay = round(rd.uniform(0.0, (1.0 - w_carb - w_clast)), 4)
                            w_sulf = round(1 - w_carb - w_clast - w_clay, 4)
                            #
                            phi_cal = round(w_carb*rd.uniform(0.975, 1.0), 4)
                            phi_dol = round(w_carb - phi_cal, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.0, 1.0), 4)
                            phi_pl = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(w_clast - phi_qz - phi_pl, 4)
                            #
                            magicnumber_2 = rd.randint(0, 1)
                            if magicnumber_2 == 0:
                                phi_kln = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_ilt = round(w_clay - phi_kln, 4)
                            else:
                                phi_ilt = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_kln = round(w_clay - phi_ilt, 4)
                            #
                            phi_py = round(w_sulf, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_qz + phi_kfs + phi_pl + phi_kln + phi_ilt + phi_py
                            #
                        elif magicnumber in [9, 10]:   # Clastic-dominated
                            w_clast = round(rd.uniform(0.05, 0.10), 4)
                            w_carb = round(rd.uniform(0.90, (1.0 - w_clast)), 4)
                            w_clay = round(rd.uniform(0.0, (1.0 - w_carb - w_clast)), 4)
                            w_sulf = round(1 - w_carb - w_clast - w_clay, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.0, 1.0), 4)
                            phi_pl = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(w_clast - phi_qz - phi_pl, 4)
                            #
                            phi_cal = round(w_carb*rd.uniform(0.975, 1.0), 4)
                            phi_dol = round(w_carb - phi_cal, 4)
                            #
                            magicnumber_2 = rd.randint(0, 1)
                            if magicnumber_2 == 0:
                                phi_kln = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_ilt = round(w_clay - phi_kln, 4)
                            else:
                                phi_ilt = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_kln = round(w_clay - phi_ilt, 4)
                            #
                            phi_py = round(w_sulf, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_qz + phi_kfs + phi_pl + phi_kln + phi_ilt + phi_py
                            #
                        elif magicnumber in [11, 12]:   # Clay-dominated
                            w_clay = round(rd.uniform(0.05, 0.10), 4)
                            w_carb = round(rd.uniform(0.90, (1.0 - w_clay)), 4)
                            w_clast = round(rd.uniform(0.0, (1.0 - w_carb - w_clay)), 4)
                            w_sulf = round(1 - w_carb - w_clast - w_clay, 4)
                            #
                            phi_cal = round(w_carb*rd.uniform(0.975, 1.0), 4)
                            phi_dol = round(w_carb - phi_cal, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.0, 1.0), 4)
                            phi_pl = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(w_clast - phi_qz - phi_pl, 4)
                            #
                            magicnumber_2 = rd.randint(0, 1)
                            if magicnumber_2 == 0:
                                phi_kln = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_ilt = round(w_clay - phi_kln, 4)
                            else:
                                phi_ilt = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_kln = round(w_clay - phi_ilt, 4)
                            #
                            phi_py = round(w_sulf, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_qz + phi_kfs + phi_pl + phi_kln + phi_ilt + phi_py
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.9 <= phi_cal <= 1.0 and 0.0 <= phi_dol <= 0.2 and 0.0 <= phi_qz <= 0.3 \
                                    and 0.0 <= phi_kfs <= 0.2 and 0.0 <= phi_pl <= 0.2 \
                                    and 0.0 <= phi_kln <= 0.3 and 0.0 <= phi_ilt <= 0.3 and 0.0 <= phi_py <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Cal"] = abs(phi_cal)
                    phi_minerals["Dol"] = abs(phi_dol)
                    phi_minerals["Qz"] = abs(phi_qz)
                    phi_minerals["Kfs"] = abs(phi_kfs)
                    phi_minerals["Pl"] = abs(phi_pl)
                    phi_minerals["Kln"] = abs(phi_kln)
                    phi_minerals["Ilt"] = abs(phi_ilt)
                    phi_minerals["Py"] = abs(phi_py)
                #
                rho_s = 0
                velocity_solid = {"vP": 0, "vS": 0}
                gamma_ray = 0.0
                photoelectricity = 0.0
                for key, value in phi_minerals.items():
                    rho_s += value*mineralogy[key]["rho"]
                    velocity_solid["vP"] += value*mineralogy[key]["vP"]
                    velocity_solid["vS"] += value*mineralogy[key]["vS"]
                    gamma_ray += value*mineralogy[key]["GR"]
                    photoelectricity += value*mineralogy[key]["PE"]
                    #
                    for element, value in mineralogy[key]["chemistry"].items():
                        if element not in elements_list:
                            elements_list.append(element)
                            w_elements[element] = 0.0
                #
                rho_s = round(rho_s, 3)
                #
                for key, value in phi_minerals.items():
                    w_result = round((phi_minerals[key] * mineralogy[key]["rho"]) / rho_s, 4)
                    w_minerals[key] = w_result
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
                    #
                    condition = True
                #
                ## Bulk Density, Porosity, Seismic Velocities
                rho_solid = round(rho_s, 3)
                vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                    rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                    rho_limits=[1800, 2800], vP_limits=[3500, 6500], vS_limits=[2000, 3500], delta=0.05,
                    porosity=porosity)
                ## Elastic Parameters
                bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                    rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                    rho=rho, vP=vP, vS=vS)
                ## Gamma Ray
                gamma_ray = round(gamma_ray, 3)
                ## Photoelectricity
                photoelectricity = round(photoelectricity, 3)

            for mineral, value in w_minerals.items():
                amounts_mineralogy[mineral].append(float(value))
            for element, value in w_elements.items():
                amounts_chemistry[element].append(float(value))

            bulk_properties["rho_s"].append(float(rho_solid))
            bulk_properties["rho"].append(float(rho))
            bulk_properties["phi"].append(float(var_porosity))
            bulk_properties["K"].append(float(bulk_modulus))
            bulk_properties["G"].append(float(shear_modulus))
            bulk_properties["E"].append(float(youngs_modulus))
            bulk_properties["nu"].append(float(poisson_ratio))
            bulk_properties["vP"].append(float(vP))
            bulk_properties["vS"].append(float(vS))
            bulk_properties["vPvS"].append(float(vPvS))
            bulk_properties["GR"].append(float(gamma_ray))
            bulk_properties["PE"].append(float(photoelectricity))

            amounts = []
            for key, value in w_elements.items():
                amounts_chemistry[key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_oxides = ["H2O", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "SO3", "K2O", "CaO", "Fe2O3"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)

            if list_oxides[0] not in amounts_compounds:
                for oxide in list_oxides:
                    amounts_compounds[oxide] = []

            for key, value in composition_oxides.items():
                amounts_compounds[key].append(value)

            n += 1
        amounts_mineralogy = dict(sorted(
            amounts_mineralogy.items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
        amounts_chemistry = dict(sorted(
            amounts_chemistry.items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
        amounts_compounds = dict(sorted(
            amounts_compounds.items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))

        ## EXPORT DATA
        #
        results = {}
        results["rock"] = "Limestone"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["compounds"] = amounts_compounds
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
            single_amounts_compounds = {}
            for mineral, value in w_minerals.items():
                single_amounts_mineralogy[mineral] = float(value)
            for element, value in w_elements.items():
                single_amounts_chemistry[element] = float(value)
            for compound, value in amounts_compounds.items():
                single_amounts_compounds[compound] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["compounds"] = single_amounts_compounds
            results["phi"] = float(var_porosity)
            results["fluid"] = "water"
            results["rho_s"] = float(rho_solid)
            results["rho"] = float(rho)
            results["vP"] = float(vP)
            results["vS"] = float(vS)
            results["vP/vS"] = float(vPvS)
            results["K"] = float(bulk_modulus)
            results["G"] = float(shear_modulus)
            results["E"] = float(youngs_modulus)
            results["nu"] = float(poisson_ratio)
            results["GR"] = float(gamma_ray)
            results["PE"] = float(photoelectricity)
        #
        return results
    #
    def create_dolostone(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
        #
        mineralogy = {"Cal": self.data_calcite, "Dol": self.data_dolomite, "Qz": self.data_quartz,
                      "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Kln": self.data_kaolinite,
                      "Ilt": data_illite, "Py": self.data_pyrite}
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        amounts_compounds = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral, dataset in mineralogy.items():
            amounts_mineralogy[dataset["mineral"]] = []
            mineral_list.append(dataset["mineral"])
            elements_mineral = list(dataset["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        while n < number:
            condition = False
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}
                #
                if composition != None:
                    phi_cal = composition["Cal"]
                    phi_dol = composition["Dol"]
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_kln = composition["Kln"]
                    phi_ilt = composition["Ilt"]
                    phi_py = composition["Py"]
                    #
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Dol"] = phi_dol
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Py"] = phi_py
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        magicnumber = rd.randint(0, 12)
                        if 0 <= magicnumber <= 8:  # Carbonate-dominated
                            w_carb = round(rd.uniform(0.8, 1.0), 4)
                            w_clast = round(rd.uniform(0.0, (1.0 - w_carb)), 4)
                            w_clay = round(rd.uniform(0.0, (1.0 - w_carb - w_clast)), 4)
                            w_sulf = round(1 - w_carb - w_clast - w_clay, 4)
                            #
                            phi_dol = round(w_carb*rd.uniform(0.9, 1.0), 4)
                            phi_cal = round(w_carb - phi_dol, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.8, 1.0), 4)
                            phi_pl = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(w_clast - phi_qz - phi_pl, 4)
                            #
                            magicnumber_2 = rd.randint(0, 1)
                            if magicnumber_2 == 0:
                                phi_kln = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_ilt = round(w_clay - phi_kln, 4)
                            else:
                                phi_ilt = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_kln = round(w_clay - phi_ilt, 4)
                            #
                            phi_py = round(w_sulf, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_qz + phi_kfs + phi_pl + phi_kln + phi_ilt + phi_py
                            #
                        elif magicnumber in [9, 10]:  # Clastic-dominated
                            w_clast = round(rd.uniform(0.15, 0.3), 4)
                            w_carb = round(rd.uniform(0.7, (1.0 - w_clast)), 4)
                            w_clay = round(rd.uniform(0.0, (1.0 - w_carb - w_clast)), 4)
                            w_sulf = round(1 - w_carb - w_clast - w_clay, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.8, 1.0), 4)
                            phi_pl = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(w_clast - phi_qz - phi_pl, 4)
                            #
                            phi_dol = round(w_carb*rd.uniform(0.9, 1.0), 4)
                            phi_cal = round(w_carb - phi_dol, 4)
                            #
                            magicnumber_2 = rd.randint(0, 1)
                            if magicnumber_2 == 0:
                                phi_kln = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_ilt = round(w_clay - phi_kln, 4)
                            else:
                                phi_ilt = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_kln = round(w_clay - phi_ilt, 4)
                            #
                            phi_py = round(w_sulf, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_qz + phi_kfs + phi_pl + phi_kln + phi_ilt + phi_py
                            #
                        elif magicnumber in [11, 12]:  # Clay-dominated
                            w_clay = round(rd.uniform(0.15, 0.3), 4)
                            w_carb = round(rd.uniform(0.7, (1.0 - w_clay)), 4)
                            w_clast = round(rd.uniform(0.0, (1.0 - w_carb - w_clay)), 4)
                            w_sulf = round(1 - w_carb - w_clast - w_clay, 4)
                            #
                            phi_dol = round(w_carb*rd.uniform(0.9, 1.0), 4)
                            phi_cal = round(w_carb - phi_dol, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.8, 1.0), 4)
                            phi_pl = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(w_clast - phi_qz - phi_pl, 4)
                            #
                            magicnumber_2 = rd.randint(0, 1)
                            if magicnumber_2 == 0:
                                phi_kln = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_ilt = round(w_clay - phi_kln, 4)
                            else:
                                phi_ilt = round(w_clay*rd.uniform(0.9, 1.0), 4)
                                phi_kln = round(w_clay - phi_ilt, 4)
                            #
                            phi_py = round(w_sulf, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_qz + phi_kfs + phi_pl + phi_kln + phi_ilt + phi_py
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.8 <= phi_dol <= 1.0 and 0.0 <= phi_cal <= 0.2 and 0.0 <= phi_qz <= 0.3 \
                                    and 0.0 <= phi_kfs <= 0.2 and 0.0 <= phi_pl <= 0.2 \
                                    and 0.0 <= phi_kln <= 0.3 and 0.0 <= phi_ilt <= 0.3 and 0.0 <= phi_py <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Cal"] = abs(phi_cal)
                    phi_minerals["Dol"] = abs(phi_dol)
                    phi_minerals["Qz"] = abs(phi_qz)
                    phi_minerals["Kfs"] = abs(phi_kfs)
                    phi_minerals["Pl"] = abs(phi_pl)
                    phi_minerals["Kln"] = abs(phi_kln)
                    phi_minerals["Ilt"] = abs(phi_ilt)
                    phi_minerals["Py"] = abs(phi_py)
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += value*mineralogy[key]["rho"]
                    for element, value in mineralogy[key]["chemistry"].items():
                        if element not in elements_list:
                            elements_list.append(element)
                            w_elements[element] = 0.0
                #
                rho_s = round(rho_s, 3)
                #
                for key, value in phi_minerals.items():
                    w_result = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, 4)
                    w_minerals[key] = w_result
                #
                if porosity == None:
                    phi_helper = round(rd.uniform(0.0, 0.4), 4)
                else:
                    phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
                #
                #rho = round((1 - phi_helper)*rho_s + phi_helper*self.data_water[2], 3)
                #
                old_index = elements_list.index("O")
                elements_list += [elements_list.pop(old_index)]
                #
                w_elements_total = 0.0
                for element in elements_list:
                    if element != "O":
                        for mineral, w_mineral in w_minerals.items():
                            if element in mineralogy[mineral]["chemistry"]:
                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], 4)
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
                    #
                    condition = True
                #
                bulk_mod = 0.0
                shear_mod = 0.0
                gamma_ray = 0.0
                photoelectricity = 0.0
                for key, value in phi_minerals.items():
                    bulk_mod += phi_minerals[key]*mineralogy[key]["K"]
                    shear_mod += phi_minerals[key]*mineralogy[key]["G"]
                    gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                    photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                    #
                    bulk_mod = round(bulk_mod, 3)
                    shear_mod = round(shear_mod, 3)
                    gamma_ray = round(gamma_ray, 3)
                    photoelectricity = round(photoelectricity, 3)
                #
                ## Bulk Density, Porosity, Seismic Velocities
                rho_solid = round(rho_s, 3)
                vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                    rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                    rho_limits=[2000, 2950], vP_limits=[3500, 7000], vS_limits=[2200, 3200], delta=0.05,
                    porosity=porosity)
                ## Elastic Parameters
                bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                    rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                    rho=rho, vP=vP, vS=vS)
            #
            for mineral, value in w_minerals.items():
                amounts_mineralogy[mineral].append(value)
            for element, value in w_elements.items():
                amounts_chemistry[element].append(value)
            #
            bulk_properties["rho_s"].append(rho_solid)
            bulk_properties["rho"].append(rho)
            bulk_properties["phi"].append(var_porosity)
            bulk_properties["K"].append(bulk_modulus)
            bulk_properties["G"].append(shear_modulus)
            bulk_properties["E"].append(youngs_modulus)
            bulk_properties["nu"].append(poisson_ratio)
            bulk_properties["vP"].append(vP)
            bulk_properties["vS"].append(vS)
            bulk_properties["vPvS"].append(vPvS)
            bulk_properties["GR"].append(gamma_ray)
            bulk_properties["PE"].append(photoelectricity)

            amounts = []
            for key, value in w_elements.items():
                amounts_chemistry[key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_oxides = ["H2O", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "SO3", "K2O", "CaO", "Fe2O3"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)

            if list_oxides[0] not in amounts_compounds:
                for oxide in list_oxides:
                    amounts_compounds[oxide] = []

            for key, value in composition_oxides.items():
                amounts_compounds[key].append(value)

            n += 1

        amounts_mineralogy = dict(sorted(
            amounts_mineralogy.items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
        amounts_chemistry = dict(sorted(
            amounts_chemistry.items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))
        amounts_compounds = dict(sorted(
            amounts_compounds.items(), key=lambda item: sum(item[1])/len(item[1]), reverse=True))

        ## EXPORT DATA
        results = {}
        results["rock"] = "Dolostone"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["compounds"] = amounts_compounds
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
            single_amounts_compounds = {}
            for mineral, value in w_minerals.items():
                single_amounts_mineralogy[mineral] = value
            for element, value in w_elements.items():
                single_amounts_chemistry[element] = value
            for compound, value in amounts_compounds.items():
                single_amounts_compounds[compound] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["compounds"] = single_amounts_compounds
            results["phi"] = var_porosity
            results["fluid"] = "water"
            results["rho_s"] = rho_solid
            results["rho"] = rho
            results["vP"] = vP
            results["vS"] = vS
            results["vP/vS"] = vPvS
            results["K"] = bulk_modulus
            results["G"] = shear_modulus
            results["E"] = youngs_modulus
            results["nu"] = poisson_ratio
            results["GR"] = gamma_ray
            results["PE"] = photoelectricity
        #
        return results

    def create_marl(self, rock="Marl", number=1, composition=None, classification="Marl", porosity=None):
        #
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
            ## Bulk Density, Porosity, Seismic Velocities
            rho_solid = round(rho_s, 3)
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                rho_limits=[2200, 2700], vP_limits=[1800, 3200], vS_limits=[1000, 2000], delta=0.05,
                porosity=porosity)
            ## Elastic Parameters
            bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                rho=rho, vP=vP, vS=vS)
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

            list_oxides = ["H2O", "CO2", "Na2O", "Al2O3", "SiO2", "K2O", "CaO"]
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

# Carbonates
class Carbonates:
    """ Class that generates geophysical and geochemical data of carbonate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1):
        if self.mineral in ["Cal", "Calcite"]:
            if number > 1:
                data = [self.create_calcite() for n in range(number)]
            else:
                data = self.create_calcite()
        elif self.mineral in ["Dol", "Dolomite"]:
            if number > 1:
                data = [self.create_dolomite() for n in range(number)]
            else:
                data = self.create_dolomite()
        elif self.mineral in ["Mgs", "Magnesite"]:
            if number > 1:
                data = [self.create_magnesite() for n in range(number)]
            else:
                data = self.create_magnesite()
        elif self.mineral in ["Sd", "Siderite"]:
            if number > 1:
                data = [self.create_siderite() for n in range(number)]
            else:
                data = self.create_siderite()
        elif self.mineral in ["Org", "Organic Matter"]:
            if number > 1:
                data = [self.create_organic_matter() for n in range(number)]
            else:
                data = self.create_organic_matter()
        elif self.mineral in ["Rdc", "Rhodochrosite"]:
            if number > 1:
                data = [self.create_rhodochrosite() for n in range(number)]
            else:
                data = self.create_rhodochrosite()
        elif self.mineral in ["Arg", "Aragonite"]:
            if number > 1:
                data = [self.create_aragonite() for n in range(number)]
            else:
                data = self.create_aragonite()
        elif self.mineral in ["Cer", "Cerussite"]:
            if number > 1:
                data = [self.create_cerussite() for n in range(number)]
            else:
                data = self.create_cerussite()
        elif self.mineral in ["Ank", "Ankerite"]:
            if number > 1:
                data = [self.create_ankerite() for n in range(number)]
            else:
                data = self.create_ankerite()
        elif self.mineral in ["Az", "Azurite"]:
            if number > 1:
                data = [self.create_azurite() for n in range(number)]
            else:
                data = self.create_azurite()
        elif self.mineral in ["Mal", "Malachite"]:
            if number > 1:
                data = [self.create_malachite() for n in range(number)]
            else:
                data = self.create_malachite()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Calcite":
                data_mineral = self.create_calcite()
            elif self.mineral == "Dolomite":
                data_mineral = self.create_dolomite()
            elif self.mineral == "Magnesite":
                data_mineral = self.create_magnesite()
            elif self.mineral == "Siderite":
                data_mineral = self.create_siderite()
            elif self.mineral == "Rhodochrosite":
                data_mineral = self.create_rhodochrosite()
            elif self.mineral == "Aragonite":
                data_mineral = self.create_aragonite()
            elif self.mineral == "Cerrusite":
                data_mineral = self.create_cerussite()
            elif self.mineral == "Ankerite":
                data_mineral = self.create_ankerite()
            elif self.mineral == "Azurite":
                data_mineral = self.create_azurite()
            elif self.mineral == "Malachite":
                data_mineral = self.create_malachite()
            elif self.mineral == "Ikaite":
                data_mineral = self.create_ikaite()
            elif self.mineral == "Smithsonite":
                data_mineral = self.create_smithsonite()
            #
            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U",
                           "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif key in ["mineral", "state", "trace elements"] and key not in dataset:
                    dataset[key] = value
                elif key in ["chemistry"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)
            #
        return dataset
    #
    def create_calcite(self):   # CaCO3
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Ca"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Fe", "Zn", "Co", "Ba", "Sr", "Pb", "Mg", "Cu", "Al", "Ni", "V", "Cr", "Mo"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        data = []
        mineral = "Cal"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + carbon[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.99, 17.06], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 76*10**9
        # Shear modulus
        G = 32*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 2*10**12
        #
        if self.data_type == False:
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_dolomite(self):   # CaMg(CO3)2
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Mg", "Ca"]
        majors_data = np.array([["C", carbon[1], 2, carbon[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Mg", magnesium[1], 1, magnesium[2]], ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Co", "Pb", "Zn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        data = []
        mineral = "Dol"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + magnesium[2] + 2*(carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.81, 16.01], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 3
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 89*10**9
        # Shear modulus
        G = 44*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 5*10**3
        #
        if self.data_type == False:
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_magnesite(self):   # Mg(CO3)
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        majors_name = ["C", "O", "Mg"]
        majors_data = np.array([["C", carbon[1], 2, carbon[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Mg", magnesium[1], 1, magnesium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Ca", "Co", "Ni"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        data = []
        mineral = "Mgs"
        #
        # Molar mass
        molar_mass_pure = magnesium[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.63, 15.03], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 105*10**9
        # Shear modulus
        G = 63*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_siderite(self):   # Fe(CO3)
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["C", "O", "Fe"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Fe", iron[1], 1, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Mg", "Ca", "Zn", "Co"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Sd"
        #
        # Molar mass
        molar_mass_pure = iron[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.69, 15.38], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 124*10**9
        # Shear modulus
        G = 51*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 70
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_organic_matter(self):
        # CHEMISTRY
        carbohydrates = Organics.carbohydrates("")
        lignin = Organics.lignin("")
        lipid = Organics.lipid("")
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        nitrogen = PeriodicSystem(name="N").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sulfur = PeriodicSystem(name="S").get_data()
        majors_name = ["H", "C", "N", "O", "S"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Org"
        #
        # Molar mass
        condition = False
        while condition == False:
            w_ch = round(rd.uniform(0.4, 0.6), 4)
            w_lg = round(rd.uniform(0.2, float(1-w_ch)), 4)
            w_lp = round(1 - w_ch - w_lg, 4)
            if w_ch+w_lg+w_lp == 1.0:
                condition = True
        molar_mass_pure = w_ch*(0.06*hydrogen[2] + 0.44*carbon[2] + 0.50*oxygen[2]) \
                          + w_lg*(0.06*hydrogen[2] + 0.63*carbon[2] + 0.003*nitrogen[2] + 0.31*oxygen[2] + 0.001*sulfur[2]) \
                          + w_lp*(0.10*hydrogen[2] + 0.80*carbon[2] + 0.10*oxygen[2])
        #
        majors_data = np.array([["H", hydrogen[1], w_ch*0.06*hydrogen[2] + w_lg*0.06*hydrogen[2] + w_lp*0.10*hydrogen[2], hydrogen[2]],
                                ["C", carbon[1], w_ch*0.44*carbon[2] + w_lg*0.63*carbon[2] + w_lp*0.80*carbon[2], carbon[2]],
                                ["N", nitrogen[1], w_lp*0.003*nitrogen[2], nitrogen[2]],
                                ["O", oxygen[1], w_ch*0.50*oxygen[2] + w_lg*0.31*oxygen[2] + w_lp*0.10*oxygen[2], oxygen[2]],
                                ["S", sulfur[1], w_lp*0.001*sulfur[2], sulfur[2]]], dtype=object)
        #
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = w_ch*carbohydrates[2] + w_lg*lignin[2] + w_lp*lipid[2]
        # Bulk modulus
        K = (w_ch*carbohydrates[3][0] + w_lg*lignin[3][0] + w_lp*lipid[3][0])*10**9
        # Shear modulus
        G = (w_ch*carbohydrates[3][1] + w_lg*lignin[3][1] + w_lp*lipid[3][1])*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        GR = 0
        # Photoelectricity
        PE = round(w_ch*carbohydrates[5][1] + w_lg*lignin[5][1] + w_lp*lipid[5][1], 3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(GR, 2), round(PE, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_rhodochrosite(self):   # Mn(CO3)
        #
        name = "Rdc"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["C", "O", "Mn"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Ca", "Mg", "Zn", "Co", "Cd"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = manganese[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.777, 15.67], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 99*10**9
        # Shear modulus
        G = 44*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = name
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_aragonite(self):   # Ca(CO3)
        #
        name = "Arg"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Ca"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Sr", "Pb", "Zn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = calcium[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.959, 7.968, 5.741], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 69*10**9
        # Shear modulus
        G = 33*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = name
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_cerussite(self):   # Pb(CO3)
        #
        name = "Cer"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["C", "O", "Pb"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Pb", lead[1], 1, lead[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = lead[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.195, 8.436, 6.152], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 56*10**9
        # Shear modulus
        G = 21*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = name
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_ankerite(self):   # Ca(Fe,Mg)(CO3)2
        #
        name = "Ank"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Mg", "Fe"]
        #
        x = round(rd.uniform(0, 1.0), 4)
        #
        majors_data = np.array([["C", carbon[1], 2, carbon[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Mg", magnesium[1], (1-x), magnesium[2]], ["Ca", calcium[1], 1, calcium[2]],
                                ["Fe", iron[1], x, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = calcium[2] + (x*iron[2] + (1-x)*magnesium[2]) + 2*(carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[4.83, 16.167], [], "trigonal"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 3
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        dataV_Mg = CrystalPhysics([[4.83, 16.167], [], "trigonal"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 3
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Fe.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        V = x*V_Fe + (1-x)*V_Mg
        V_m = x*V_m_Fe + (1-x)*V_m_Mg
        rho = x*rho_Fe + (1-x)*rho_Mg
        rho_e = x*rho_e_Fe + (1-x)*rho_e_Mg
        # Bulk modulus
        K_Fe = 73*10**9
        K_Mg = 89*10**9
        K = x*K_Fe + (1-x)*K_Mg
        # Shear modulus
        G_Fe = 32*10**9
        G_Mg = 44*10**9
        G = x*G_Fe + (1-x)*G_Mg
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = name
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_azurite(self):   # Cu3(CO3)2(OH)2
        #
        name = "Az"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["H", "C", "O", "Cu"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["C", carbon[1], 2, carbon[2]],
                                ["O", oxygen[1], 8, oxygen[2]], ["Cu", copper[1], 3, copper[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = 3*copper[2] + 2*(carbon[2] + 3*oxygen[2]) + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.008, 5.844, 10.336], [92.333], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 112.35*10**9
        # Shear modulus
        G = 49.33*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = name
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_malachite(self):   # Cu2(CO3)(OH)2
        #
        name = "Mal"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["H", "C", "O", "Cu"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["C", carbon[1], 2, carbon[2]],
                                ["O", oxygen[1], 8, oxygen[2]], ["Cu", copper[1], 3, copper[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Zn", "Co", "Ni"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = 2*copper[2] + (carbon[2] + 3*oxygen[2]) + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.502, 11.974, 3.24], [98.75], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 118.77*10**9
        # Shear modulus
        G = 50.06*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = name
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_ikaite(self):   # CaCO3 * 6*H2O
        #
        name = "Ika"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "C", "O", "Ca"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["C", carbon[1], 2, carbon[2]],
                                ["O", oxygen[1], 8, oxygen[2]], ["Ca", calcium[1], 3, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = calcium[2] + (carbon[2] + 3*oxygen[2]) + 6*(2*hydrogen[2] + oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.87, 8.23, 11.02], [110.2], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = np.mean([30.9524, 40.3962])*10**9
        # Shear modulus
        G = np.mean([16.4415, 18.714])*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_smithsonite(self):  # ZnCO3
        #
        name = "Smt"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["C", "O", "Zn"]
        #
        majors_data = np.array(
            [["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]], ["Zn", zinc[1], 1, zinc[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Co", "Cu", "Mn", "Ca", "Cd", "Mg", "In"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = zinc[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.653, 15.028], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 117*10**9
        # Shear modulus
        G = 49*10**9
        # Young's modulus
        E = (9 * K * G) / (3 * K + G)
        # Poisson's ratio
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
        # vP/vS
        vPvS = ((K + 4 / 3 * G) / G) ** 0.5
        # P-wave velocity
        vP = ((K + 4 / 3 * G) / rho) ** 0.5
        # S-wave velocity
        vS = (G / rho) ** 0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe * rho_e * 10 ** (-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G * 10 ** (-9), 4)
        results["K"] = round(K * 10 ** (-9), 4)
        results["E"] = round(E * 10 ** (-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
#
class CustomCarbonates:
    #
    def __init__(self, fluid, actualThickness, output_type=False, porosity=None):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.output_type = output_type
        self.porosity = porosity
    #
    def create_custom_rock_01(self, amounts=None):
        #
        self.amounts = amounts
        #
        # Mineralogy + Fluids
        org = Carbonates(impurity="pure", dict=True).create_organic_matter()
        quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        calcite = Carbonates(impurity="pure", dict=True).create_calcite()
        dolomite = Carbonates(impurity="pure", dict=True).create_dolomite()
        siderite = Carbonates(impurity="pure", dict=True).create_siderite()
        anhydrite = Sulfates(data_type=True).create_anhydrite()
        gypsum = Sulfates(data_type=True).create_gypsum()
        pyrite = Sulfides(impurity="pure", dict=True).create_pyrite()
        illite = Pyllosilicates(impurity="pure", dict=True).create_illite()
        #
        mineralogy = [org, quartz, alkalifeldspar, plagioclase, calcite, dolomite, siderite, anhydrite, gypsum, pyrite,
                      illite]
        #
        water = fluids.Water.water("")
        #
        data = []
        results = {}
        #
        cond = False
        composition = []
        while cond == False:
            if self.amounts == None:
                w_org = round(0.1/100, 4)
                w_qz = round(1.3/100, 4)
                w_kfs = round(0/100, 4)
                w_pl = round(0/100, 4)
                w_cal = round(8.1/100, 4)
                w_dol = round(87/100, 4)
                w_sd = round(0/100, 4)
                w_anh = round(0/100, 4)
                w_gyp = round(0/100, 4)
                w_py = round(0.1/100, 4)
                w_ilt = round(1 - w_org - w_qz - w_kfs - w_pl - w_cal - w_dol - w_sd - w_anh - w_gyp - w_py, 4)
            elif type(self.amounts) is list:
                w_org = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_qz = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_dol = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_sd = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_anh = round(abs(np.random.normal(self.amounts[7], 0.025)), 4)
                w_gyp = round(abs(np.random.normal(self.amounts[8], 0.025)), 4)
                w_py = round(abs(np.random.normal(self.amounts[9], 0.025)), 4)
                w_ilt = round(1 - w_org - w_qz - w_kfs - w_pl - w_cal - w_dol - w_sd - w_anh - w_gyp - w_py, 4)
            #
            if w_org >= 0.0 and w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_cal >= 0.0 and w_dol >= 0.0 \
                    and w_sd >= 0.0 and w_anh >= 0.0 and w_gyp >= 0.0 and w_py >= 0.0 and w_ilt >= 0.0:
                sumMin = round(w_org + w_qz + w_kfs + w_pl + w_cal + w_dol + w_sd + w_anh + w_gyp + w_py + w_ilt, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_org*org["chemistry"]["H"] + w_gyp*gypsum["chemistry"]["H"] + w_ilt*illite["chemistry"]["H"], 4)
            w_C = round(w_org*org["chemistry"]["C"] + w_cal*calcite["chemistry"]["C"] + w_dol*dolomite["chemistry"]["C"] + w_sd*siderite["chemistry"]["C"], 4)
            w_N = round(w_org*org["chemistry"]["N"], 4)
            w_Na = round(w_kfs*alkalifeldspar["chemistry"]["Na"] + w_pl*plagioclase["chemistry"]["Na"], 4)
            w_Mg = round(w_dol*dolomite["chemistry"]["Mg"] + w_ilt*illite["chemistry"]["Mg"], 4)
            w_Al = round(w_kfs*alkalifeldspar["chemistry"]["Al"] + w_pl*plagioclase["chemistry"]["Al"] + w_ilt*illite["chemistry"]["Al"], 4)
            w_Si = round(w_qz*quartz["chemistry"]["Si"] + w_kfs*alkalifeldspar["chemistry"]["Si"] + w_pl*plagioclase["chemistry"]["Si"] + w_ilt*illite["chemistry"]["Si"], 4)
            w_S = round(w_org*org["chemistry"]["S"] + w_py*pyrite["chemistry"]["S"], 4)
            w_K = round(w_kfs*alkalifeldspar["chemistry"]["K"] + w_ilt*illite["chemistry"]["K"], 4)
            w_Ca = round(w_pl*plagioclase["chemistry"]["Ca"] + w_cal*calcite["chemistry"]["Ca"] + w_dol*dolomite["chemistry"]["Ca"] + w_anh*anhydrite["chemistry"]["Ca"] + w_gyp*gypsum["chemistry"]["Ca"], 4)
            w_Fe = round(w_sd*siderite["chemistry"]["Fe"] + w_py*pyrite["chemistry"]["Fe"] + w_ilt*illite["chemistry"]["Fe"], 4)
            w_O = round(1 - w_H - w_C - w_N - w_Na - w_Mg - w_Al - w_Si - w_S - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_C + w_N + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Org", "Qz", "Kfs", "Pl", "Cal", "Dol", "Sd", "Anh", "Gyp", "Py", "Ilt"]))
                concentrations = [w_H, w_C, w_N, w_O, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe]
                amounts = [w_org, w_qz, w_kfs, w_pl, w_cal, w_dol, w_sd, w_anh, w_gyp, w_py, w_ilt]
            else:
                cond = False
        #
        element_list = ["H", "C", "N", "O", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Fe"]
        mineral_list = ["Org", "Qz", "Kfs", "Pl", "Cal", "Dol", "Sd", "Anh", "Gyp", "Py", "Ilt"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_org*org["rho"] + w_qz*quartz["rho"] + w_kfs*alkalifeldspar["rho"] + w_pl*plagioclase["rho"]
                    + w_cal*calcite["rho"] + w_dol*dolomite["rho"] + w_sd*siderite["rho"] + w_anh*anhydrite["rho"]
                    + w_gyp*gypsum["rho"] + w_py*pyrite["rho"] + w_ilt*illite["rho"]) / 1000
        rhoSolid = 0.975*rhoSolid
        X = [w_org, w_qz, w_kfs, w_pl, w_cal, w_dol, w_sd, w_anh, w_gyp, w_py, w_ilt]
        K_list = [mineralogy[i]["K"] for i in range(len(mineralogy))]
        G_list = [mineralogy[i]["G"] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = 0.35*K_geo
        G_solid = 0.275*G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
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
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_org*org["GR"] + w_qz*quartz["GR"] + w_kfs*alkalifeldspar["GR"] + w_pl*plagioclase["GR"] \
             + w_cal*calcite["GR"] + w_dol*dolomite["GR"] + w_sd*siderite["GR"] + w_anh*anhydrite["GR"] \
             + w_gyp*gypsum["GR"] + w_py*pyrite["GR"] + w_ilt*illite["GR"]
        PE = w_org*org["PE"] + w_qz*quartz["PE"] + w_kfs*alkalifeldspar["PE"] + w_pl*plagioclase["PE"] \
             + w_cal*calcite["PE"] + w_dol*dolomite["PE"] + w_sd*siderite["PE"] + w_anh*anhydrite["PE"] \
             + w_gyp*gypsum["PE"] + w_py*pyrite["PE"] + w_ilt*illite["PE"]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_org*org["nu"] + w_qz*quartz["nu"] + w_kfs*alkalifeldspar["nu"] + w_pl*plagioclase["nu"] \
                                + w_cal*calcite["nu"] + w_dol*dolomite["nu"] + w_sd*siderite["nu"] + w_anh*anhydrite["nu"] \
                                + w_gyp*gypsum["nu"] + w_py*pyrite["nu"] + w_ilt*illite["nu"]
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
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append("water")
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        #
        if dict == False:
            return data
        else:
            return results