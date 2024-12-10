#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		siliciclastics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		10.12.2024

#-----------------------------------------------

## MODULES
import datetime
import sys

import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals, oxides, fluids
from modules.geophysics import Elasticity as elast
from modules import oxides, carbonates, silicates
from modules.chemistry import PeriodicSystem, OxideCompounds
from modules.oxides import Oxides
from modules.carbonates import Carbonates
from modules.silicates import Phyllosilicates
from modules.silicates import Tectosilicates
from modules.sulfides import Sulfides
from modules.organics import Organics
from modules.phosphates import Phosphates
from modules.fluids import Water, Hydrocarbons
from modules.petrophysics import SeismicVelocities

class Geophysics:
    #
    def __init__(self, var_data):
        self.var_data = var_data
    #
    def calculate_seismic_velocity(self):
        results = {"vP": 0, "vS": 0}
        var_porosity = self.var_data["porosity"]
        #
        for var_type in results.keys():
            inv_v = 0
            for key, velocity in self.var_data["v"].items():
                fraction = self.var_data["Phi"][key]
                inv_v += fraction/velocity[var_type]
                #inv_v += fraction/(velocity[var_type]**(1 - var_porosity))
                #inv_v += ((1 - var_porosity)*fraction)/velocity[var_type]
            results[var_type] = round((1/inv_v)**(1 - var_porosity), 3)
            #results[var_type] = round((1/inv_v), 3)
        #
        return results

class Soil:
    #
    def __init__(self):
        pass
    #
    def create_simple_soil(self, w_C=None, amounts=None, grainsize_list=False):
        self.w_C = w_C
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        illite = minerals.phyllosilicates.illite("")
        kaolinite = minerals.phyllosilicates.kaolinite("")
        organic = minerals.natives.organic_matter("")
        #
        mineralogy = [quartz, illite, kaolinite, organic]
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        air = fluids.Gas.air("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_C == None and self.amounts == None:
                w_org = 0.05
                w_qz = round(abs(rd.uniform(0.5, 0.95)), 4)
                w_ilt = round(abs(rd.uniform(0.0, (1-w_org-w_qz))), 4)
                w_kln = round(abs(1-w_qz-w_ilt-w_org), 4)
            elif self.w_C != None:
                w_org = round(self.w_C, 4)
                w_mineral = round(1-w_org, 4)
                w_qz = round(abs(w_mineral*rd.uniform(0.25, 1)), 4)
                w_ilt = round(abs(w_mineral*rd.uniform(0, (1-w_qz))), 4)
                w_kln = round(abs(w_mineral*(1-w_qz-w_ilt)), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_org = round(1-w_qz-w_ilt-w_kln, 4)
            #
            if 0.0 <= w_qz <= 1.0 and 0.0 <= w_ilt <= 1.0 and 0.0 <= w_kln <= 1.0 and 0.0 <= w_org <= 1.0:
                sumMin = round(w_qz + w_ilt + w_kln + w_org, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ilt*illite[6][0] + w_kln*kaolinite[6][0], 4)
            w_C = round(w_org, 4)
            w_O = round(w_qz*quartz[6][0] + w_ilt*illite[6][1] + w_kln*kaolinite[6][1], 4)
            w_Mg = round(w_ilt*illite[6][2], 4)
            w_Al = round(w_ilt*illite[6][3] + w_kln*kaolinite[6][2], 4)
            w_Si = round(w_qz*quartz[6][1] + w_ilt*illite[6][4] + w_kln*kaolinite[6][3], 4)
            w_K = round(w_ilt*illite[6][5], 4)
            w_Fe = round(w_ilt*illite[6][6], 4)
            sumConc = w_H + w_C + w_O + w_Mg + w_Al + w_Si + w_K + w_Fe
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Qz", "Ilt", "Kln", "Org"]))
                concentrations = [w_H, w_C, w_O, w_Mg, w_Al, w_Si, w_K, w_Fe]
                amounts = [w_qz, w_ilt, w_kln, w_org]
            else:
                cond = False
        #
        grainsize = []
        n_Grains = 100
        if w_qz > 0:
            grainsize.append([rd.randint(2, 2000) for i in range(int(round(w_qz*n_Grains, 0)))])
        if w_ilt > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_ilt*n_Grains, 0)))], 2)))
        if w_kln > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_kln*n_Grains, 0)))], 2)))
        if w_org > 0:
            grainsize.append([rd.randint(2, 63) for i in range(int(round(w_org*n_Grains, 0)))])
        #
        rhoSolid = (w_qz*quartz[2] + w_ilt*illite[2] + w_kln*kaolinite[2] + w_org*organic[2])/1000
        X = [w_qz, w_ilt, w_kln, w_org]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        #vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        #vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        vP_solid = w_qz*quartz[4][0] + w_ilt*illite[4][0] + w_kln*kaolinite[4][0] + w_org*organic[4][0]
        vS_solid = w_qz*quartz[4][1] + w_ilt*illite[4][1] + w_kln*kaolinite[4][1] + w_org*organic[4][1]
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = rd.uniform(0.5, 0.65)
        #
        rho = (1 - phi) * rhoSolid + phi * (0.5*water[2]+0.5*air[2]/1000) / 1000
        vP = ((1-phi)*vP_solid + phi*(0.5*water[4][0] + 0.5*air[4][0]))/3
        vS = ((1 - phi) * vS_solid)/3
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - (0.5*water[2]+0.5*air[2]/1000) / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_ilt*illite[5][0] + w_kln*kaolinite[5][0] + w_org*organic[5][0]
        PE = w_qz*quartz[5][1] + w_ilt*illite[5][1] + w_kln*kaolinite[5][1] + w_org*organic[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_ilt*illite[3][3] + w_kln*kaolinite[3][3] + w_org*organic[3][3]
        #
        data.append(composition)
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6), round(air[2]/1000, 3)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2), round(air[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append(["water", "air"])
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        if grainsize_list == True:
            data.append(grainsize)
        #
        return data
    #
    def create_simple_sand(self, w_C=None, amounts=None, grainsize_list=False):
        self.w_C = w_C
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        illite = minerals.phyllosilicates.illite("")
        kaolinite = minerals.phyllosilicates.kaolinite("")
        organic = minerals.natives.organic_matter("")
        #
        mineralogy = [quartz, illite, kaolinite, organic]
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        air = fluids.Gas.air("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_C == None and self.amounts == None:
                w_org = 0.025
                w_qz = round(abs(rd.uniform(0.85, 0.975)), 4)
                w_ilt = round(abs(rd.uniform(0.0, (1-w_org-w_qz))), 4)
                w_kln = round(abs(1-w_qz-w_ilt-w_org), 4)
            elif self.w_C != None:
                w_org = round(self.w_C, 4)
                w_mineral = round(1-w_org, 4)
                w_qz = round(abs(w_mineral*rd.uniform(0.25, 1)), 4)
                w_ilt = round(abs(w_mineral*rd.uniform(0, (1-w_qz))), 4)
                w_kln = round(abs(w_mineral*(1-w_qz-w_ilt)), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_org = round(1-w_qz-w_ilt-w_kln, 4)
            #
            if 0.0 <= w_qz <= 1.0 and 0.0 <= w_ilt <= 1.0 and 0.0 <= w_kln <= 1.0 and 0.0 <= w_org <= 1.0:
                sumMin = round(w_qz + w_ilt + w_kln + w_org, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ilt*illite[6][0] + w_kln*kaolinite[6][0], 4)
            w_C = round(w_org, 4)
            w_O = round(w_qz*quartz[6][0] + w_ilt*illite[6][1] + w_kln*kaolinite[6][1], 4)
            w_Mg = round(w_ilt*illite[6][2], 4)
            w_Al = round(w_ilt*illite[6][3] + w_kln*kaolinite[6][2], 4)
            w_Si = round(w_qz*quartz[6][1] + w_ilt*illite[6][4] + w_kln*kaolinite[6][3], 4)
            w_K = round(w_ilt*illite[6][5], 4)
            w_Fe = round(w_ilt*illite[6][6], 4)
            sumConc = w_H + w_C + w_O + w_Mg + w_Al + w_Si + w_K + w_Fe
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Qz", "Ilt", "Kln", "Org"]))
                concentrations = [w_H, w_C, w_O, w_Mg, w_Al, w_Si, w_K, w_Fe]
                amounts = [w_qz, w_ilt, w_kln, w_org]
            else:
                cond = False
        #
        grainsize = []
        n_Grains = 100
        if w_qz > 0:
            grainsize.append([rd.randint(2, 2000) for i in range(int(round(w_qz*n_Grains, 0)))])
        if w_ilt > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_ilt*n_Grains, 0)))], 2)))
        if w_kln > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_kln*n_Grains, 0)))], 2)))
        if w_org > 0:
            grainsize.append([rd.randint(2, 63) for i in range(int(round(w_org*n_Grains, 0)))])
        #
        rhoSolid = (w_qz*quartz[2] + w_ilt*illite[2] + w_kln*kaolinite[2] + w_org*organic[2])/1000
        X = [w_qz, w_ilt, w_kln, w_org]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        #vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        #vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        vP_solid = w_qz*quartz[4][0] + w_ilt*illite[4][0] + w_kln*kaolinite[4][0] + w_org*organic[4][0]
        vS_solid = w_qz*quartz[4][1] + w_ilt*illite[4][1] + w_kln*kaolinite[4][1] + w_org*organic[4][1]
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = rd.uniform(0.45, 0.55)
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = ((1-phi)*vP_solid + phi*water[4][0])/3
        vS = ((1 - phi) * vS_solid)/3
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2]/1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_ilt*illite[5][0] + w_kln*kaolinite[5][0] + w_org*organic[5][0]
        PE = w_qz*quartz[5][1] + w_ilt*illite[5][1] + w_kln*kaolinite[5][1] + w_org*organic[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_ilt*illite[3][3] + w_kln*kaolinite[3][3] + w_org*organic[3][3]
        #
        data.append(composition)
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append("water")
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        if grainsize_list == True:
            data.append(grainsize)
        #
        return data

#########################
## SILICICLASTIC ROCKS ##
#########################
class SiliciclasticRocks:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        self.data_pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()
        self.data_uraninite = Oxides(impurity="pure", data_type=True).create_uraninite()
        self.data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        #
        self.data_water = Water.water("")
        self.data_oil = Hydrocarbons.oil("")
        self.data_gas = Hydrocarbons.natural_gas("")
    #
    def create_sandstone(self, rock="Sandstone", number=1, composition=None, classification="Sandstone", porosity=None):
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
        results_container["phi"] = []
        results_container["phi_true"] = []
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

        n = 0
        while n < number:
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                          "Kln": self.data_kaolinite, "Hem": self.data_hematite}
            minerals_list = list(mineralogy.keys())

            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []

            condition = False
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}

                if composition != None:
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_kln = composition["Kln"]
                    phi_hem = composition["Hem"]

                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Hem"] = phi_hem

                else:
                    condition_2 = False
                    while condition_2 == False:
                        if classification == "Sandstone":
                            qz_limits = [0.7, 1.0]
                            kfs_limits = [0.0, 0.25]
                            pl_limits = [0.0, 0.25]
                            kln_limits = [0.0, 0.25]
                            hem_limits = [0.0, 0.05]

                        phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                        phi_kfs = round(rd.uniform(kfs_limits[0], (1 - phi_qz)), 4)
                        phi_pl = round(rd.uniform(pl_limits[0], (1 - phi_qz - phi_kfs)), 4)
                        phi_kln = round(rd.uniform(kln_limits[0], (1 - phi_qz - phi_kfs - phi_pl)), 4)
                        phi_hem = round(1 - phi_qz - phi_kfs - phi_pl - phi_kln, 4)

                        phi_total = phi_qz + phi_kfs + phi_pl + phi_kln + phi_hem

                        if np.isclose(phi_total, 1.0000) == True:
                            if qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                    and kln_limits[0] <= phi_kln <= kln_limits[1] \
                                    and hem_limits[0] <= phi_hem <= hem_limits[1]:
                                condition_2 = True

                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Hem"] = phi_hem

                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                    for element, value in mineralogy[key]["chemistry"].items():
                        if element not in elements_list:
                            elements_list.append(element)
                            w_elements[element] = 0.0

                if elements_list[0] not in results_container["chemistry"]:
                    for element in elements_list:
                        results_container["chemistry"][element] = []

                rho_s = round(rho_s, 3)
                for key, value in phi_minerals.items():
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4

                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)

                if self.fluid == "water":
                    data_fluid = self.data_water

                old_index = elements_list.index("O")
                elements_list += [elements_list.pop(old_index)]

                w_elements_total = 0.0
                for element in elements_list:
                    if element != "O":
                        for mineral, w_mineral in w_minerals.items():
                            if element in mineralogy[mineral]["chemistry"]:
                                if element == "U":
                                    n_digits = 4
                                else:
                                    n_digits = 4

                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        w_elements[element] = round(w_elements[element], 4)

                total_w_minerals = round(sum(w_minerals.values()), 4)
                total_w_elements = round(sum(w_elements.values()), 4)
                if total_w_minerals == 1.0 and total_w_elements == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    condition = True

            velocity_solid = {"vP": 0, "vS": 0}
            gamma_ray = 0.0
            photoelectricity = 0.0
            for key, value in phi_minerals.items():
                velocity_solid["vP"] += phi_minerals[key]*mineralogy[key]["vP"]
                velocity_solid["vS"] += phi_minerals[key]*mineralogy[key]["vS"]
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]

            ## Bulk Density, Porosity, Seismic Velocities
            rho_solid = round(rho_s, 3)
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                rho_limits=[1800, 2800], vP_limits=[2800, 4800], vS_limits=[1500, 2500], delta=0.05,
                porosity=porosity)
            phi_neutron = round((1900/rho)*0.15, 4)
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

            list_oxides = ["H2O", "Na2O", "Al2O3", "SiO2", "K2O", "CaO", "Fe2O3"]
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
            results_container["phi"].append(phi_neutron)
            results_container["phi_true"].append(var_porosity)
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
            n += 1

        return results_container

    def create_siltstone(self, rock="Siltstone", number=1, composition=None, classification="Siltstone", porosity=None):
        results_container = {}
        results_container["rock"] = rock
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
        results_container["phi"] = []
        results_container["phi_true"] = []
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

        n = 0
        while n < number:
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_apatite = Phosphates(data_type=True).create_aptite()
            mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                          "Kln": self.data_kaolinite, "Ap": data_apatite}
            minerals_list = list(mineralogy.keys())

            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []

            condition = False
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}

                if composition != None:
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_kln = composition["Kln"]
                    phi_ap = composition["Ap"]

                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ap"] = phi_ap

                else:
                    condition_2 = False
                    while condition_2 == False:
                        if classification == "Siltstone":
                            qz_limits = [0.7, 1.0]
                            kfs_limits = [0.0, 0.25]
                            pl_limits = [0.0, 0.25]
                            kln_limits = [0.0, 0.25]
                            ap_limits = [0.0, 0.025]

                        phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                        phi_kfs = round(rd.uniform(kfs_limits[0], (1 - phi_qz)), 4)
                        phi_pl = round(rd.uniform(pl_limits[0], (1 - phi_qz - phi_kfs)), 4)
                        phi_kln = round(rd.uniform(kln_limits[0], (1 - phi_qz - phi_kfs - phi_pl)), 4)
                        phi_ap = round(1 - phi_qz - phi_kfs - phi_pl - phi_kln, 4)

                        phi_total = phi_qz + phi_kfs + phi_pl + phi_kln + phi_ap

                        if np.isclose(phi_total, 1.0000) == True:
                            if qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                    and kln_limits[0] <= phi_kln <= kln_limits[1] \
                                    and ap_limits[0] <= phi_ap <= ap_limits[1]:
                                condition_2 = True

                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ap"] = phi_ap

                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                    for element, value in mineralogy[key]["chemistry"].items():
                        if element not in elements_list:
                            elements_list.append(element)
                            w_elements[element] = 0.0

                if elements_list[0] not in results_container["chemistry"]:
                    for element in elements_list:
                        results_container["chemistry"][element] = []

                rho_s = round(rho_s, 3)
                for key, value in phi_minerals.items():
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4

                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)

                if self.fluid == "water":
                    data_fluid = self.data_water

                old_index = elements_list.index("O")
                elements_list += [elements_list.pop(old_index)]

                w_elements_total = 0.0
                for element in elements_list:
                    if element != "O":
                        for mineral, w_mineral in w_minerals.items():
                            if element in mineralogy[mineral]["chemistry"]:
                                if element == "U":
                                    n_digits = 4
                                else:
                                    n_digits = 4

                                value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                w_elements[element] += value
                                w_elements_total += value
                                w_elements[element] = round(w_elements[element], n_digits)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        w_elements[element] = round(w_elements[element], 4)

                total_w_minerals = round(sum(w_minerals.values()), 4)
                total_w_elements = round(sum(w_elements.values()), 4)
                if total_w_minerals == 1.0 and total_w_elements == 1.0:
                    for key, value in w_minerals.items():
                        w_minerals[key] = abs(value)
                    for key, value in w_elements.items():
                        w_elements[key] = abs(value)
                    condition = True

            velocity_solid = {"vP": 0, "vS": 0}
            gamma_ray = 0.0
            photoelectricity = 0.0
            for key, value in phi_minerals.items():
                velocity_solid["vP"] += phi_minerals[key]*mineralogy[key]["vP"]
                velocity_solid["vS"] += phi_minerals[key]*mineralogy[key]["vS"]
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]

            ## Bulk Density, Porosity, Seismic Velocities
            rho_solid = round(rho_s, 3)
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                rho_limits=[1800, 2800], vP_limits=[2800, 4800], vS_limits=[1500, 2500], delta=0.05,
                porosity=porosity)
            phi_neutron = round((1900/rho)*0.15, 4)
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

            list_oxides = ["H2O", "F", "Na2O", "Al2O3", "SiO2", "P2O5", "Cl", "K2O", "CaO"]
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
            results_container["phi"].append(phi_neutron)
            results_container["phi_true"].append(var_porosity)
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
            n += 1

        return results_container

    def create_mudstone_alt(self, rock="Mudstone", number=1, composition=None, classification="Mudstone",
                            porosity=None):
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
        helper = [[], []]
        while n < number:
            data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite_simple()
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_organics = Organics(data_type=True).create_organic_matter()
            #
            mineralogy = {
                "Kln": self.data_kaolinite, "Ilt": data_illite, "Qz": self.data_quartz, "Kfs": data_alkalifeldspar,
                "Pl": data_plagioclase, "Org": data_organics, "Py": self.data_pyrite}
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
                    phi_kln = composition["Kln"]
                    phi_ilt = composition["Ilt"]
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_org = composition["Org"]
                    phi_py = composition["Py"]
                    #
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Org"] = phi_org
                    phi_minerals["Py"] = phi_py
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        if classification == "Mudstone":
                            kln_limits = [0.4, 0.5]
                            ilt_limits = [0.0, 0.1]
                            qz_limits = [0.1, 0.3]
                            kfs_limits = [0.0, 0.1]
                            pl_limits = [0.0, 0.1]
                            org_limits = [0.0, 0.1]
                            py_limits = [0.0, 0.05]
                        #
                        phi_kln = round(rd.uniform(kln_limits[0], kln_limits[1]), 4)
                        #
                        condition_ilt = False
                        while condition_ilt == False:
                            phi_ilt = round(rd.uniform(ilt_limits[0], (1 - phi_kln)), 4)
                            if ilt_limits[0] <= phi_ilt <= ilt_limits[1]:
                                condition_ilt = True
                        #
                        condition_qz = False
                        while condition_qz == False:
                            phi_qz = round(rd.uniform(qz_limits[0], (1 - phi_kln - phi_ilt)), 4)
                            if qz_limits[0] <= phi_qz <= qz_limits[1]:
                                condition_qz = True
                        #
                        condition_kfs = False
                        while condition_kfs == False:
                            phi_kfs = round(rd.uniform(kfs_limits[0], (1 - phi_kln - phi_ilt - phi_qz)), 4)
                            if kfs_limits[0] <= phi_kfs <= kfs_limits[1]:
                                condition_kfs = True
                        #
                        condition_pl = False
                        while condition_pl == False:
                            phi_pl = round(rd.uniform(pl_limits[0], (1 - phi_kln - phi_ilt - phi_qz - phi_kfs)), 4)
                            if pl_limits[0] <= phi_pl <= pl_limits[1]:
                                condition_pl = True
                        #
                        condition_org = False
                        while condition_org == False:
                            phi_org = round(rd.uniform(
                                org_limits[0], (1 - phi_kln - phi_ilt - phi_qz - phi_kfs - phi_pl)), 4)
                            if org_limits[0] <= phi_org <= org_limits[1]:
                                condition_org = True
                        #
                        phi_py = round(1 - phi_kln - phi_ilt - phi_qz - phi_kfs - phi_pl - phi_org, 4)
                        #
                        phi_total = phi_kln + phi_ilt + phi_qz + phi_kfs + phi_pl + phi_org + phi_py
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if kln_limits[0] <= phi_kln <= kln_limits[1] \
                                    and ilt_limits[0] <= phi_ilt <= ilt_limits[1] \
                                    and qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                    and org_limits[0] <= phi_org <= org_limits[1] \
                                    and py_limits[0] <= phi_py <= py_limits[1]:
                                condition_2 = True
                    #
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Org"] = phi_org
                    phi_minerals["Py"] = phi_py
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

            velocity_solid = {"vP": 0, "vS": 0}
            gamma_ray = 0.0
            photoelectricity = 0.0
            for key, value in phi_minerals.items():
                velocity_solid["vP"] += phi_minerals[key]*mineralogy[key]["vP"]
                velocity_solid["vS"] += phi_minerals[key]*mineralogy[key]["vS"]
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
            #
            ## Bulk Density, Porosity, Seismic Velocities
            rho_solid = round(rho_s, 3)
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                rho_limits=[1800, 2800], vP_limits=[3000, 4500], vS_limits=[2000, 2600], delta=0.05,
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

            list_oxides = ["H2O", "CO2", "N2O5", "Na2O", "Al2O3", "SiO2", "SO3", "K2O", "CaO", "Fe2O3"]
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

    def create_conglomerate(self, number, porosity=None):
        #
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        assemblage = [self.data_quartz, data_alkalifeldspar, data_plagioclase, data_biotite, self.data_calcite,
                      self.data_hematite]
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
                    if value >= 0.0 and 0.5 <= value <= 0.8:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Cal":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.1, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.1 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Bt":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
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
                phi_helper = round(rd.uniform(0.0, 0.2), 4)
            else:
                phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
            #
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
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
                if element in data_biotite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Bt"][n] * data_biotite["chemistry"][element], 4)
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
                    elif mineral == "Bt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_biotite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_biotite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_biotite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_biotite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_biotite["PE"], 3)
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

                amounts = []
                for key, value in amounts_chemistry.items():
                    chem_data = PeriodicSystem(name=key).get_data()
                    amounts.append([key, chem_data[1], value[-1]])

                list_oxides = ["H2O", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "K2O", "CaO", "Fe2O3"]
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

        results = {}
        results["rock"] = "Conglomerate"
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
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            for compound, value in amounts_compounds.items():
                single_amounts_compounds[compound] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["compounds"] = single_amounts_compounds
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

    def create_shale_alt(self, number=1, composition=None, porosity=None):
        results_container = {}
        results_container["rock"] = "Shale"
        results_container["mineralogy"] = {}
        results_container["chemistry"] = {}
        results_container["compounds"] = {}
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
            data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
            data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite_simple()
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_organics = Organics(data_type=True).create_organic_matter()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            mineralogy = {"Ilt": data_illite, "Mnt": data_montmorillonite, "Chl": data_chlorite,
                          "Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Bt": data_biotite, "Org": data_organics,
                          "Urn": self.data_uraninite, "Py": self.data_pyrite}
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
                    phi_ilt = composition["Ilt"]
                    phi_mnt = composition["Mnt"]
                    phi_chl = composition["Chl"]
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_bt = composition["Bt"]
                    phi_org = composition["Org"]
                    phi_urn = composition["Urn"]
                    phi_py = composition["Py"]
                    #
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Mnt"] = phi_mnt
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Bt"] = phi_bt
                    phi_minerals["Org"] = phi_org
                    phi_minerals["Urn"] = phi_urn
                    phi_minerals["Py"] = phi_py
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        w_misc = round(rd.uniform(0.0, 0.2), 4)
                        w_silic = round(rd.uniform(0.0, 0.4), 4)
                        w_clay = round(1 - w_misc - w_silic, 4)
                        #
                        ## Others
                        upper_limit_urn = 0.000025
                        phi_urn = round(rd.uniform(0.000001, upper_limit_urn), 6)
                        phi_org = round(w_misc*rd.uniform(0.0, (1 - phi_urn)), 6)
                        phi_py = round(w_misc*rd.uniform(0.0, (1 - phi_urn - phi_org)), 6)
                        phi_bt = round(w_misc - phi_urn - phi_org - phi_py, 6)
                        #
                        ## Siliciclastics
                        phi_qz = round(w_silic*rd.uniform(0.2, 0.8), 6)
                        phi_kfs = round(w_silic - phi_qz, 6)
                        #
                        ## Clays
                        phi_ilt = round(w_clay*rd.uniform(0.5, 1.0), 6)
                        phi_mnt = round(w_clay*rd.uniform(0.0, (1 - phi_ilt)), 6)
                        phi_chl = round(w_clay - phi_ilt - phi_mnt, 6)
                        #
                        phi_total = phi_urn + phi_org + phi_bt + phi_py + phi_qz + phi_kfs + phi_ilt + phi_mnt + phi_chl
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.0 <= phi_urn <= upper_limit_urn and 0.0 <= phi_org <= w_misc \
                                    and 0.0 <= phi_bt <= w_misc and 0.0 <= phi_py <= 0.1 \
                                    and 0.0 <= phi_qz <= w_silic and 0.1 <= phi_kfs <= w_silic \
                                    and 0.0 <= phi_ilt <= w_clay and 0.0 <= phi_mnt <= w_clay \
                                    and 0.0 <= phi_chl <= w_clay:
                                condition_2 = True
                        #
                    phi_minerals["Ilt"] = phi_ilt
                    phi_minerals["Mnt"] = phi_mnt
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Bt"] = phi_bt
                    phi_minerals["Org"] = phi_org
                    phi_minerals["Urn"] = phi_urn
                    phi_minerals["Py"] = phi_py
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
                        n_digits = 6
                    else:
                        n_digits = 6
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
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
            rho_solid = 0
            gamma_ray = 0.0
            photoelectricity = 0.0
            for key, value in phi_minerals.items():
                rho_solid += phi_minerals[key] * mineralogy[key]["rho"]
                gamma_ray += phi_minerals[key] * mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key] * mineralogy[key]["PE"]
            #
            ## Bulk Density, Porosity, Seismic Velocities
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=data_fluid[2]).calculate_seismic_velocities(
                rho_limits=[1800, 2900], vP_limits=[2000, 5000], vS_limits=[1000, 2000], delta=0.05, porosity=porosity)
            ## Elastic Parameters
            bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                rho=rho, vP=vP, vS=vS)
            ## Gamma Ray
            gamma_ray = round(gamma_ray, 3)
            ## Photoelectricity
            photoelectricity = round(photoelectricity, 3)
            phi_neutron = (2400/rho)*0.39
            # Composition data
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            amounts = []
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_oxides = ["H2O", "CO2", "N2O5", "Na2O", "MgO", "Al2O3", "SiO2", "SO3", "K2O", "CaO", "Mn2O3", "Fe2O3",
                           "NiO", "UO2"]
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

            ## Results
            results_container["phi"].append(var_porosity)
            results_container["rho_s"].append(rho_solid)
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
    #
    def create_greywacke_huckenholz(self, rock="Greywacke", number=1, composition=None, enrichment_fsp="Pl",
                                    enrichment_mca="Bt", porosity=None):
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
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            mineralogy = {
                "Qz": self.data_quartz, "Pl": data_plagioclase, "Kfs": data_alkalifeldspar, "Bt": data_biotite,
                "Ms": data_muscovite, "Chl": data_chlorite, "Cal": self.data_calcite, "Py": self.data_pyrite}
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
                    phi_pl = composition["Pl"]
                    phi_kfs = composition["Kfs"]
                    phi_bt = composition["Bt"]
                    phi_ms = composition["Ms"]
                    phi_chl = composition["Chl"]
                    phi_cal = composition["Cal"]
                    phi_py = composition["Py"]
                    #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Bt"] = phi_bt
                    phi_minerals["Ms"] = phi_ms
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Py"] = phi_py
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        fsp_dominance = round(rd.uniform(0.75, 1.0), 2)
                        mca_dominance = round(rd.uniform(0.75, 1.0), 2)
                        #
                        if enrichment_fsp == "Pl":
                            qz_limits = [0.25, 0.55]
                            #
                            pl_limits = [round(fsp_dominance*0.25, 2), round(fsp_dominance*0.45, 2)]
                            kfs_limits = [round((1 - fsp_dominance)*0.25, 2), round((1 - fsp_dominance)*0.47, 2)]
                            #
                            if enrichment_mca == "Bt":
                                bt_limits = [round(mca_dominance*0.0, 2), round(mca_dominance*0.20, 2)]
                                ms_limits = [round((1 - mca_dominance)*0.0, 2), round((1 - mca_dominance)*0.21, 2)]
                            else:
                                ms_limits = [round(mca_dominance*0.0, 2), round(mca_dominance*0.20, 2)]
                                bt_limits = [round((1 - mca_dominance)*0.0, 2), round((1 - mca_dominance)*0.21, 2)]
                            #
                            chl_limits = [0.0, 0.25]
                            cal_limits = [0.0, 0.06]
                            py_limits = [0.0, 0.03]
                            #
                        else:
                            qz_limits = [0.25, 0.55]
                            #
                            kfs_limits = [round(fsp_dominance*0.25, 2), round(fsp_dominance*0.45, 2)]
                            pl_limits = [round((1 - fsp_dominance)*0.25, 2), round((1 - fsp_dominance)*0.47, 2)]
                            #
                            if enrichment_mca == "Bt":
                                bt_limits = [round(mca_dominance*0.0, 2), round(mca_dominance*0.20, 2)]
                                ms_limits = [round((1 - mca_dominance)*0.0, 2), round((1 - mca_dominance)*0.21, 2)]
                            else:
                                ms_limits = [round(mca_dominance*0.0, 2), round(mca_dominance*0.20, 2)]
                                bt_limits = [round((1 - mca_dominance)*0.0, 2), round((1 - mca_dominance)*0.21, 2)]
                            #
                            chl_limits = [0.0, 0.22]
                            cal_limits = [0.0, 0.05]
                            py_limits = [0.0, 0.03]
                        #
                        phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                        phi_pl = round(rd.uniform(pl_limits[0], (1 - phi_qz)), 4)
                        phi_kfs = round(rd.uniform(kfs_limits[0], (1 - phi_qz - phi_pl)), 4)
                        phi_bt = round(rd.uniform(bt_limits[0], (1 - phi_qz - phi_pl - phi_kfs)), 4)
                        phi_ms = round(rd.uniform(ms_limits[0], (1 - phi_qz - phi_pl - phi_kfs - phi_bt)), 4)
                        phi_chl = round(rd.uniform(chl_limits[0], (1 - phi_qz - phi_pl - phi_kfs - phi_bt - phi_ms)), 4)
                        phi_cal = round(rd.uniform(cal_limits[0],
                                                   (1 - phi_qz - phi_pl - phi_kfs - phi_bt - phi_ms - phi_chl)), 4)
                        phi_py = round(1 - phi_qz - phi_pl - phi_kfs - phi_bt - phi_ms - phi_chl - phi_cal, 4)
                        #
                        phi_total = round(phi_qz + phi_pl + phi_kfs + phi_bt + phi_ms + phi_chl + phi_cal + phi_py, 4)
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                    and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                    and bt_limits[0] <= phi_bt <= bt_limits[1] \
                                    and ms_limits[0] <= phi_ms <= ms_limits[1] \
                                    and chl_limits[0] <= phi_chl <= chl_limits[1] \
                                    and cal_limits[0] <= phi_cal <= cal_limits[1] \
                                    and py_limits[0] <= phi_py <= py_limits[1]:
                                condition_2 = True
                        #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Bt"] = phi_bt
                    phi_minerals["Ms"] = phi_ms
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Py"] = phi_py
                #
                rho_s = 0
                velocities_minerals = {}
                for key, value in phi_minerals.items():
                    rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                    #
                    velocities_minerals[key] = {}
                    velocities_minerals[key]["vP"] = mineralogy[key]["vP"]
                    velocities_minerals[key]["vS"] = mineralogy[key]["vS"]
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
                #
                rho_solid = rho_s
                for key, value in phi_minerals.items():
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4
                    #
                    w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                #
                if self.fluid == "water":
                    data_fluid = self.data_water
                elif self.fluid == "oil":
                    data_fluid = self.data_oil
                elif self.fluid == "gas":
                    data_fluid = self.data_gas
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
            ## Bulk Density, Porosity, Seismic Velocities
            rho_solid = round(rho_s, 3)
            vP, vS, vPvS, rho, var_porosity = SeismicVelocities(
                rho_solid=rho_solid, rho_fluid=self.data_water[2]).calculate_seismic_velocities(
                rho_limits=[2000, 2800], vP_limits=[2500, 6000], vS_limits=[2000, 3000], delta=0.05,
                porosity=porosity)
            ## Elastic Parameters
            bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                rho=rho, vP=vP, vS=vS)
            # Composition data
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            amounts = []
            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)
                chem_data = PeriodicSystem(name=key).get_data()
                amounts.append([key, chem_data[1], value])

            list_oxides = ["H2O", "CO2", "F", "Na2O", "MgO", "Al2O3", "SiO2", "SO3", "K2O", "CaO", "Mn2O3", "Fe2O3",
                           "NiO"]
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
    #
#
## TEST
# print(Sandstone(fluid="water", actualThickness=0).create_sandstone(number=100))
# print(Sandstone(fluid="water", actualThickness=0).create_conglomerate(number=10))