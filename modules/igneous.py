#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		igneous.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		24.11.2022

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals, geochemistry, oxides
from modules import fluids
from modules.geophysics import Elasticity as elast
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Phyllosilicates, Inosilicates, Nesosilicates
from modules.fluids import Water, Hydrocarbons

class plutonic:
    #
    def __init__(self,):
        pass
    #
    def createGranite(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemQuartz = minerals.oxides.quartz("")
        chemAlkaliFeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        chemPlagioclase = minerals.feldspars.plagioclase(self, "Na")
        chemBiotite = minerals.Biotites.biotite_group(self, "Biotite")
        chemMuscovite = minerals.phyllosilicates.muscovite("")
        chem_actinolite = minerals.inosilicates.actinolite("")
        chem_tremolite = minerals.inosilicates.tremolite("")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = minerals.oxides.water("")
        #
        granite = []
        #
        cond = False
        composition = []
        while cond == False:
            xQuartz = rd.randint(24,27)/100
            xAlkalifeldspar = rd.randint(14,51)/100
            xPlagioclase = rd.randint(10,32)/100
            xBiotite = rd.randint(3,8)/100
            xMuscovite = rd.randint(5,8)/100
            xAmphibole = rd.randint(6,14)/100
            xActinolite2 = rd.randint(0,100)/100
            xTremolite2 = 1-xActinolite2
            xActinolite = xAmphibole*xActinolite2
            xTremolite = xAmphibole*xTremolite2
            sumMin = round(xQuartz + xAlkalifeldspar + xPlagioclase + xBiotite + xMuscovite + xActinolite + xTremolite, 2)
            if sumMin == 1:
                cond = True
                composition.extend([["Qz", round(xQuartz,2), round(chemQuartz[1],2)], ["Kfs", round(xAlkalifeldspar,2), round(chemAlkaliFeldspar[1][0],2), chemAlkaliFeldspar[1][1]], ["Pl", round(xPlagioclase,2), round(chemPlagioclase[1][0],2), chemPlagioclase[1][1]], ["Bt", round(xBiotite,2), round(chemBiotite[1][0],2), chemBiotite[1][1], chemBiotite[1][2]], ["Ms", round(xMuscovite,2), round(chemMuscovite[1],2)], ["Act", round(xActinolite,2), round(chem_actinolite[1][0],2), chem_actinolite[1][1]], ["Tr", round(xTremolite,2), round(chem_tremolite[1],2)]])
            else:
                cond = False
        xQuartz = composition[0][1]
        xAlkalifeldspar = composition[1][1]
        xPlagioclase = composition[2][1]
        xBiotite = composition[3][1]
        xMuscovite = composition[4][1]
        xActinolite = composition[5][1]
        xTremolite = composition[6][1]
        granite.append(composition)
        mineralogy = [chemQuartz, chemAlkaliFeldspar, chemPlagioclase, chemBiotite, chemMuscovite, chem_actinolite, chem_tremolite]
        #
        rhoSolid = (xQuartz*chemQuartz[2] + xAlkalifeldspar*chemAlkaliFeldspar[2] + xPlagioclase*chemPlagioclase[2] + xBiotite*chemBiotite[2] + xMuscovite*chemMuscovite[2] + xActinolite*chem_actinolite[2] + xTremolite*chem_tremolite[2]) / 1000
        X = [xQuartz, xAlkalifeldspar, xPlagioclase, xBiotite, xMuscovite, xActinolite, xTremolite]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = randint(0, 2)/100
        rho = (1 - phi) * rhoSolid + phi * chemWater[2] / 1000
        vP = (1-phi)*vP_solid + phi*chemWater[5]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkaliFeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xBiotite*chemBiotite[5][0] + xMuscovite*chemMuscovite[5][0] + xActinolite*chem_actinolite[5][0] + xTremolite*chem_tremolite[5][0]
        PE = xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkaliFeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xBiotite*chemBiotite[5][1] + xMuscovite*chemMuscovite[5][1] + xActinolite*chem_actinolite[5][1] + xTremolite*chem_tremolite[5][1]
        #poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        #poisson_elastic = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
        poisson_mineralogical = xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkaliFeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xBiotite*chemBiotite[3][3] + xMuscovite*chemMuscovite[3][3] + xActinolite*chem_actinolite[3][3] + xTremolite*chem_tremolite[3][3]
        #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
        #
        granite.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        granite.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        granite.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[4], 2)])
        granite.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        granite.append("water")
        granite.append([GR, PE])
        #granite.append([["Qz", round(chemQuartz[1][0], 2)], ["Kfs", round(chemAlkaliFeldspar[1][0], 2), chemAlkaliFeldspar[1][1]], ["Pl", round(chemPlagioclase[1][0], 2), chemPlagioclase[1][1]], ["Bt", round(chemBiotite[1][0], 2), chemBiotite[1][1], chemBiotite[1][2]], ["Ms", round(chemMuscovite[1][0], 2)], ["Act", round(chem_actinolite[1][0], 2), chem_actinolite[1][1]], ["Tr", round(chem_tremolite[1][0], 2)]])
        #
        #  granite = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return granite
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
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)
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
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)
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
                                         enrichment_pl=None, upper_streckeisen=True):
        results_container = {}
        results_container["rock"] = rock
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
                        elif rock == "Gabbro":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.52, 1.0]
                            bt_limits = [0.0, 0.05]
                        elif rock == "Diorite":
                            qz_limits = [0.0, 0.2]
                            kfs_limits = [0.0, 0.35]
                            pl_limits = [0.52, 1.0]
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
                        elif rock == "Foid-bearing Syenite":
                            nph_limits = [0.0, 0.1]
                            kfs_limits = [0.58, 1.0]
                            pl_limits = [0.0, 0.42]
                            bt_limits = [0.0, 0.05]
                        #
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
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)
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
    def create_granite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                        if 0.2 <= phi_qz <= 0.6 and 0.15 <= phi_kfs <= 0.8 and 0.0 <= phi_pl <= 0.52 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Granite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_granodiorite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                        if 0.2 <= phi_qz <= 0.6 and 0.05 <= phi_kfs <= 0.28 and 0.25 <= phi_pl <= 0.72 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Granodiorite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_tonalite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.2, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.35, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.2 <= phi_qz <= 0.6 and 0.0 <= phi_kfs <= 0.08 and 0.35 <= phi_pl <= 0.8 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Tonalite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_gabbro_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Ca")
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_bt = round(rd.uniform(0.35, 0.5), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_bt)), 4)
                    phi_pl = round(rd.uniform(0.52, (1.0 - phi_bt - phi_kfs)), 4)
                    phi_qz = round(1 - phi_bt - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.0 <= phi_kfs <= 0.35 and 0.52 <= phi_pl <= 1.0 and 0.35 <= phi_bt <= 0.5:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Gabbro"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_diorite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Na")
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_bt = round(rd.uniform(0.0, 0.35), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_bt)), 4)
                    phi_pl = round(rd.uniform(0.52, (1.0 - phi_bt - phi_kfs)), 4)
                    phi_qz = round(1 - phi_bt - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.0 <= phi_kfs <= 0.35 and 0.52 <= phi_pl <= 1.0 and 0.0 <= phi_bt <= 0.35:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Diorite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_monzonite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.28, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.28, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.28 <= phi_kfs <= 0.65 and 0.28 <= phi_pl <= 0.65 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Monzonite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_syenite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.52, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.52 <= phi_kfs <= 1.0 and 0.0 <= phi_pl <= 0.35 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Syenite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_granitoid_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.6, 0.9), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.6 <= phi_qz <= 0.9 and 0.0 <= phi_kfs <= 0.4 and 0.0 <= phi_pl <= 0.4 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Granitoid"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_quarzolite_streckeisen(self, number=1, composition=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_qz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.9, 1.0), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.9 <= phi_qz <= 1.0 and 0.0 <= phi_kfs <= 0.1 and 0.0 <= phi_pl <= 0.1 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity) * rho_s + self.porosity * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Quarzolite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foidbearing_syenite_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.0, 0.1), 4)
                    phi_kfs = round(rd.uniform(0.58, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_nph <= 0.1 and 0.58 <= phi_kfs <= 1.0 and 0.0 <= phi_pl <= 0.42 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid-bearing Syenite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foidbearing_monzonite_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.0, 0.1), 4)
                    phi_kfs = round(rd.uniform(0.32, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.32, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_nph <= 0.1 and 0.32 <= phi_kfs <= 0.65 and 0.32 <= phi_pl <= 0.65 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid-bearing Monzonite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foidbearing_monzodiorite_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Na")
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.0, 0.1), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.58, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_nph <= 0.1 and 0.0 <= phi_kfs <= 0.35 and 0.58 <= phi_pl <= 1.0 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid-bearing Monzodiorite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foidbearing_monzogabbro_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Ca")
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.0, 0.1), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.58, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_nph <= 0.1 and 0.0 <= phi_kfs <= 0.35 and 0.58 <= phi_pl <= 1.0 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid-bearing Monzogabbro"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foid_monzosyenite_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.1, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.2, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.1 <= phi_nph <= 0.6 and 0.2 <= phi_kfs <= 0.9 and 0.0 <= phi_pl <= 0.45 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid Monzosyenite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foid_monzodiorite_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Na")
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.1, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.2, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.1 <= phi_nph <= 0.6 and 0.0 <= phi_kfs <= 0.45 and 0.2 <= phi_pl <= 0.9 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid Monzodiorite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foid_monzogabbro_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Ca")
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.1, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.2, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.1 <= phi_nph <= 0.6 and 0.0 <= phi_kfs <= 0.45 and 0.2 <= phi_pl <= 0.9 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid Monzodiorite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
    def create_foidolite_streckeisen(self, number=1, composition=None):
        data_nepheline = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {
            "Nph": data_nepheline, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_nph = composition["Nph"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_nph = round(rd.uniform(0.6, 1.0), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_nph)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_nph - phi_kfs)), 4)
                    phi_bt = round(1 - phi_nph - phi_kfs - phi_pl, 4)
                    phi_total = phi_nph + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.6 <= phi_nph <= 1.0 and 0.0 <= phi_kfs <= 0.4 and 0.0 <= phi_pl <= 0.4 \
                                and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Nph"] = phi_nph
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
            if self.porosity == None:
                self.porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - self.porosity)*rho_s + self.porosity*self.data_water[2]/1000, 3)
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
        results["rock"] = "Foid Monzodiorite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = self.porosity
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
            results["phi"] = self.porosity
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
    #
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
    #
    def create_simple_granite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        quartz = oxides.Oxides(impurity="pure").create_quartz()
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.15, 0.8)), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.52)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.15, 0.8)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.52)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.15, 0.8)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            elif isinstance(self.amounts, dict):
                w_qz = round(abs(np.random.normal(self.amounts["Qz"], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts["Kfs"], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts["Pl"], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts["Bt"], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts["Ms"], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts["Act"], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts["Tr"], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0][2] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1][2] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.20 <= phi_V[0] <= 0.6 and 0.15 <= phi_V[1] <= 0.8 and 0.0 <= phi_V[2] <= 0.52:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
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
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Granite"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_syenite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.52, 1)), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.25, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.52, 1)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.52, 1)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.52 <= phi_V[1] <= 1.0 and 0.0 <= phi_V[2] <= 0.35:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Syenite"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_monzonite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_pl = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.25, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.28 <= phi_V[1] <= 0.65 and 0.28 <= phi_V[2] <= 0.65:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Monzonite"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_gabbro(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        olivine = minerals.nesosilicates.olivine(self, "olivine")
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                #w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_kfs = 0.0
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_pyr = round(w_acc*rd.uniform(0.5, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*(1-w_pyr-w_amph), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.75, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*(1-w_pyr-w_amph), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.75, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            elif isinstance(self.amounts, dict):
                w_qz = round(abs(np.random.normal(self.amounts["Qz"], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts["Kfs"], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts["Pl"], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts["Bt"], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts["Ms"], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts["Act"], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts["Tr"], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.0 <= phi_V[1] <= 0.35 and 0.52 <= phi_V[2] <= 1.0:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = 0.9*K_geo
        G_solid = 0.9*G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Gabbro"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_diorite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(rd.uniform(0, 0.1), 4)
                w_kfs = round(rd.uniform(0, 0.1), 4)
                w_pl = round(rd.uniform(0.4, 0.5), 4)
                w_mica = round(rd.uniform(0, 0.1), 4)
                w_bt = round(w_mica*rd.uniform(0.95, 1), 4)
                w_ms = round(w_mica - w_bt, 4)
                w_amph = round(rd.uniform(0.2, 0.275), 4)
                w_act = round(w_amph*rd.uniform(0.5, 1), 4)
                w_tr = round(w_amph - w_act, 4)
                w_pyr = round(1 - w_qz - w_kfs - w_pl - w_bt - w_ms - w_act - w_tr, 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            elif isinstance(self.amounts, dict):
                w_qz = round(abs(np.random.normal(self.amounts["Qz"], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts["Kfs"], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts["Pl"], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts["Bt"], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts["Ms"], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts["Act"], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts["Tr"], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.0 <= phi_V[1] <= 0.35 and 0.52 <= phi_V[2] <= 1.0:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = 0.75*K_geo
        G_solid = 0.75*G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Diorite"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_granodiorite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.03, 0.28)), 4)
                w_pl = round(abs(rd.uniform(0.26, 0.72)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.03, 0.28)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.26, 0.72)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.03, 0.28)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.20 <= phi_V[0] <= 0.6 and 0.03 <= phi_V[1] <= 0.28 and 0.26 <= phi_V[2] <= 0.72:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Granodiorite"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_tonalite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_pl = round(abs(rd.uniform(0.35, 0.8)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.35, 0.8)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            w_O = round(1 - w_H - w_F - w_Na - w_Mg - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.20 <= phi_V[0] <= 0.6 and 0.0 <= phi_V[1] <= 0.1 and 0.35 <= phi_V[2] <= 0.8:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
            results["rock"] = "Tonalite"
            #
            element_list = ["H", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]
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
    def create_simple_quartzrich_granitoid(self, w_Na=None, w_K=None, w_Ca=None, amounts=None):
        #
        self.w_Na = w_Na
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        muscovite = minerals.phyllosilicates.muscovite("")
        schorl = minerals.Tourmalines.schorl("")
        topaz = minerals.nesosilicates.topaz("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, muscovite, schorl, topaz]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Na == None and self.w_K == None and self.w_Ca == None and self.amounts == None:
                w_qz = round(rd.uniform(0.6, 0.9), 4)
                w_fsp = round(rd.uniform(0, (1-w_qz)), 4)
                w_kfs = round(w_fsp*rd.uniform(0.0, 0.4), 4)
                w_pl = round(w_fsp-w_kfs, 4)
                w_acc = round((1-w_qz-w_fsp), 4)
                w_ms = round(w_acc*rd.uniform(0.25, 1), 4)
                w_trm = round(w_acc*rd.uniform(0, (1-w_ms)), 4)
                w_tpz = round(1-w_qz-w_fsp-w_ms-w_trm, 4)
            elif self.w_Na != None:
                w_qz = round(abs(rd.uniform(0.6, 0.9)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.4)), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.4)), 4)
                w_fsp = w_kfs + w_pl
                w_trm = round((self.w_Na - w_kfs*alkalifeldspar[6][1] - w_pl*plagioclase[6][1])/(schorl[6][3]), 4)
                w_ms = round(abs(rd.uniform(0.0, (1-w_qz-w_fsp-w_trm))), 4)
                w_tpz = round(abs(1-w_qz-w_fsp-w_trm), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0, 0.05)), 4)
                w_ms = round(abs(w_acc*rd.uniform(0.5, 1)), 4)
                w_trm = round(abs(w_acc*rd.uniform(0, (1-w_ms))), 4)
                w_tpz = round(abs(w_acc*(1-w_ms-w_trm)), 4)
                w_kfs = round((self.w_K - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.4)), 4)
                w_fsp = w_kfs + w_pl
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_qz = round(abs(rd.uniform(0.6, 0.9)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.4)), 4)
                w_pl = round((self.w_Ca)/(plagioclase[6][4]), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_ms = round(abs(w_acc*rd.uniform(0.5, 1)), 4)
                w_trm = round(abs(w_acc*rd.uniform(0, (1-w_ms))), 4)
                w_tpz = round(abs(w_acc*(1-w_ms-w_trm)), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_trm = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_tpz = round(1-w_qz-w_kfs-w_pl-w_ms-w_trm, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_ms >= 0.0 and w_trm >= 0.0 and w_tpz >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_ms + w_trm + w_tpz, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ms*muscovite[6][0] + w_trm*schorl[6][0] + w_tpz*topaz[6][0], 4)
            w_B = round(w_trm*schorl[6][1], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_ms*muscovite[6][1] + w_trm*schorl[6][2] + w_tpz*topaz[6][1], 4)
            w_F = round(w_ms*muscovite[6][2] + w_tpz*topaz[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1] + w_trm*schorl[6][3], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_ms*muscovite[6][3] + w_trm*schorl[6][4] + w_tpz*topaz[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_ms*muscovite[6][4] + w_trm*schorl[6][5] + w_tpz*topaz[6][4], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4], 4)
            w_Fe = round(w_trm*schorl[6][6], 4)
            w_O = round(1 - w_H - w_B - w_Fe - w_Na - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_B + w_O + w_F + w_Na + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Trm", w_trm, round(schorl[1], 2)], ["Tpz", w_tpz, round(topaz[1][0], 2), round(topaz[1][1], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Ms", "Trm", "Tpz"]))
                concentrations = [w_H, w_B, w_O, w_F, w_Na, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_ms, w_trm, w_tpz]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.6 <= phi_V[0] <= 0.9 and 0.0 <= phi_V[1] <= 0.4 and 0.0 <= phi_V[2] <= 0.4:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_ms*muscovite[2] + w_trm*schorl[2] + w_tpz*topaz[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_ms, w_trm, w_tpz]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_ms*muscovite[5][0] + w_trm*schorl[5][0] + w_tpz*topaz[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_ms*muscovite[5][1] + w_trm*schorl[5][1] + w_tpz*topaz[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_ms*muscovite[3][3] + w_trm*schorl[3][3] + w_tpz*topaz[3][3]
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
            results["rock"] = "Qz-rich Granitoid"
            #
            element_list = ["H", "B", "O", "F", "Na", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Ms", "Trm", "Tpz"]
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
    def create_simple_quartzolite(self, w_Na=None, w_K=None, w_Ca=None, amounts=None):
        #
        self.w_Na = w_Na
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        muscovite = minerals.phyllosilicates.muscovite("")
        schorl = minerals.Tourmalines.schorl("")
        topaz = minerals.nesosilicates.topaz("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, muscovite, schorl, topaz]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Na == None and self.w_K == None and self.w_Ca == None and self.amounts == None:
                w_qz = round(rd.uniform(0.9, 1), 4)
                w_fsp = round(rd.uniform(0, (1-w_qz)), 4)
                w_kfs = round(w_fsp*rd.uniform(0.0, 0.1), 4)
                w_pl = round(w_fsp-w_kfs, 4)
                w_acc = round((1-w_qz-w_fsp), 4)
                w_ms = round(w_acc*rd.uniform(0.25, 1), 4)
                w_trm = round(w_acc*rd.uniform(0, (1-w_ms)), 4)
                w_tpz = round(1-w_qz-w_fsp-w_ms-w_trm, 4)
            elif self.w_Na != None:
                w_qz = round(abs(rd.uniform(0.9, 1)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_fsp = w_kfs + w_pl
                w_trm = round((self.w_Na - w_kfs*alkalifeldspar[6][1] - w_pl*plagioclase[6][1])/(schorl[6][3]), 4)
                w_ms = round(abs(rd.uniform(0.0, (1-w_qz-w_fsp-w_trm))), 4)
                w_tpz = round(abs(1-w_qz-w_fsp-w_trm), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0, 0.05)), 4)
                w_ms = round(abs(w_acc*rd.uniform(0.5, 1)), 4)
                w_trm = round(abs(w_acc*rd.uniform(0, (1-w_ms))), 4)
                w_tpz = round(abs(w_acc*(1-w_ms-w_trm)), 4)
                w_kfs = round((self.w_K - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_fsp = w_kfs + w_pl
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_qz = round(abs(rd.uniform(0.9, 1)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_pl = round((self.w_Ca)/(plagioclase[6][4]), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_ms = round(abs(w_acc*rd.uniform(0.5, 1)), 4)
                w_trm = round(abs(w_acc*rd.uniform(0, (1-w_ms))), 4)
                w_tpz = round(abs(w_acc*(1-w_ms-w_trm)), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_trm = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_tpz = round(1-w_qz-w_kfs-w_pl-w_ms-w_trm, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_ms >= 0.0 and w_trm >= 0.0 and w_tpz >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_ms + w_trm + w_tpz, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ms*muscovite[6][0] + w_trm*schorl[6][0] + w_tpz*topaz[6][0], 4)
            w_B = round(w_trm*schorl[6][1], 4)
            #w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_ms*muscovite[6][1] + w_trm*schorl[6][2] + w_tpz*topaz[6][1], 4)
            w_F = round(w_ms*muscovite[6][2] + w_tpz*topaz[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1] + w_trm*schorl[6][3], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_ms*muscovite[6][3] + w_trm*schorl[6][4] + w_tpz*topaz[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_ms*muscovite[6][4] + w_trm*schorl[6][5] + w_tpz*topaz[6][4], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4], 4)
            w_Fe = round(w_trm*schorl[6][6], 4)
            w_O = round(1 - w_H - w_B - w_F - w_Na - w_Al - w_Si - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_B + w_O + w_F + w_Na + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Trm", w_trm, round(schorl[1], 2)], ["Tpz", w_tpz, round(topaz[1][0], 2), round(topaz[1][1], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Ms", "Trm", "Tpz"]))
                concentrations = [w_H, w_B, w_O, w_F, w_Na, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_ms, w_trm, w_tpz]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.9 <= phi_V[0] <= 1.0 and 0.0 <= phi_V[1] <= 0.1 and 0.0 <= phi_V[2] <= 0.1:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_ms*muscovite[2] + w_trm*schorl[2] + w_tpz*topaz[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_ms, w_trm, w_tpz]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
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
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_ms*muscovite[5][0] + w_trm*schorl[5][0] + w_tpz*topaz[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_ms*muscovite[5][1] + w_trm*schorl[5][1] + w_tpz*topaz[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_ms*muscovite[3][3] + w_trm*schorl[3][3] + w_tpz*topaz[3][3]
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
            results["rock"] = "Quartzolite"
            #
            element_list = ["H", "B", "O", "F", "Na", "Al", "Si", "K", "Ca", "Fe"]
            mineral_list = ["Qz", "Kfs", "Pl", "Ms", "Trm", "Tpz"]
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
class Volcanic:
    #
    def __init__(self, fluid, actualThickness):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_water = Water.water("")
    #
    def create_volcanic_rock(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                phi_pl= composition["Pl"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_phi_mineral = False
                while condition_phi_mineral == False:
                    phi_qz = round(rd.uniform(0.0, 1.0), 4)
                    phi_kfs = round(rd.uniform(0.0, (1 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs , 4)
                    #
                    phi_total = phi_qz + phi_kfs + phi_pl
                    if np.isclose(phi_total, 1.0000) == True:
                        condition_phi_mineral = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            #
            rho_s = 0
            for key, value in phi_minerals.items():
                rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                for element, value in mineralogy[key]["chemistry"].items():
                    if element not in elements_list:
                        elements_list.append(element)
                        w_elements[element] = 0.0
            rho_s = round(rho_s, 3)
            for key, value in phi_minerals.items():
                w_minerals[key] = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, 4)
            #
            if porosity == None:
                porosity = round(rd.uniform(0.0, 0.2), 4)
            rho = round((1 - porosity)*rho_s + porosity*self.data_water[2]/1000, 3)
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
        youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
        poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
        vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
        vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
        vPvS = round(vP/vS, 3)
        #
        results = {}
        results["rock"] = "Volcanic"
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
    #
    def create_rhyolite(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.2, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.35, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.2 <= phi_qz <= 0.6 and 0.35 <= phi_kfs <= 0.90 and 0.1 <= phi_pl <= 0.65:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Rhyolite"
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
    #
    def create_alkaline_rhyolite(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.2, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.4, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.2 <= phi_qz <= 0.6 and 0.4 <= phi_kfs <= 0.8 and 0.0 <= phi_pl <= 0.1:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Alkaline-Rhyolite"
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
    #
    def create_dacite(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.2, 0.6), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.25, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.2 <= phi_qz <= 0.6 and 0.0 <= phi_kfs <= 0.28 and 0.25 <= phi_pl <= 0.8 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Dacite"
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
    #
    def create_basalt(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.52, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.0 <= phi_kfs <= 0.35 and 0.52 <= phi_pl <= 1.0 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Basalt"
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
    #
    def create_andesite(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                phi_bt = composition["Bt"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.52, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.0 <= phi_kfs <= 0.35 and 0.52 <= phi_pl <= 1.0 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Andesite"
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
    #
    def create_trachyte(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.05), 4)
                    phi_kfs = round(rd.uniform(0.62, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.05 and 0.62 <= phi_kfs <= 0.9 and 0.1 <= phi_pl <= 0.35:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Trachyte"
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
    #
    def create_alkaline_trachyte(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.05), 4)
                    phi_kfs = round(rd.uniform(0.85, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.05 and 0.85 <= phi_kfs <= 1.0 and 0.0 <= phi_pl <= 0.1:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Alkaline-Trachyte"
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
    #
    def create_quartzitic_trachyte(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.05, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.52, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.05 <= phi_qz <= 0.2 and 0.52 <= phi_kfs <= 0.85 and 0.1 <= phi_pl <= 0.35:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Quartzitic-Trachyte"
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
    #
    def create_quartzitic_alkaline_trachyte(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.05, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.7, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.05 <= phi_qz <= 0.2 and 0.7 <= phi_kfs <= 0.95 and 0.0 <= phi_pl <= 0.1:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Quartzitic-Alkaline-Trachyte"
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
    #
    def create_latite(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.05), 4)
                    phi_kfs = round(rd.uniform(0.32, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.05 and 0.32 <= phi_kfs <= 0.65 and 0.32 <= phi_pl <= 0.65:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Latite"
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
    #
    def create_quartzitic_latite(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase}
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
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.05, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.32, (1.0 - phi_qz)), 4)
                    phi_pl = round(1 - phi_qz - phi_kfs, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.05 <= phi_qz <= 0.2 and 0.28 <= phi_kfs <= 0.58 and 0.28 <= phi_pl <= 0.58:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Latite"
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
    #
    def create_rhyolite_generalized(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
                        if 0.2 <= phi_qz <= 0.6 and 0.15 <= phi_kfs <= 0.8 and 0.0 <= phi_pl <= 0.52 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Rhyolite"
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
    #
    def create_trachyte_generalized(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.52, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.52 <= phi_kfs <= 1.0 and 0.0 <= phi_pl <= 0.35 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Trachyte"
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
    #
    def create_latite_generalized(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite}
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
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.0, 0.2), 4)
                    phi_kfs = round(rd.uniform(0.28, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.28, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(1 - phi_qz - phi_kfs - phi_pl, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.0 <= phi_qz <= 0.2 and 0.28 <= phi_kfs <= 0.65 and 0.28 <= phi_pl <= 0.65 and 0.0 <= phi_bt <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
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
                porosity = round(rd.uniform(0.0, 0.2), 4)
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
        results["rock"] = "Latite"
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
    #
    def create_simple_basalt(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_pyr = round(w_acc*rd.uniform(0.5, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*(1-w_pyr-w_amph), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*(1-w_pyr-w_amph), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_pyr = round(w_acc*rd.uniform(0, 1), 4)
                w_aug = w_pyr
                w_amph = round(w_acc*rd.uniform(0, (1-w_pyr)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Bt", "Ms", "Act", "Tr", "Aug"]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.0 <= phi_V[1] <= 0.35 and 0.52 <= phi_V[2] <= 1.0:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = 0.775*np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = 0.775*np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if self.actualThickness <= 1000:
            phi = rd.uniform(0.08, 0.10)
        elif self.actualThickness > 1000 and self.actualThickness <= 2000:
            phi = rd.uniform(0.06, 0.08)
        elif self.actualThickness > 2000 and self.actualThickness <= 3000:
            phi = rd.uniform(0.04, 0.06)
        elif self.actualThickness > 3000 and self.actualThickness <= 4000:
            phi = rd.uniform(0.02, 0.04)
        elif self.actualThickness > 4000:
            phi = rd.uniform(0.0, 0.02)
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
        return data
    #
    def create_simple_andesite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.52, 1.0)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0.5, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.0 <= phi_V[1] <= 0.35 and 0.52 <= phi_V[2] <= 1.0:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if self.actualThickness <= 1000:
            phi = rd.uniform(0.125, 0.15)
        elif self.actualThickness > 1000 and self.actualThickness <= 2000:
            phi = rd.uniform(0.10, 0.125)
        elif self.actualThickness > 2000 and self.actualThickness <= 3000:
            phi = rd.uniform(0.075, 0.10)
        elif self.actualThickness > 3000 and self.actualThickness <= 4000:
            phi = rd.uniform(0.05, 0.075)
        elif self.actualThickness > 4000:
            phi = rd.uniform(0.025, 0.05)
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
        return data
    #
    def create_simple_latite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_pl = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.28, 0.65)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.28 <= phi_V[1] <= 0.65 and 0.28 <= phi_V[2] <= 0.65:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
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
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
        return data
    #
    def create_simple_trachyte(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.52, 1)), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.52, 1)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.35)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.0, 0.2)), 4)
                w_kfs = round(abs(rd.uniform(0.52, 1)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.0 <= phi_V[0] <= 0.2 and 0.52 <= phi_V[1] <= 1.0 and 0.0 <= phi_V[2] <= 0.35:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
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
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
        return data
    #
    def create_simple_rhyolite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.15, 0.8)), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.52)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.15, 0.8)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.0, 0.52)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.15, 0.8)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.20 <= phi_V[0] <= 0.6 and 0.15 <= phi_V[1] <= 0.8 and 0.0 <= phi_V[2] <= 0.52:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
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
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
        return data
    #
    def create_simple_dacite(self, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        plagioclase = minerals.feldspars.plagioclase(self, "Na")
        biotite = minerals.Biotites.biotite_group(self, "Biotite")
        muscovite = minerals.phyllosilicates.muscovite("")
        actinolite = minerals.inosilicates.actinolite("")
        tremolite = minerals.inosilicates.tremolite("")
        augite = minerals.inosilicates.augite("")
        #
        mineralogy = [quartz, alkalifeldspar, plagioclase, biotite, muscovite, actinolite, tremolite, augite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.28)), 4)
                w_pl = round(abs(rd.uniform(0.25, 0.8)), 4)
                w_fsp = w_kfs + w_pl
                w_acc = round((1-w_qz-w_fsp), 4)
                w_mica = round(w_acc*rd.uniform(0.5, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
            elif self.w_Mg != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Mg - w_act*actinolite[6][2] - w_tr*tremolite[6][2] - w_aug*augite[6][1])/(biotite[6][3]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.28)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif self.w_K != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_mica = round(w_acc*rd.uniform(0, 1), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0.5, 1)), 4)
                w_ms = round(abs(w_acc*w_mica*rd.uniform(0, (1-w_bt))), 4)
                w_amph = round(w_acc*rd.uniform(0, (1-w_mica)), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(abs(w_acc*(1-w_mica-w_amph)), 4)
                w_aug = w_pyr
                w_kfs = round((self.w_K - w_bt*biotite[6][6] - w_ms*muscovite[6][5])/(alkalifeldspar[6][4]), 4)
                w_pl = round(abs(rd.uniform(0.25, 0.8)), 4)
                w_fsp = w_kfs + w_pl
                #
                w_qz = round(abs(1-w_acc-w_fsp), 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_pl = round((self.w_Ca - w_act*actinolite[6][4] - w_tr*tremolite[6][4] - w_aug*augite[6][3])/(plagioclase[6][4]), 4)
                w_mica = round(w_acc*(1-w_amph-w_pyr), 4)
                w_bt = round(abs(w_acc*w_mica*rd.uniform(0, 1)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(1-w_acc-w_qz-w_pl), 4)
            elif self.w_Fe != None:
                w_acc = round(abs(rd.uniform(0.0, 0.05)), 4)
                w_amph = round(w_acc*rd.uniform(0, 1), 4)
                w_act = round(abs(w_acc*w_amph*rd.uniform(0, 1)), 4)
                w_tr = round(abs(w_acc*w_amph*rd.uniform(0, (1-w_act))), 4)
                w_pyr = round(w_acc*rd.uniform(0, (1-w_amph)), 4)
                w_aug = w_pyr
                w_bt = round((self.w_Fe - w_act*actinolite[6][5] - w_aug*augite[6][4])/(biotite[6][7]), 4)
                w_mica = round(w_acc*rd.uniform(0, (1-w_amph-w_pyr)), 4)
                w_ms = round(abs(w_acc*(w_mica-w_bt)), 4)
                #
                w_qz = round(abs(rd.uniform(0.2, 0.6)), 4)
                w_kfs = round(abs(rd.uniform(0.0, 0.28)), 4)
                w_pl = round(abs(1-w_acc-w_qz-w_kfs), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ms = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_act = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_tr = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_aug = round(1-w_qz-w_kfs-w_pl-w_bt-w_ms-w_act-w_tr, 4)
            #
            if w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_act >= 0.0 and w_tr >= 0.0 and w_aug >= 0.0:
                sumMin = round(w_qz + w_kfs + w_pl + w_bt + w_ms + w_act + w_tr + w_aug, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_bt*biotite[6][0] + w_ms*muscovite[6][0] + w_act*actinolite[6][0] + w_tr*tremolite[6][0], 4)
            w_O = round(w_qz*quartz[6][0] + w_kfs*alkalifeldspar[6][0] + w_pl*plagioclase[6][0] + w_bt*biotite[6][1] + w_ms*muscovite[6][1] + w_act*actinolite[6][1] + w_tr*tremolite[6][1] + w_aug*augite[6][0], 4)
            w_F = round(w_bt*biotite[6][2] + w_ms*muscovite[6][2], 4)
            w_Na = round(w_kfs*alkalifeldspar[6][1] + w_pl*plagioclase[6][1], 4)
            w_Mg = round(w_bt*biotite[6][3] + w_act*actinolite[6][2] + w_tr*tremolite[6][2] + w_aug*augite[6][1], 4)
            w_Al = round(w_kfs*alkalifeldspar[6][2] + w_pl*plagioclase[6][2] + w_bt*biotite[6][4] + w_ms*muscovite[6][3], 4)
            w_Si = round(w_qz*quartz[6][1] + w_kfs*alkalifeldspar[6][3] + w_pl*plagioclase[6][3] + w_bt*biotite[6][5] + w_ms*muscovite[6][4] + w_act*actinolite[6][3] + w_tr*tremolite[6][3] + w_aug*augite[6][2], 4)
            w_K = round(w_kfs*alkalifeldspar[6][4] + w_bt*biotite[6][6] + w_ms*muscovite[6][5], 4)
            w_Ca = round(w_pl*plagioclase[6][4] + w_act*actinolite[6][4] + w_tr*tremolite[6][4] + w_aug*augite[6][3], 4)
            w_Fe = round(w_bt*biotite[6][7] + w_act*actinolite[6][5] + w_aug*augite[6][4], 4)
            sumConc = round(w_H + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Kfs", w_kfs, round(alkalifeldspar[1][0], 2), round(alkalifeldspar[1][1], 2)], ["Pl", w_pl, round(plagioclase[1][0], 2), round(plagioclase[1][1], 2)], ["Bt", w_bt, round(biotite[1][0], 2), round(biotite[1][1], 2), round(biotite[1][2], 2)], ["Ms", w_ms, round(muscovite[1], 2)], ["Act", w_act, round(actinolite[1][0], 2), round(actinolite[1][1], 2)], ["Tr", w_tr, round(tremolite[1], 2)], ["Aug", w_aug, round(augite[1][0], 2), round(augite[1][1], 2), round(augite[1][2], 2), round(augite[1][3], 2)]))
                concentrations = [w_H, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
                amounts = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr, w_aug]
                phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                #print(np.around(phi_V, 4))
                if 0.20 <= phi_V[0] <= 0.6 and 0.0 <= phi_V[1] <= 0.28 and 0.25 <= phi_V[2] <= 0.8:
                    cond = True
                else:
                    composition = []
                    cond = False
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_qz*quartz[2] + w_kfs*alkalifeldspar[2] + w_pl*plagioclase[2] + w_bt*biotite[2] + w_ms*muscovite[2] + w_act*actinolite[2] + w_tr*tremolite[2] + w_aug*augite[2]) / 1000
        X = [w_qz, w_kfs, w_pl, w_bt, w_ms, w_act, w_tr]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
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
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_kfs*alkalifeldspar[5][0] + w_pl*plagioclase[5][0] + w_bt*biotite[5][0] + w_ms*muscovite[5][0] + w_act*actinolite[5][0] + w_tr*tremolite[5][0] + w_aug*augite[5][0]
        PE = w_qz*quartz[5][1] + w_kfs*alkalifeldspar[5][1] + w_pl*plagioclase[5][1] + w_bt*biotite[5][1] + w_ms*muscovite[5][1] + w_act*actinolite[5][1] + w_tr*tremolite[5][1] + w_aug*augite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_kfs*alkalifeldspar[3][3] + w_pl*plagioclase[3][3] + w_bt*biotite[3][3] + w_ms*muscovite[3][3] + w_act*actinolite[3][3] + w_tr*tremolite[3][3] + w_aug*augite[3][3]
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
        return data
#
class volcanic:
    #
    def __init__(self,):
        pass
    #
    def createBasalt(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chem_plagioclase = minerals.feldspars.plagioclase(self, "Ca")
        chem_enstatite = minerals.inosilicates.enstatite("")
        chem_ferrosilite = minerals.inosilicates.ferrosilite("")
        chem_biotite = minerals.Biotites.biotite_group(self, "Biotite")
        chem_actinolite = minerals.inosilicates.actinolite("")
        chem_tremolite = minerals.inosilicates.tremolite("")
        chem_olivine = minerals.nesosilicates.olivine(self, "Olivine")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = minerals.oxides.water("")
        #
        basalt = []
        #
        cond = False
        composition = []
        while cond == False:
            xPlagioclase = rd.randint(13,54)/100
            xPyroxene = rd.randint(9,74)/100
            xEnstatite2 = rd.randint(0,100)/100
            xFerrosilite2 = 1-xEnstatite2
            xEnstatite = xPyroxene*xEnstatite2
            xFerrosilite = xPyroxene*xFerrosilite2
            xBiotite = rd.randint(0,4)/100
            xAmphibole = rd.randint(0,23)/100
            xActinolite2 = rd.randint(0,100)/100
            xTremolite2 = 1-xActinolite2
            xActinolite = xAmphibole*xActinolite2
            xTremolite = xAmphibole*xTremolite2
            xOlivine = rd.randint(0,14)/100
            sumMin = round(xPlagioclase + (xEnstatite+xFerrosilite) + xBiotite + (xActinolite+xTremolite) + xOlivine, 2)
            if sumMin == 1:
                cond = True
                composition.extend([["Pl", round(xPlagioclase,2), round(chem_plagioclase[1][0],2), chem_plagioclase[1][1]], ["En", round(xEnstatite,2), round(chem_enstatite[1],2)], ["Fs", round(xFerrosilite,2), round(chem_ferrosilite[1],2)], ["Bt", round(xBiotite,2), round(chem_biotite[1][0],2), round(chem_biotite[1][1],2), round(chem_biotite[1][2],2)], ["Act", round(xActinolite,2), round(chem_actinolite[1][0],2), chem_actinolite[1][1]], ["Tr", round(xTremolite,2), round(chem_tremolite[1],2)], ["Ol", round(xOlivine,2), round(chem_olivine[1][0],2), round(chem_olivine[1][1],2), round(chem_olivine[1][2],2), round(chem_olivine[1][3],2)]])
            else:
                cond = False
        xPlagioclase = composition[0][1]
        xEnstatite = composition[1][1]
        xFerrosilite = composition[2][1]
        xBiotite = composition[3][1]
        xActinolite = composition[4][1]
        xTremolite = composition[5][1]
        xOlivine = composition[6][1]
        mineralogy = [chem_plagioclase, chem_enstatite, chem_ferrosilite, chem_biotite, chem_actinolite, chem_tremolite, chem_olivine]
        #
        basalt.append(composition)
        #
        rhoSolid = 1.15*(xPlagioclase*chem_plagioclase[2] + xEnstatite*chem_enstatite[2] + xFerrosilite*chem_ferrosilite[2] + xBiotite*chem_biotite[2] + xActinolite*chem_actinolite[2] + xTremolite*chem_tremolite[2] + xOlivine*chem_olivine[2]) / 1000
        X = [xPlagioclase, xEnstatite, xFerrosilite, xBiotite, xActinolite, xTremolite, xOlivine]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = randint(0, 15)/100
        rho = (1 - phi) * rhoSolid + phi * chemWater[2] / 1000
        vP = (1-phi)*vP_solid + phi*chemWater[5]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = xPlagioclase*chem_plagioclase[5][0] + xEnstatite*chem_enstatite[5][0] + xFerrosilite*chem_ferrosilite[5][0] + xBiotite*chem_biotite[5][0] + xActinolite*chem_actinolite[5][0] + xTremolite*chem_tremolite[5][0] + xOlivine*chem_olivine[5][0]
        PE = xPlagioclase*chem_plagioclase[5][1] + xEnstatite*chem_enstatite[5][1] + xFerrosilite*chem_ferrosilite[5][1] + xBiotite*chem_biotite[5][1] + xActinolite*chem_actinolite[5][1] + xTremolite*chem_tremolite[5][1] + xOlivine*chem_olivine[5][1]
        #poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        #poisson_elastic = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
        poisson_mineralogical = xPlagioclase*chem_plagioclase[3][3] + xEnstatite*chem_enstatite[3][3] + xFerrosilite*chem_ferrosilite[3][3] + xBiotite*chem_biotite[3][3] + xActinolite*chem_actinolite[3][3] + xTremolite*chem_tremolite[3][3] + xOlivine*chem_olivine[3][3]
        #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
        #
        basalt.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        basalt.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        basalt.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[4], 2)])
        basalt.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        basalt.append("water")
        basalt.append([GR, PE])
        #
        #  basalt = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return basalt
#
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
    #