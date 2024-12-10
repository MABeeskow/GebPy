#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		metamorphics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		10.12.2024

# -----------------------------------------------

## MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem, OxideCompounds
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Nesosilicates, Sorosilicates, Phyllosilicates, Inosilicates
from modules.fluids import Water
from modules.geophysics import Elasticity as elast

#######################
## METAMORPHIC ROCKS ##
#######################
class GranuliteFacies:
    #
    def __init__(self, fluid="water", actual_thickness=100, porosity=None):
        self.fluid = fluid
        self.actualThickness = actual_thickness
        self.porosity = porosity
        #
        ## Minerals
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_kyanite = Nesosilicates(impurity="pure", data_type=True).create_kyanite()
        self.data_sillimanite = Nesosilicates(impurity="pure", data_type=True).create_sillimanite()
        self.data_tremolite = Inosilicates(impurity="pure", data_type=True).create_tremolite()
        #
        ## Fluids
        self.data_water = Water.water("")
    #
    def create_granulite(self, number=1, composition=None, classification="felsic"):
        #
        results_container = {}
        results_container["rock"] = "Granulite"
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
            data_garnet_al = Nesosilicates(impurity="pure", data_type=True).create_aluminium_garnet()
            data_clinopyroxene = Inosilicates(
                mineral="Clinopyroxene", data_type=True, traces_list=[]).generate_dataset(
                number=1)
            data_orthopyroxene = Inosilicates(
                mineral="Orthopyroxene", data_type=True, traces_list=[]).generate_dataset(
                number=1)
            #
            mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                          "Grt": data_garnet_al, "Opx": data_orthopyroxene, "Cpx": data_clinopyroxene}
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
                    phi_grt = composition["Grt"]
                    phi_opx = composition["Opx"]
                    phi_cpx = composition["Cpx"]
                    #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Grt"] = phi_grt
                    phi_minerals["Opx"] = phi_opx
                    phi_minerals["Cpx"] = phi_cpx
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        if classification == "felsic":
                            qz_limits = [0.1, 0.26]
                            kfs_limits = [0.2, 0.52]
                            pl_limits = [0.2, 0.52]
                            grt_limits = [0.02, 0.08]
                            opx_limits = [0.0, 0.05]
                            cpx_limits = [0.0, 0.05]
                            #
                            phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                            phi_kfs = round(rd.uniform(kfs_limits[0], (1 - phi_qz)), 4)
                            phi_pl = round(rd.uniform(pl_limits[0], (1 - phi_qz - phi_kfs)), 4)
                            phi_grt = round(rd.uniform(grt_limits[0], (1 - phi_qz - phi_kfs - phi_pl)), 4)
                            phi_opx = round(rd.uniform(opx_limits[0], (1 - phi_qz - phi_kfs - phi_pl - phi_grt)), 4)
                            phi_cpx = round(1 - phi_qz - phi_kfs - phi_pl - phi_grt - phi_opx, 4)
                            #
                        elif classification == "mafic":
                            qz_limits = [0.0, 0.12]
                            kfs_limits = [0.0, 0.04]
                            pl_limits = [0.05, 0.25]
                            grt_limits = [0.02, 0.08]
                            opx_limits = [0.2, 0.52]
                            cpx_limits = [0.2, 0.52]
                            #
                            phi_qz = round(rd.uniform(qz_limits[0], qz_limits[1]), 4)
                            phi_opx = round(rd.uniform(opx_limits[0], (1 - phi_qz)), 4)
                            phi_cpx = round(rd.uniform(cpx_limits[0], (1 - phi_qz - phi_opx)), 4)
                            phi_grt = round(rd.uniform(grt_limits[0], (1 - phi_qz - phi_opx - phi_cpx)), 4)
                            phi_pl = round(rd.uniform(pl_limits[0], (1 - phi_qz - phi_opx - phi_cpx - phi_grt)), 4)
                            phi_kfs = round(1 - phi_qz - phi_opx - phi_cpx - phi_grt - phi_pl, 4)
                        #
                        phi_total = phi_qz + phi_kfs + phi_pl + phi_grt + phi_opx + phi_cpx
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and kfs_limits[0] <= phi_kfs <= kfs_limits[1] \
                                    and pl_limits[0] <= phi_pl <= pl_limits[1] \
                                    and grt_limits[0] <= phi_grt <= grt_limits[1] \
                                    and opx_limits[0] <= phi_opx <= opx_limits[1] \
                                    and cpx_limits[0] <= phi_cpx <= cpx_limits[1]:
                                condition_2 = True
                        #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Grt"] = phi_grt
                    phi_minerals["Opx"] = phi_opx
                    phi_minerals["Cpx"] = phi_cpx
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    try:
                        rho_s += phi_minerals[key]*mineralogy[key]["rho"]
                    except:
                        rho_s += phi_minerals[key]*mineralogy[key]["rho"][0]
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
                                    n_digits = 4
                                else:
                                    n_digits = 4
                                #
                                try:
                                    value = round(w_mineral*mineralogy[mineral]["chemistry"][element], n_digits)
                                except:
                                    value = round(w_mineral*mineralogy[mineral]["chemistry"][element][0], n_digits)
                                #
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
            anisotropic_factor = round(rd.uniform(0.415, 0.425), 2)
            #
            bulk_mod = K_geo/anisotropic_factor
            shear_mod = G_geo/anisotropic_factor
            #
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 6)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 6)

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
class MetamorphicRocks:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        #
        ## Minerals
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_kyanite = Nesosilicates(impurity="pure", data_type=True).create_kyanite()
        self.data_sillimanite = Nesosilicates(impurity="pure", data_type=True).create_sillimanite()
        self.data_tremolite = Inosilicates(impurity="pure", data_type=True).create_tremolite()
        #
        ## Fluids
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
class GreenschistFacies:
    #
    def __init__(self, fluid="water", actual_thickness=100, porosity=None):
        self.fluid = fluid
        self.actualThickness = actual_thickness
        self.porosity = porosity
        #
        ## Minerals
        self.data_brucite = Oxides(impurity="pure", data_type=True).create_brucite()
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_kyanite = Nesosilicates(impurity="pure", data_type=True).create_kyanite()
        self.data_sillimanite = Nesosilicates(impurity="pure", data_type=True).create_sillimanite()
        self.data_diopside = Inosilicates(impurity="pure", data_type=True).create_diopside()
        self.data_tremolite = Inosilicates(impurity="pure", data_type=True).create_tremolite()
        self.data_epidote = Sorosilicates(impurity="pure", data_type=True).create_epidote()
        self.data_chrysotile = Phyllosilicates(impurity="pure", data_type=True).create_chrysotile()
        self.data_pyrophyllite = Phyllosilicates(impurity="pure", data_type=True).create_pyrophyllite()
        self.data_talc = Phyllosilicates(impurity="pure", data_type=True).create_talc()
        #
        ## Fluids
        self.data_water = Water.water("")
    #
    def create_greenschist_basaltic_alt(self, number=1, composition=None):
        results_container = {}
        results_container["rock"] = "Greenschist"
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
            data_actinolite = Inosilicates(impurity="pure", data_type=True).create_actinolite()
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            #
            mineralogy = {"Chl": data_chlorite, "Act": data_actinolite, "Pl": data_plagioclase, "Ep": self.data_epidote}
            #
            minerals_list = list(mineralogy.keys())
            #
            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []
            #
            if self.porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.4), 4)
            else:
                phi_helper = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
            #
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
                    phi_chl = composition["Chl"]
                    phi_act = composition["Act"]
                    phi_pl = composition["Pl"]
                    phi_ep = composition["Ep"]
                    #
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Act"] = phi_act
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Ep"] = phi_ep
                else:
                    condition_2 = False
                    while condition_2 == False:
                        phi_act = round(rd.uniform(0.20, 0.30), 4)
                        phi_chl = round(rd.uniform(0.25, (1.0 - phi_act)), 4)
                        phi_pl = round(rd.uniform(0.2, (1.0 - phi_act - phi_chl)), 4)
                        phi_ep = round(1 - phi_act - phi_chl - phi_pl, 4)
                        phi_total = phi_act + phi_chl + phi_pl + phi_ep
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.2 <= phi_act <= 0.3 and 0.25 <= phi_chl <= 0.5 and 0.2 <= phi_pl <= 0.55 \
                                    and 0.0 <= phi_ep <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Act"] = phi_act
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Ep"] = phi_ep
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
            bulk_mod = 0
            shear_mod = 0
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                bulk_mod += phi_minerals[key]*mineralogy[key]["K"]
                shear_mod += phi_minerals[key]*mineralogy[key]["G"]
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            rho = round((1 - phi_helper)*rho_s + phi_helper*data_fluid[2]/1000, 3)
            bulk_mod = round(bulk_mod, 3)
            shear_mod = round(shear_mod, 3)
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)

            ## RESULTS
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

            results_container["phi"].append(phi_helper)
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
    def create_greenschist_ultramafic_alt(self, number=1, composition=None):
        results_container = {}
        results_container["rock"] = "Greenschist"
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
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            mineralogy = {"Chl": data_chlorite, "Ctl": self.data_chrysotile, "Tlc": self.data_talc,
                          "Tr": self.data_tremolite, "Di": self.data_diopside, "Bru": self.data_brucite}
            #
            minerals_list = list(mineralogy.keys())
            #
            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []
            #
            if self.porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.1), 4)
            else:
                phi_helper = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
            #
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
                    phi_chl = composition["Chl"]
                    phi_ctl = composition["Ctl"]
                    phi_tlc = composition["Tlc"]
                    phi_tr = composition["Tr"]
                    phi_di = composition["Di"]
                    phi_bru = composition["Bru"]
                    #
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Ctl"] = phi_ctl
                    phi_minerals["Tlc"] = phi_tlc
                    phi_minerals["Tr"] = phi_tr
                    phi_minerals["Di"] = phi_di
                    phi_minerals["Bru"] = phi_bru
                    #
                else:
                    condition_2 = False
                    #
                    while condition_2 == False:
                        phi_chl = round(rd.uniform(0.25, 0.50), 4)
                        phi_ctl = round(rd.uniform(0.25, (1.0 - phi_chl)), 4)
                        phi_tlc = round(rd.uniform(0.15, (1.0 - phi_chl - phi_ctl)), 4)
                        phi_tr = round(rd.uniform(0.1, (1.0 - phi_chl - phi_ctl - phi_tlc)), 4)
                        phi_di = round(rd.uniform(0.0, (1.0 - phi_chl - phi_ctl - phi_tlc - phi_tr)), 4)
                        phi_bru = round(1 - phi_chl - phi_ctl - phi_tlc - phi_tr - phi_di, 4)
                        #
                        phi_total = phi_chl + phi_ctl + phi_tlc + phi_tr + phi_di + phi_bru
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.25 <= phi_chl <= 0.5 and 0.25 <= phi_ctl <= 0.5 and 0.15 <= phi_tlc <= 0.2 \
                                    and 0.1 <= phi_tr <= 0.2 and 0.0 <= phi_di <= 0.05 and 0.0 <= phi_bru <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Ctl"] = phi_ctl
                    phi_minerals["Tlc"] = phi_tlc
                    phi_minerals["Tr"] = phi_tr
                    phi_minerals["Di"] = phi_di
                    phi_minerals["Bru"] = phi_bru
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
                n_minerals = len(minerals_list)
                w_minerals_total = 0
                for index, (key, value) in enumerate(phi_minerals.items()):
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4
                    #
                    if index < n_minerals - 1:
                        result = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                        w_minerals[key] = result
                        w_minerals_total += result
                    else:
                        w_minerals[key] = round(1 - w_minerals_total, n_digits)
                #
                total_w_minerals = round(sum(w_minerals.values()), 4)
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
                                #
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        #
                        w_elements[element] = round(w_elements[element], 4)
                #
                total_w_elements = round(sum(w_elements.values()), 4)
                #
                if np.isclose(total_w_minerals, 1.00) == True and np.isclose(total_w_elements, 1.00) == True:
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
            bulk_mod = 0
            shear_mod = 0
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                bulk_mod += phi_minerals[key]*mineralogy[key]["K"]
                shear_mod += phi_minerals[key]*mineralogy[key]["G"]
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            rho = round((1 - phi_helper)*rho_s + phi_helper*data_fluid[2]/1000, 3)
            bulk_mod = round(bulk_mod, 3)
            shear_mod = round(shear_mod, 3)
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)

            ## RESULTS
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

            results_container["phi"].append(phi_helper)
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
    def create_greenschist_pelitic_alt(self, number=1, composition=None):
        results_container = {}
        results_container["rock"] = "Greenschist"
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
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
            data_garnet_al = Nesosilicates(impurity="pure", data_type=True).create_aluminium_garnet()
            #
            mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                          "Chl": data_chlorite, "Ms": data_muscovite, "Grt": data_garnet_al,
                          "Prl": self.data_pyrophyllite}
            #
            minerals_list = list(mineralogy.keys())
            #
            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []
            #
            if self.porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.1), 4)
            else:
                phi_helper = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
            #
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
                    phi_chl = composition["Chl"]
                    phi_ms = composition["Ms"]
                    phi_grt = composition["Grt"]
                    phi_prl = composition["Prl"]
                    #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Ms"] = phi_ms
                    phi_minerals["Grt"] = phi_grt
                    phi_minerals["Prl"] = phi_prl
                    #
                else:
                    condition_2 = False
                    #
                    while condition_2 == False: # Qz (0.25, 0.7), Kfs (0.05, 0.15), Pl (0.05, 0.15), Chl (0.15, 0.3), Ms (0, 0.1), Grt (0, 0.05), Prl (0, 0.05)
                        phi_qz = round(rd.uniform(0.25, 0.7), 4)
                        phi_chl = round(rd.uniform(0.15, (1.0 - phi_qz)), 4)
                        phi_kfs = round(rd.uniform(0.05, (1.0 - phi_qz - phi_chl)), 4)
                        phi_pl = round(rd.uniform(0.05, (1.0 - phi_qz - phi_chl - phi_kfs)), 4)
                        phi_ms = round(rd.uniform(0.0, (1.0 - phi_qz - phi_chl - phi_kfs - phi_pl)), 4)
                        phi_grt = round(rd.uniform(0.0, (1.0 - phi_qz - phi_chl - phi_kfs - phi_pl - phi_ms)), 4)
                        phi_prl = round(1 - phi_qz - phi_chl - phi_kfs - phi_pl - phi_ms - phi_grt, 4)
                        #
                        phi_total = phi_qz + phi_chl + phi_kfs + phi_pl + phi_ms + phi_grt + phi_prl
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.25 <= phi_qz <= 0.7 and 0.15 <= phi_chl <= 0.3 and 0.05 <= phi_kfs <= 0.15 \
                                    and 0.05 <= phi_pl <= 0.15 and 0.0 <= phi_ms <= 0.1 and 0.0 <= phi_grt <= 0.05 \
                                    and 0.0 <= phi_prl <= 0.05:
                                condition_2 = True
                        #
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Chl"] = phi_chl
                    phi_minerals["Ms"] = phi_ms
                    phi_minerals["Grt"] = phi_grt
                    phi_minerals["Prl"] = phi_prl
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
                n_minerals = len(minerals_list)
                w_minerals_total = 0
                for index, (key, value) in enumerate(phi_minerals.items()):
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4
                    #
                    if index < n_minerals - 1:
                        result = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                        w_minerals[key] = result
                        w_minerals_total += result
                    else:
                        w_minerals[key] = round(1 - w_minerals_total, n_digits)
                #
                total_w_minerals = round(sum(w_minerals.values()), 4)
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
                                #
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        #
                        w_elements[element] = round(w_elements[element], 4)
                #
                total_w_elements = round(sum(w_elements.values()), 4)
                #
                if np.isclose(total_w_minerals, 1.00) == True and np.isclose(total_w_elements, 1.00) == True:
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
            bulk_mod = 0
            shear_mod = 0
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                bulk_mod += phi_minerals[key]*mineralogy[key]["K"]
                shear_mod += phi_minerals[key]*mineralogy[key]["G"]
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            rho = round((1 - phi_helper)*rho_s + phi_helper*data_fluid[2]/1000, 3)
            bulk_mod = round(bulk_mod, 3)
            shear_mod = round(shear_mod, 3)
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)

            ## RESULTS
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

            results_container["phi"].append(phi_helper)
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
                        value = round(rd.uniform(0.0, 0.05), 4) # Bru (0.0, 0.05), Chl (0.25, 0.50), Ctl (0.25, 0.50), Tlc (0.15, 0.20), Tr (0.10, 0.20), Di (0.0, 0.05)
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
                if mineral == "Qz": # Qz (0.25, 0.7), Kfs (0.05, 0.15), Pl (0.05, 0.15), Chl (0.15, 0.3), Ms (0, 0.1), Grt (0, 0.05), Prl (0, 0.05)
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
class AmphiboliteFacies:
    #
    def __init__(self, fluid="water", actual_thickness=100, porosity=None):
        self.fluid = fluid
        self.actualThickness = actual_thickness
        self.porosity = porosity
        #
        ## Minerals
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_kyanite = Nesosilicates(impurity="pure", data_type=True).create_kyanite()
        self.data_sillimanite = Nesosilicates(impurity="pure", data_type=True).create_sillimanite()
        self.data_tremolite = Inosilicates(impurity="pure", data_type=True).create_tremolite()
        #
        ## Fluids
        self.data_water = Water.water("")
    #
    def create_amphibolite_ortho(self, number=1, composition=None):
        results_container = {}
        results_container["rock"] = "Ortho-Amphibolite"
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
            data_actinolite = Inosilicates(impurity="pure", data_type=True).create_actinolite()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            mineralogy = {"Tr": self.data_tremolite, "Act": data_actinolite, "Pl": data_plagioclase, "Bt": data_biotite,
                          "Qz": self.data_quartz, "Chl": data_chlorite}
            #
            minerals_list = list(mineralogy.keys())
            #
            if minerals_list[0] not in results_container["mineralogy"]:
                for mineral in minerals_list:
                    results_container["mineralogy"][mineral] = []
            #
            if self.porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.1), 4)
            else:
                phi_helper = round(rd.uniform(self.porosity[0], self.porosity[1]), 4)
            #
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
                    phi_tr = composition["Tr"]
                    phi_act = composition["Act"]
                    phi_pl = composition["Pl"]
                    phi_bt = composition["Bt"]
                    phi_qz = composition["Qz"]
                    phi_chl = composition["Chl"]
                    #
                    phi_minerals["Tr"] = phi_tr
                    phi_minerals["Act"] = phi_act
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Chl"] = phi_chl
                    #
                else:
                    condition_2 = False
                    #
                    while condition_2 == False:
                        phi_tr = round(rd.uniform(0.0, 0.5), 4)
                        phi_act = round(rd.uniform(0.0, (1.0 - phi_tr)), 4)
                        phi_pl = round(rd.uniform(0.15, (1.0 - phi_tr - phi_act)), 4)
                        phi_bt = round(rd.uniform(0.0, (1.0 - phi_tr - phi_act - phi_pl)), 4)
                        phi_qz = round(rd.uniform(0.0, (1.0 - phi_tr - phi_act - phi_pl - phi_bt)), 4)
                        phi_chl = round(1 - phi_tr - phi_act - phi_pl - phi_bt - phi_qz, 4)
                        #
                        phi_total = phi_tr + phi_act + phi_pl + phi_bt + phi_qz + phi_chl
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.0 <= phi_tr <= 0.5 and 0.0 <= phi_act <= 0.5 and 0.15 <= phi_pl <= 0.4 \
                                    and 0.0 <= phi_bt <= 0.2 and 0.0 <= phi_qz <= 0.2 and 0.0 <= phi_chl <= 0.1:
                                condition_2 = True
                        #
                    phi_minerals["Tr"] = phi_tr
                    phi_minerals["Act"] = phi_act
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Bt"] = phi_bt
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Chl"] = phi_chl
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
                n_minerals = len(minerals_list)
                w_minerals_total = 0
                for index, (key, value) in enumerate(phi_minerals.items()):
                    if key == "Urn":
                        n_digits = 4
                    else:
                        n_digits = 4
                    #
                    if index < n_minerals - 1:
                        result = round((phi_minerals[key]*mineralogy[key]["rho"])/rho_s, n_digits)
                        w_minerals[key] = result
                        w_minerals_total += result
                    else:
                        w_minerals[key] = round(1 - w_minerals_total, n_digits)
                #
                total_w_minerals = round(sum(w_minerals.values()), 4)
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
                                #
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        #
                        w_elements[element] = round(w_elements[element], 4)
                #
                total_w_elements = round(sum(w_elements.values()), 4)
                #
                if np.isclose(total_w_minerals, 1.00) == True and np.isclose(total_w_elements, 1.00) == True:
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
            bulk_mod = 0
            shear_mod = 0
            phi_list = []
            for key, value in phi_minerals.items():
                gamma_ray += phi_minerals[key]*mineralogy[key]["GR"]
                photoelectricity += phi_minerals[key]*mineralogy[key]["PE"]
                #
                gamma_ray = round(gamma_ray, 3)
                photoelectricity = round(photoelectricity, 3)
                #
                bulk_mod += phi_minerals[key]*mineralogy[key]["K"]
                shear_mod += phi_minerals[key]*mineralogy[key]["G"]
                K_list.append(round(phi_minerals[key]*mineralogy[key]["K"], 3))
                G_list.append(round(phi_minerals[key]*mineralogy[key]["G"], 3))
                phi_list.append(phi_minerals[key])
            #
            rho = round((1 - phi_helper)*rho_s + phi_helper*data_fluid[2]/1000, 3)
            bulk_mod = round(bulk_mod, 3)
            shear_mod = round(shear_mod, 3)
            youngs_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 3)
            poisson_rat = round((3*bulk_mod - 2*shear_mod)/(6*bulk_mod + 2*shear_mod), 4)
            vP = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho))**0.5, 3)
            vS = round(((shear_mod*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)

            ## RESULTS
            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)

            for key, value in w_elements.items():
                results_container["chemistry"][key].append(value)

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

            results_container["phi"].append(phi_helper)
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
## TEST
# print("Test: Granulite")
# print(MetamorphicRocks().create_granulite(number=1))
# print(MetamorphicRocks().create_granulite(number=10)
# print("Test: Greenschist")
# print(MetamorphicRocks().create_greenschist(number=1))
# print(MetamorphicRocks().create_greenschist(number=10))
# print("Test: Amphibolite")
# print(MetamorphicRocks().create_amphibolite_ortho(number=1))