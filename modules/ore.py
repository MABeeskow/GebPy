#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		ore.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		10.12.2024

# -----------------------------------------------

# MODULES
import time
import numpy as np
import random as rd
from modules.silicates import Phyllosilicates, Tectosilicates, Nesosilicates, Sorosilicates, Inosilicates, Cyclosilicates
from modules.oxides import Oxides
from modules.carbonates import Carbonates
from modules.sulfides import Sulfides
from modules.phosphates import Phosphates
from modules import fluids
from modules.fluids import Water
from modules.geophysics import Elasticity as elast
from modules.chemistry import PeriodicSystem, OxideCompounds
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry
from modules.pyllosilicates import Pyllosilicates
from modules.petrophysics import SeismicVelocities

class OreRocks:
    #
    def __init__(self, fluid, actual_thickness, porosity):
        ## Settings
        self.fluid = fluid
        self.actual_thickness = actual_thickness
        self.porosity = porosity
        #
        ## Minerals
        self.data_quartz = Oxides(data_type=True).create_quartz()
        self.data_magnetite = Oxides(data_type=True).create_magnetite()
        self.data_hematite = Oxides(data_type=True).create_hematite()
        self.data_goethite = Oxides(data_type=True).create_goethite()
        self.data_gibbsite = Oxides(data_type=True).create_gibbsite()
        self.data_boehmite = Oxides(data_type=True).create_boehmite()
        self.data_anatase = Oxides(data_type=True).create_anatase()
        self.data_calcite = Carbonates(data_type=True).create_calcite()
        self.data_siderite = Carbonates(data_type=True).create_siderite()
        self.data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        self.data_pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()
        #
        ## Fluids
        #
        self.data_water = Water.water("")
    #
    def create_siliciclastic_itabirite(self, rock="Itabirite", number=1, composition=None, classification="Itabirite"):
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

        n = 0

        classification = rock

        while n < number:
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            mineralogy = {"Hem": self.data_hematite, "Goe": self.data_goethite, "Kln": self.data_kaolinite,
                          "Qz": self.data_quartz, "Mag": self.data_magnetite, "Gbs": self.data_gibbsite,
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
                    phi_hem = composition["Hem"]
                    phi_goe = composition["Goe"]
                    phi_kln = composition["Kln"]
                    phi_qz = composition["Qz"]
                    phi_mag = composition["Mag"]
                    phi_gbs = composition["Gbs"]
                    phi_bt = composition["Bt"]
                    #
                    phi_minerals["Hem"] = phi_hem
                    phi_minerals["Goe"] = phi_goe
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Mag"] = phi_mag
                    phi_minerals["Gbs"] = phi_gbs
                    phi_minerals["Bt"] = phi_bt
                    #
                else:
                    condition_2 = False
                    while condition_2 == False:
                        if classification == "Itabirite":
                            gbs_limits = [0.0, 0.05]
                            qz_limits = [0.0, 0.4]
                            mag_limits = [0.0, 0.05]
                            hem_limits = [0.3, 0.7]
                            goe_limits = [0.05, 0.26]
                            kln_limits = [0.0, 0.14]
                            bt_limits = [0.0, 0.05]
                        elif classification == "Banded Iron Formation":
                            gbs_limits = [0.0, 0.05]
                            qz_limits = [0.4, 0.75]
                            mag_limits = [0.15, 0.4]
                            hem_limits = [0.05, 0.20]
                            goe_limits = [0.0, 0.05]
                            kln_limits = [0.0, 0.15]
                            bt_limits = [0.0, 0.05]
                        elif classification == "Compact Hematite":
                            gbs_limits = [0.0, 0.05]
                            qz_limits = [0.0, 0.05]
                            mag_limits = [0.0, 0.05]
                            hem_limits = [0.7, 0.8]
                            goe_limits = [0.05, 0.15]
                            kln_limits = [0.0, 0.1]
                            bt_limits = [0.0, 0.05]
                        elif classification == "Friable Hematite":
                            gbs_limits = [0.0, 0.02]
                            qz_limits = [0.0, 0.02]
                            mag_limits = [0.0, 0.02]
                            hem_limits = [0.85, 0.95]
                            goe_limits = [0.04, 0.06]
                            kln_limits = [0.0, 0.02]
                            bt_limits = [0.0, 0.02]
                        elif classification == "Goethite Hematite":
                            gbs_limits = [0.05, 0.09]
                            qz_limits = [0.0, 0.02]
                            mag_limits = [0.0, 0.02]
                            hem_limits = [0.66, 0.72]
                            goe_limits = [0.18, 0.22]
                            kln_limits = [0.0, 0.02]
                            bt_limits = [0.0, 0.02]
                        elif classification == "Al-rich Itabirite":
                            gbs_limits = [0.0, 0.02]
                            qz_limits = [0.0, 0.08]
                            mag_limits = [0.0, 0.05]
                            hem_limits = [0.6, 0.67]
                            goe_limits = [0.1, 0.18]
                            kln_limits = [0.1, 0.18]
                            bt_limits = [0.0, 0.02]
                        elif classification == "Compact Quartz Itabirite":
                            gbs_limits = [0.0, 0.02]
                            qz_limits = [0.42, 0.5]
                            mag_limits = [0.0, 0.05]
                            hem_limits = [0.42, 0.5]
                            goe_limits = [0.0, 0.05]
                            kln_limits = [0.0, 0.02]
                            bt_limits = [0.0, 0.02]
                        elif classification == "Friable Quartz Itabirite":
                            gbs_limits = [0.0, 0.02]
                            qz_limits = [0.3, 0.38]
                            mag_limits = [0.0, 0.05]
                            hem_limits = [0.54, 0.62]
                            goe_limits = [0.0, 0.05]
                            kln_limits = [0.0, 0.02]
                            bt_limits = [0.0, 0.02]
                        elif classification == "Goethite Itabirite":
                            gbs_limits = [0.0, 0.05]
                            qz_limits = [0.18, 0.28]
                            mag_limits = [0.0, 0.05]
                            hem_limits = [0.32, 0.4]
                            goe_limits = [0.3, 0.38]
                            kln_limits = [0.0, 0.02]
                            bt_limits = [0.0, 0.02]
                        #
                        phi_hem = round(rd.uniform(hem_limits[0], hem_limits[1]), 4)
                        phi_goe = round(rd.uniform(goe_limits[0], (1 - phi_hem)), 4)
                        phi_kln = round(rd.uniform(kln_limits[0], (1 - phi_hem - phi_goe)), 4)
                        phi_qz = round(rd.uniform(qz_limits[0], (1 - phi_hem - phi_goe - phi_kln)), 4)
                        phi_mag = round(rd.uniform(mag_limits[0], (1 - phi_hem - phi_goe - phi_kln - phi_qz)), 4)
                        phi_gbs = round(rd.uniform(gbs_limits[0], (1 - phi_hem - phi_goe - phi_kln - phi_qz
                                                                   - phi_mag)), 4)
                        phi_bt = round(1 - phi_hem - phi_goe - phi_kln - phi_qz - phi_mag - phi_gbs, 4)
                        #
                        phi_total = phi_hem + phi_goe + phi_kln + phi_qz + phi_mag + phi_gbs + phi_bt
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if hem_limits[0] <= phi_hem <= hem_limits[1] \
                                    and goe_limits[0] <= phi_goe <= goe_limits[1] \
                                    and kln_limits[0] <= phi_kln <= kln_limits[1] \
                                    and qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and mag_limits[0] <= phi_mag <= mag_limits[1] \
                                    and gbs_limits[0] <= phi_gbs <= gbs_limits[1] \
                                    and bt_limits[0] <= phi_bt <= bt_limits[1]:
                                condition_2 = True
                        #
                    phi_minerals["Hem"] = phi_hem
                    phi_minerals["Goe"] = phi_goe
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Mag"] = phi_mag
                    phi_minerals["Gbs"] = phi_gbs
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
                        n_digits = 4
                    else:
                        n_digits = 4
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
            #
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

    def create_bandedironformation(self, rock="Banded Iron Formation", number=1, composition=None, mode="random"):
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
        w_last = None
        list_ore_values = np.around(np.linspace(0.25, 0.5, number), 4)
        list_ore_values_reverse = sorted(list_ore_values, reverse=True)
        while n < number:
            data_apatite = Phosphates(data_type=True).create_aptite()
            mineralogy = {
                "Qz": self.data_quartz, "Mag": self.data_magnetite, "Hem": self.data_hematite,
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
                    phi_mag = composition["Mag"]
                    phi_hem = composition["Hem"]
                    phi_kln = composition["Kln"]
                    phi_ap = composition["Ap"]

                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Mag"] = phi_mag
                    phi_minerals["Hem"] = phi_hem
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Ap"] = phi_ap

                else:
                    condition_2 = False
                    while condition_2 == False:
                        if mode == "random":
                            amount_ore = round(rd.uniform(0.25, 0.5), 4)
                        elif mode == "increasing ore":
                            amount_ore = list_ore_values[n]
                        elif mode == "decreasing ore":
                            amount_ore = list_ore_values_reverse[n]

                        amount_residuals = round(1 - amount_ore, 4)
                        amount_mag = round(rd.uniform(0.6, 0.8), 4)
                        amount_qz = round(rd.uniform(0.6, 1.0), 4)

                        qz_limits = [0.4, 0.75]
                        mag_limits = [0.15, 0.4]
                        hem_limits = [0.05, 0.20]
                        kln_limits = [0.0, 0.15]
                        ap_limits = [0.0, 0.025]
                        # Ore minerals
                        phi_mag = round(amount_ore*amount_mag, 4)
                        phi_hem = round(amount_ore - phi_mag, 4)
                        # Other minerals
                        phi_qz = round(amount_residuals*amount_qz, 4)
                        phi_ap = round(amount_residuals*rd.uniform(0, (1 - amount_qz)), 4)
                        phi_kln = round(1 - phi_qz - phi_mag - phi_hem - phi_ap, 4)

                        phi_total = phi_qz + phi_mag + phi_hem + phi_kln + phi_ap

                        if np.isclose(phi_total, 1.0000) == True:
                            if (qz_limits[0] <= phi_qz <= qz_limits[1] \
                                    and mag_limits[0] <= phi_mag <= mag_limits[1] \
                                    and hem_limits[0] <= phi_hem <= hem_limits[1] \
                                    and kln_limits[0] <= phi_kln <= kln_limits[1]
                                    and ap_limits[0] <= phi_ap <= ap_limits[1]):
                                condition_2 = True

                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Mag"] = phi_mag
                    phi_minerals["Hem"] = phi_hem
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
                rho_limits=[3000, 4500], vP_limits=[2500, 5000], vS_limits=[1500, 3000], delta=0.05,
                porosity=self.porosity)
            phi_neutron = round((1900/rho)*0.15, 4)
            ## Elastic Parameters
            bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                rho=rho, vP=vP, vS=vS)
            ## Gamma Ray
            gamma_ray = round(gamma_ray, 3)
            ## Photoelectricity
            photoelectricity = round(photoelectricity, 3)

            for key, value in w_minerals.items():
                results_container["mineralogy"][key].append(value)
            # for key, value in w_elements.items():
            #     results_container["chemistry"][key].append(value)
            # w_last = w_minerals

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

    def create_bauxite(self, rock="Bauxite", number=1, composition=None, mode="random"):
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

        mineralogy = {
            "Gbs": self.data_gibbsite, "Bhm": self.data_boehmite, "Goe": self.data_goethite, "Hem": self.data_hematite,
            "Kln": self.data_kaolinite, "Anat": self.data_anatase}
        minerals_list = list(mineralogy.keys())
        if minerals_list[0] not in results_container["mineralogy"]:
            for mineral in minerals_list:
                results_container["mineralogy"][mineral] = []

        n = 0
        w_last = None
        list_ore_values = np.around(np.linspace(0.25, 0.5, number), 4)
        list_ore_values_reverse = sorted(list_ore_values, reverse=True)
        while n < number:
            condition = False
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}

                if composition != None:
                    phi_gbs = composition["Gbs"]
                    phi_bhm = composition["Bhm"]
                    phi_goe = composition["Goe"]
                    phi_hem = composition["Hem"]
                    phi_kln = composition["Kln"]
                    phi_anat = composition["Anat"]

                    phi_minerals["Gbs"] = phi_gbs
                    phi_minerals["Bhm"] = phi_bhm
                    phi_minerals["Goe"] = phi_goe
                    phi_minerals["Hem"] = phi_hem
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Anat"] = phi_anat

                else:
                    condition_2 = False
                    while condition_2 == False:
                        if mode == "random":
                            amount_al_ore = round(rd.uniform(0.35, 0.65), 4)
                        elif mode == "increasing ore":
                            amount_al_ore = list_ore_values[n]
                        elif mode == "decreasing ore":
                            amount_al_ore = list_ore_values_reverse[n]

                        amount_fe_ore = round(rd.uniform(0.10, (1 - amount_al_ore)), 4)
                        amount_residuals = round(1 - amount_al_ore - amount_fe_ore, 4)
                        amount_gbs = round(rd.uniform(0.2, 0.8), 4)
                        amount_goe = round(rd.uniform(0.7, 0.9), 4)
                        amount_kln = round(rd.uniform(0.4, 0.9), 4)

                        # Al ore minerals
                        phi_gbs = round(amount_al_ore*amount_gbs, 4)
                        phi_bhm = round(amount_al_ore - phi_gbs, 4)
                        # Fe ore minerals
                        phi_goe = round(amount_fe_ore*amount_goe, 4)
                        phi_hem = round(amount_fe_ore - phi_goe, 4)
                        # Other minerals

                        phi_kln = round(amount_residuals*amount_kln, 4)
                        phi_anat = round(1 - phi_gbs - phi_bhm - phi_goe - phi_hem - phi_kln, 4)

                        phi_total = phi_gbs + phi_bhm + phi_goe + phi_hem + phi_kln + phi_anat

                        gbs_limits = [0.2, 0.8]
                        bhm_limits = [0.2, 0.6]
                        goe_limits = [0.2, 0.5]
                        hem_limits = [0.0, 0.25]
                        kln_limits = [0.0, 0.3]
                        anat_limits = [0.0, 0.10]

                        if np.isclose(phi_total, 1.0000) == True:
                            if (gbs_limits[0] <= phi_gbs <= gbs_limits[1] \
                                    and bhm_limits[0] <= phi_bhm <= bhm_limits[1] \
                                    and goe_limits[0] <= phi_goe <= goe_limits[1] \
                                    and hem_limits[0] <= phi_hem <= hem_limits[1] \
                                    and kln_limits[0] <= phi_kln <= kln_limits[1] \
                                    and anat_limits[0] <= phi_kln <= anat_limits[1]):
                                condition_2 = True

                    phi_minerals["Gbs"] = phi_gbs
                    phi_minerals["Bhm"] = phi_bhm
                    phi_minerals["Goe"] = phi_goe
                    phi_minerals["Hem"] = phi_hem
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Anat"] = phi_anat

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
                rho_limits=[2500, 4000], vP_limits=[3000, 5000], vS_limits=[1500, 3000], delta=0.05,
                porosity=self.porosity)
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

class Ores:
    #
    def __init__(self, fluid, actualThickness, porosity=None, data_type=False):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.porosity = porosity
        self.data_type = data_type
        self.data_water = Water.water("")
    #
    def create_kupferschiefer(self, w_Cu=None, amounts=None):
        #
        results = {}
        results["rock"] = "Kupferschiefer"
        #
        self.w_Cu = w_Cu
        self.amounts = amounts
        #
        # Mineralogy
        quartz = Oxides(data_type=True).create_quartz()
        calcite = Carbonates(data_type=True).create_calcite()
        pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()
        chalcopyrite = Sulfides(impurity="pure", data_type=True).create_chalcopyrite()
        galena = Sulfides(impurity="pure", data_type=True).create_galena()
        bornite = Sulfides(impurity="pure", data_type=True).create_bornite()
        sphalerite = Sulfides(impurity="pure", data_type=True).create_sphalerite()
        chalcocite = Sulfides(impurity="pure", data_type=True).create_chalcocite()
        covellite = Sulfides(impurity="pure", data_type=True).create_covellite()
        digenite = Sulfides(impurity="pure", data_type=True).create_digenite()
        illite = Pyllosilicates(impurity="pure", dict=True).create_illite()
        kaolinite = Pyllosilicates(impurity="pure", dict=True).create_kaolinite()
        montmorillonite = Pyllosilicates(impurity="pure", dict=True).create_montmorillonite()
        #
        mineralogy = [illite, kaolinite, montmorillonite, quartz, calcite, chalcopyrite, bornite, chalcocite, covellite,
                      digenite, galena, sphalerite, pyrite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Cu == None and self.amounts == None:
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ore = round(rd.uniform(0.4, (1-w_clay)), 4)
                w_misc = round(rd.uniform(0.0, (1-w_clay-w_ore)), 4)
                #
                w_ilt = round(w_clay*rd.uniform(0.5, 1), 4)
                w_kln = round(w_clay*rd.uniform(0, (1-w_ilt)), 4)
                w_mnt = round(w_clay-w_ilt-w_kln, 4)
                #
                magicnumber = rd.randint(1, 7)
                if magicnumber == 1:
                    w_cv = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_bn = round(w_ore*rd.uniform(0, (1-w_cv)), 4)
                    w_cc = round(w_ore*rd.uniform(0, (1-w_cv-w_bn)), 4)
                    w_ccp = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc)), 4)
                    w_dg = round(w_ore-w_cv-w_bn-w_cc-w_ccp, 4)
                    #
                    w_gn = 0.0
                    w_py = 0.0
                    w_sp = 0.0
                elif magicnumber == 2:
                    w_bn = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_ccp = round(w_ore*rd.uniform(0, (1-w_bn)), 4)
                    w_sp = round(w_ore*rd.uniform(0, (1-w_bn-w_ccp)), 4)
                    w_dg = round(w_ore-w_bn-w_ccp-w_sp, 4)
                    #
                    w_cv = 0.0
                    w_cc = 0.0
                    w_gn = 0.0
                    w_py = 0.0
                    w_dg = 0.0
                elif magicnumber == 3:
                    w_cc = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_gn = round(w_ore*rd.uniform(0, (1-w_cc)), 4)
                    w_sp = round(w_ore-w_cc-w_gn, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_ccp = 0.0
                    w_py = 0.0
                    w_dg = 0.0
                elif magicnumber == 4:
                    w_gn = round(w_ore*rd.uniform(0.33, 0.67), 4)
                    w_ccp = round(w_ore*rd.uniform(0.05, (1-w_gn)), 4)
                    w_py = round(w_ore*rd.uniform(0, (1-w_gn-w_ccp)), 4)
                    w_sp = round(w_ore-w_gn-w_ccp-w_py, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_cc = 0.0
                    w_dg = 0.0
                elif magicnumber == 5:
                    w_py = round(w_ore*rd.uniform(0.33, 0.67), 4)
                    w_gn = round(w_ore*rd.uniform(0, (1-w_py)), 4)
                    w_sp = round(w_ore*rd.uniform(0, (1-w_py-w_gn)), 4)
                    w_ccp = round(w_ore-w_py-w_gn-w_sp, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_cc = 0.0
                    w_dg = 0.0
                elif magicnumber == 6:
                    w_sp = round(w_ore*rd.uniform(0.33, 0.67), 4)
                    w_ccp = round(w_ore*rd.uniform(0.05, (1-w_sp)), 4)
                    w_gn = round(w_ore*rd.uniform(0, (1-w_sp-w_ccp)), 4)
                    w_py = round(w_ore-w_sp-w_ccp-w_gn, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_cc = 0.0
                    w_dg = 0.0
                elif magicnumber == 7:
                    w_dg = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_cv = round(w_ore*rd.uniform(0, (1-w_dg)), 4)
                    w_bn = round(w_ore*rd.uniform(0, (1-w_dg-w_cv)), 4)
                    w_cc = round(w_ore*rd.uniform(0, (1-w_dg-w_cv-w_bn)), 4)
                    w_ccp = round(w_ore-w_dg-w_cv-w_bn-w_cc, 4)
                    #
                    w_gn = 0.0
                    w_py = 0.0
                    w_sp = 0.0
                #
                w_qz = round(w_misc*rd.uniform(0, 1), 4)
                w_cal = round(1 - w_clay - w_ore - w_qz, 4)
            elif self.w_Cu != None:
                w_clay = round(rd.uniform(0.33, 0.67), 4)
                w_ore = round(rd.uniform(0.33, (1-w_clay)), 4)
                w_misc = round(rd.uniform(0.0, (1-w_clay-w_ore)), 4)
                #
                w_ilt = round(w_clay*rd.uniform(0, 1), 4)
                w_kln = round(w_clay*rd.uniform(0, (1-w_ilt)), 4)
                w_mnt = round(w_clay-w_ilt-w_kln, 4)
                #
                w_cv = round(w_ore*rd.uniform(0, 1), 4)
                w_bn = round(w_ore*rd.uniform(0, (1-w_cv)), 4)
                w_cc = round(w_ore*rd.uniform(0, (1-w_cv-w_bn)), 4)
                w_ccp = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc)), 4)
                w_gn = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc-w_ccp)), 4)
                w_py = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc-w_ccp-w_gn)), 4)
                w_dg = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc-w_ccp-w_gn-w_py)), 4)
                w_sp = round(w_ore-w_cv-w_bn-w_cc-w_ccp-w_gn-w_py-w_dg, 4)
                #
                w_qz = round(w_misc*rd.uniform(0, 1), 4)
                w_cal = 1 - w_clay - w_ore - w_qz
            elif type(self.amounts) is list:
                w_hl = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_anh = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_gp = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_syl = round(1-w_hl-w_anh-w_gp, 4)
            #
            if w_ilt >= 0.0 and w_kln >= 0.0 and w_mnt >= 0.0 and w_cv >= 0.0 and w_bn >= 0.0 and w_cc >= 0.0 \
                    and w_ccp >= 0.0 and w_gn >= 0.0 and w_py >= 0.0 and w_dg >= 0.0 and w_sp >= 0.0 and w_qz >= 0.0 \
                    and w_cal >= 0.0:
                sumMin = round(w_ilt + w_kln + w_mnt + w_cv + w_bn + w_cc + w_ccp + w_gn + w_py + w_dg + w_sp + w_qz + w_cal, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ilt*illite["chemistry"]["H"] + w_kln*kaolinite["chemistry"]["H"] + w_mnt*montmorillonite["chemistry"]["H"], 4)
            w_C = round(w_cal*calcite["chemistry"]["C"], 4)
            w_Na = round(w_mnt*montmorillonite["chemistry"]["Na"], 4)
            w_Mg = round(w_ilt*illite["chemistry"]["Mg"] + w_mnt*montmorillonite["chemistry"]["Mg"], 4)
            w_Al = round(w_ilt*illite["chemistry"]["Al"] + w_kln*kaolinite["chemistry"]["Al"] + w_mnt*montmorillonite["chemistry"]["Al"], 4)
            w_Si = round(w_ilt*illite["chemistry"]["Si"] + w_kln*kaolinite["chemistry"]["Si"] + w_mnt*montmorillonite["chemistry"]["Si"] + w_qz*quartz["chemistry"]["Si"], 4)
            w_S = round(w_cv*covellite["chemistry"]["S"] + w_bn*bornite["chemistry"]["S"] + w_cc*chalcocite["chemistry"]["S"] + w_ccp*chalcopyrite["chemistry"]["S"] + w_gn*galena["chemistry"]["S"] + w_py*pyrite["chemistry"]["S"] + w_sp*sphalerite["chemistry"]["S"] + w_dg*digenite["chemistry"]["S"], 4)
            w_K = round(w_ilt*illite["chemistry"]["K"], 4)
            w_Ca = round(w_mnt*montmorillonite["chemistry"]["Ca"] + w_cal*calcite["chemistry"]["Ca"], 4)
            w_Fe = round(w_ilt*illite["chemistry"]["Fe"] + w_bn*bornite["chemistry"]["Fe"] + w_ccp*chalcopyrite["chemistry"]["Fe"] + w_py*pyrite["chemistry"]["Fe"], 4)
            w_Cu = round(w_cv*covellite["chemistry"]["Cu"] + w_bn*bornite["chemistry"]["Cu"] + w_cc*chalcocite["chemistry"]["Cu"] + w_ccp*chalcopyrite["chemistry"]["Cu"] + w_dg*digenite["chemistry"]["Cu"], 4)
            w_Zn = round(w_sp*sphalerite["chemistry"]["Zn"], 4)
            w_Pb = round(w_gn*galena["chemistry"]["Pb"], 4)
            w_O = round(1 - w_H - w_C - w_Na - w_Mg - w_Al - w_Si - w_S - w_K - w_Ca - w_Fe - w_Cu - w_Zn - w_Pb, 4)
            sumConc = round(w_H + w_C + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe + w_Cu + w_Zn + w_Pb, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Ilt", "Kln", "Mnt", "Qz", "Cal", "Ccp", "Bn", "Cc", "Cv", "Dg", "Gn", "Sp", "Py"]))
                concentrations = [w_H, w_C, w_O, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe, w_Cu, w_Zn, w_Pb]
                amounts = [w_ilt, w_kln, w_mnt, w_qz, w_cal, w_ccp, w_bn, w_cc, w_cv, w_dg, w_gn, w_sp, w_py]
            else:
                cond = False
        #
        element_list = ["H", "C", "O", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Fe", "Cu", "Zn", "Pb"]
        mineral_list = ["Ilt", "Kln", "Mnt", "Qz", "Cal", "Ccp", "Bn", "Cc", "Cv", "Dg", "Gn", "Sp", "Py"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_ilt*illite["rho"] + w_kln*kaolinite["rho"] + w_mnt*montmorillonite["rho"] + w_qz*quartz["rho"]
                    + w_cal*calcite["rho"] + w_ccp*chalcopyrite["rho"] + w_bn*bornite["rho"] + w_cc*chalcocite["rho"]
                    + w_cv*covellite["rho"] + w_dg*digenite["rho"] + w_gn*galena["rho"] + w_sp*sphalerite["rho"]
                    + w_py*pyrite["rho"]) / 1000
        X = [w_ilt, w_kln, w_mnt, w_qz, w_cal, w_ccp, w_bn, w_cc, w_cv, w_dg, w_gn, w_sp, w_py]
        K_list = [mineralogy[i]["K"] for i in range(len(mineralogy))]
        G_list = [mineralogy[i]["G"] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        #
        if self.porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.uniform(0.0, 0.1)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.uniform(0.0, 0.075)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.uniform(0.0, 0.05)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 4000:
                phi = rd.uniform(0.0, 0.0125)
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
        GR = w_ilt*illite["GR"] + w_kln*kaolinite["GR"] + w_mnt*montmorillonite["GR"] + w_qz*quartz["GR"] \
             + w_cal*calcite["GR"] + w_ccp*chalcopyrite["GR"] + w_bn*bornite["GR"] + w_cc*chalcocite["GR"] \
             + w_cv*covellite["GR"] + w_dg*digenite["GR"] + w_gn*galena["GR"] + w_sp*sphalerite["GR"] + w_py*pyrite["GR"]
        PE = w_ilt*illite["PE"] + w_kln*kaolinite["PE"] + w_mnt*montmorillonite["PE"] + w_qz*quartz["PE"] \
             + w_cal*calcite["PE"] + w_ccp*chalcopyrite["PE"] + w_bn*bornite["PE"] + w_cc*chalcocite["PE"] \
             + w_cv*covellite["PE"] + w_dg*digenite["PE"] + w_gn*galena["PE"] + w_sp*sphalerite["PE"] + w_py*pyrite["PE"]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_ilt*illite["nu"] + w_kln*kaolinite["nu"] + w_mnt*montmorillonite["nu"] + w_qz*quartz["nu"] \
                                + w_cal*calcite["nu"] + w_ccp*chalcopyrite["nu"] + w_bn*bornite["nu"] + w_cc*chalcocite["nu"] \
                                + w_cv*covellite["nu"] + w_dg*digenite["nu"] + w_gn*galena["nu"] + w_sp*sphalerite["nu"] + w_py*pyrite["nu"]
        #
        if self.data_type == False:
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
        else:
            #
            results["rho"] = round(rho*1000, 4)
            results["rho_s"] = round(rhoSolid*1000, 4)
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
    def create_compact_hematite_ore(self, number, porosity=None):
        #
        data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()           # fixed
        data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()       # fixed
        data_magnetite = Oxides(impurity="pure", data_type=True).create_magnetite()     # fixed
        data_pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()         # fixed
        data_rutile = Oxides(impurity="pure", data_type=True).create_rutile()           # fixed
        #
        assemblage_minerals = [data_quartz, data_hematite, data_magnetite, data_pyrite, data_rutile]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage_minerals:
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
                        value = round(rd.uniform(0.10, 0.35), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.35:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.45, 0.80), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.45 <= value <= 0.80:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Mag":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.10, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Py":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.10), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.10:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Rt":
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
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
            amounts_helper.clear()
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
            for element in elements:
                amounts_helper[element] = 0
                if element in data_quartz["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Qz"][n]*data_quartz["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_hematite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Hem"][n]*data_hematite["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_magnetite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Mag"][n]*data_magnetite["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_pyrite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Py"][n]*data_pyrite["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_rutile["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Rt"][n]*data_rutile["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                #
                n_elements += 1
            #
            value_total = 0
            for element, value in amounts_helper.items():
                amounts_helper[element] = round(value, 4)
                value_total += int(round(amounts_helper[element]*10000, 4))
            value_total = int(value_total)
            if value_total == 10000:
                condition_sum = True
            else:
                condition_sum = False
            if sum(amounts_helper.values()) == 1.0 or condition_sum == True or w_total == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_quartz["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_quartz["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_hematite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_hematite["PE"], 3)
                    elif mineral == "Mag":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_magnetite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["PE"], 3)
                    elif mineral == "Py":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_pyrite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["PE"], 3)
                    elif mineral == "Rt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_rutile["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_rutile["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_rutile["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_rutile["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_rutile["PE"], 3)
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
            #else:
            #    break
        #
        results = {}
        results["rock"] = "Compact Hematite Ore"
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
    def create_banded_iron_formation(self, number, porosity=None):
        #
        data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()  # fixed
        data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()  # fixed
        data_magnetite = Oxides(impurity="pure", data_type=True).create_magnetite()  # fixed
        data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()  # fixed
        data_goethite = Oxides(impurity="pure", data_type=True).create_goethite()  # fixed
        #
        assemblage_minerals = [data_quartz, data_hematite, data_magnetite, data_kaolinite, data_goethite]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage_minerals:
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
                    if value >= 0.0 and 0.10 <= value <= 0.35:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.25, 0.70), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.70:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Mag":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.10, 0.25), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kln":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.10, 0.25), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.25:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Goe":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.05, 0.25), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.25:
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
            amounts_helper.clear()
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
            for element in elements:
                amounts_helper[element] = 0
                if element in data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Qz"][n] * data_quartz["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_hematite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Hem"][n] * data_hematite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_magnetite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Mag"][n] * data_magnetite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_kaolinite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Kln"][n] * data_kaolinite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_goethite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Goe"][n] * data_goethite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                #
                n_elements += 1
            #
            value_total = 0
            for element, value in amounts_helper.items():
                amounts_helper[element] = round(value, 4)
                value_total += int(round(amounts_helper[element] * 10000, 4))
            value_total = int(value_total)
            if value_total == 10000:
                condition_sum = True
            else:
                condition_sum = False
            if sum(amounts_helper.values()) == 1.0 or condition_sum == True or w_total == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_quartz["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_quartz["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_hematite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_hematite["PE"], 3)
                    elif mineral == "Mag":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_magnetite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["PE"], 3)
                    elif mineral == "Kln":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_kaolinite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["PE"], 3)
                    elif mineral == "Goe":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_goethite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_goethite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_goethite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_goethite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_goethite["PE"], 3)
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
            #else:
            #    break
        #
        results = {}
        results["rock"] = "Banded Iron Formation"
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
# print("Test: Compact Hematite Ore")
# start = time.process_time()
# data = []
# for i in range(100):
#     data.append(Ores(fluid="water", actualThickness=0).create_compact_hematite_ore(number=1))
# # print(Ores(fluid="water", actualThickness=0).create_compact_hematite_ore(number=1))
# print(data)
# end1 = time.process_time() - start
# print(end1)
# start = time.process_time()
# print(Ores(fluid="water", actualThickness=0).create_compact_hematite_ore(number=100))
# end2 = time.process_time() - start
# print(end2)
# print(end2/end1)