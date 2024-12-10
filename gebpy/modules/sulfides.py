#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		sulfides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		09.10.2023

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry, TraceElements

# Sulfides
class Sulfides():
    """ Class that generates geophysical and geochemical data of sulfide minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral

        # Chemistry
        self.boron = ["B", 5, 10.806]
        self.carbon = ["C", 6, 12.011]
        self.oxygen = ["O", 8, 15.999]
        self.sodium = ["Na", 11, 22.990]
        self.aluminium = ["Al", 13, 26.982]
        self.silicon = ["Si", 14, 28.085]
        self.sulfur = ["S", 16, 32.059]
        self.chlorine = ["Cl", 17, 35.450]
        self.potassium = ["K", 19, 39.098]
        self.calcium = ["Ca", 20, 40.078]
        self.manganese = ["Mn", 25, 54.938]
        self.iron = ["Fe", 26, 55.845]
        self.copper = ["Cu", 29, 63.546]
        self.arsenic = ["As", 33, 74.922]
        self.antimony = ["Sb", 51, 121.76]
        self.tungsten = ["W", 74, 183.84]

    def get_data(self, number=1): # ["Pyrite", "Chalcopyrite", "Bornite", "Covellite", "Molybdenite", "Sphalerite", "Galena", "Fahlore"]
        if self.mineral in ["Py", "Pyrite"]:
            if number > 1:
                data = [self.create_pyrite() for n in range(number)]
            else:
                data = self.create_pyrite()
        elif self.mineral in ["Ccp", "Chalcopyrite"]:
            if number > 1:
                data = [self.create_chalcopyrite() for n in range(number)]
            else:
                data = self.create_chalcopyrite()
        elif self.mineral in ["Bn", "Bornite"]:
            if number > 1:
                data = [self.create_bornite() for n in range(number)]
            else:
                data = self.create_bornite()
        elif self.mineral in ["Cv", "Covellite"]:
            if number > 1:
                data = [self.create_covellite() for n in range(number)]
            else:
                data = self.create_covellite()
        elif self.mineral in ["Mol", "Molybdenite"]:
            if number > 1:
                data = [self.create_molybdenite() for n in range(number)]
            else:
                data = self.create_molybdenite()
        elif self.mineral in ["Sp", "Sphalerite"]:
            if number > 1:
                data = [self.create_sphalerite() for n in range(number)]
            else:
                data = self.create_sphalerite()
        elif self.mineral in ["Gn", "Galena"]:
            if number > 1:
                data = [self.create_galena() for n in range(number)]
            else:
                data = self.create_galena()
        elif self.mineral in ["Fh", "Fahlore"]:
            if number > 1:
                data = [self.create_fahlore() for n in range(number)]
            else:
                data = self.create_fahlore()
        elif self.mineral in ["Cbt", "Cobaltite"]:
            if number > 1:
                data = [self.create_cobaltite() for n in range(number)]
            else:
                data = self.create_cobaltite()
        elif self.mineral in ["Marmatite"]:
            if number > 1:
                data = [self.create_marmatite() for n in range(number)]
            else:
                data = self.create_marmatite()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Pyrite":
                data_mineral = self.create_pyrite()
            elif self.mineral == "Chalcopyrite":
                data_mineral = self.create_chalcopyrite()
            elif self.mineral == "Galena":
                data_mineral = self.create_galena()
            elif self.mineral == "Acanthite":
                data_mineral = self.create_acanthite()
            elif self.mineral == "Chalcocite":
                data_mineral = self.create_chalcocite()
            elif self.mineral == "Bornite":
                data_mineral = self.create_bornite()
            elif self.mineral == "Pyrrhotite":
                data_mineral = self.create_pyrrhotite()
            elif self.mineral == "Millerite":
                data_mineral = self.create_millerite()
            #
            # Sphalerite Group
            elif self.mineral == "Sphalerite":
                data_mineral = self.create_sphalerite()
            elif self.mineral == "Marmatite":
                data_mineral = self.create_marmatite()
            #
            elif self.mineral == "Pentlandite":
                data_mineral = self.create_pentlandite()
            elif self.mineral == "Covellite":
                data_mineral = self.create_covellite()
            elif self.mineral == "Cinnabar":
                data_mineral = self.create_cinnabar()
            elif self.mineral == "Realgar":
                data_mineral = self.create_realgar()
            elif self.mineral == "Orpiment":
                data_mineral = self.create_orpiment()
            elif self.mineral == "Stibnite":
                data_mineral = self.create_stibnite()
            elif self.mineral == "Marcasite":
                data_mineral = self.create_marcasite()
            elif self.mineral == "Molybdenite":
                data_mineral = self.create_molybdenite()
            elif self.mineral == "Fahlore":
                data_mineral = self.create_fahlore()
            elif self.mineral == "Gallite":
                data_mineral = self.create_gallite()
            elif self.mineral == "Roquesite":
                data_mineral = self.create_roquesite()
            elif self.mineral == "Lenaite":
                data_mineral = self.create_lenaite()
            elif self.mineral == "Laforetite":
                data_mineral = self.create_laforetite()
            elif self.mineral == "Vaesite":
                data_mineral = self.create_vaesite()
            elif self.mineral == "Cattierite":
                data_mineral = self.create_cattierite()
            elif self.mineral == "Cobaltite":
                data_mineral = self.create_cobaltite()
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
    def create_cinnabar(self):
        #
        name = "Ci"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        mercury = PeriodicSystem(name="Hg").get_data()
        majors_name = ["S", "Hg"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Hg", mercury[1], 1, mercury[2]]], dtype=object)
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
        molar_mass_pure = mercury[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.149, 9.495], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 3, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 8*10**9
        # Shear modulus
        G = 7*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -45734     # J/mol
        thermodynamics["Enthalpy"] = -53346         # J/mol
        thermodynamics["Entropy"] = 82.425          # J/(mol K)
        thermodynamics["Heat Capacity"] = 48.41     # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_pyrite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["S", "Fe"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
        # Trace elements
        elements_traces = ["Ni", "Co", "As", "Cu", "Zn", "Ag", "Au", "Tl", "Se", "V"]
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ni", "Co", "As", "Cu", "Zn", "Ag", "Au", "Tl", "Se", "V"]
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
        mineral = "Py"
        #
        # Molar mass
        molar_mass_pure = iron[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.417], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 146*10**9
        # Shear modulus
        G = 135*10**9
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
        p = 3*10**(-1)
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -159572     # J/mol
        thermodynamics["Enthalpy"] = -171050         # J/mol
        thermodynamics["Entropy"] = 52.930         # J/(mol K)
        thermodynamics["Heat Capacity"] = 62.19    # J/(mol K)
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
            results["mineral"] = "Py"
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            results["trace elements"] = elements_traces
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_bornite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["S", "Fe", "Cu"]
        majors_data = np.array([["S", sulfur[1], 4, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["Cu", copper[1], 5, copper[2]]], dtype=object)
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
                minors = ["Ag", "Ge", "Bi", "In", "Pb"]
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
        mineral = "Bn"
        #
        # Molar mass
        molar_mass_pure = 5*copper[2] + iron[2] + 4*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.95, 21.862, 10.95], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 16, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 92.29*10**9
        # Shear modulus
        G = 36.85*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -394621    # J/mol
        thermodynamics["Enthalpy"] = -371600        # J/mol
        thermodynamics["Entropy"] = 398.500         # J/(mol K)
        thermodynamics["Heat Capacity"] = 245.60    # J/(mol K)
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
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_galena(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["S", "Pb"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Pb", lead[1], 1, lead[2]]], dtype=object)
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
                minors = ["Ag", "Cu", "Fe", "Bi"]
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
        mineral = "Gn"
        #
        # Molar mass
        molar_mass_pure = lead[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.936], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 52*10**9
        # Shear modulus
        G = 29*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -96676     # J/mol
        thermodynamics["Enthalpy"] = -98320         # J/mol
        thermodynamics["Entropy"] = 91.340          # J/(mol K)
        thermodynamics["Heat Capacity"] = 49.44     # J/(mol K)
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
            results["mineral"] = "Gn"
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_chalcopyrite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["S", "Fe", "Cu"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["Cu", copper[1], 1, copper[2]]], dtype=object)
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
                minors = ["Ag", "Au", "In", "Tl", "Se", "Te"]
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
        mineral = "Ccp"
        #
        # Molar mass
        molar_mass_pure = copper[2] + iron[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.28, 10.41], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 56*10**9
        # Shear modulus
        G = 19*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -194996    # J/mol
        thermodynamics["Enthalpy"] = -194900        # J/mol
        thermodynamics["Entropy"] = 124.900         # J/(mol K)
        thermodynamics["Heat Capacity"] = 96.65     # J/(mol K)
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
            results["mineral"] = "Ccp"
            results["state"] = var_state
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_molybdenite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        molybdenum = PeriodicSystem(name="Mo").get_data()
        majors_name = ["S", "Mo"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Mo", molybdenum[1], 1, molybdenum[2]]], dtype=object)
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
        name = "Mol"
        #
        # Molar mass
        molar_mass_pure = molybdenum[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.16, 12.3], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 2, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 25*10**9
        # Shear modulus
        G = 17*10**9
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
            results["V"] = round(V, 4)
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
    def create_sphalerite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["S", "Zn"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Zn", zinc[1], 1, zinc[2]]], dtype=object)
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
                minors = ["Mn", "Cd", "Hg", "In", "Tl", "Ga", "Ge", "Sb", "Sn", "Pb", "Ag", "Co"]
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
        mineral = "Sp"
        #
        # Molar mass
        molar_mass_pure = zinc[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.406], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 68*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -198600     # J/mol
        thermodynamics["Enthalpy"] = -203107         # J/mol
        thermodynamics["Entropy"] = 58.655         # J/(mol K)
        thermodynamics["Heat Capacity"] = 45.76    # J/(mol K)
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
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_stibnite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        antimony = PeriodicSystem(name="Sb").get_data()
        majors_name = ["S", "Sb"]
        majors_data = np.array([["S", sulfur[1], 3, sulfur[2]], ["Sb", antimony[1], 2, antimony[2]]], dtype=object)
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
        name = "Stbn"
        #
        # Molar mass
        molar_mass_pure = 2*antimony[2] + 3*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.229,11.31, 3.893], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 16*10**9
        # Shear modulus
        G = 9*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -149980     # J/mol
        thermodynamics["Enthalpy"] = -151529         # J/mol
        thermodynamics["Entropy"] = 182.004         # J/(mol K)
        thermodynamics["Heat Capacity"] = 225.19    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_arsenopyrite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["S", "Fe", "As"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["As", arsenic[1], 1, arsenic[2]]], dtype=object)
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
                minors = ["Ag", "Au", "Co", "Sn", "Ni", "Sb", "Bi", "Cu", "Pb"]
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
        name = "Apy"
        #
        # Molar mass
        molar_mass_pure = iron[2] + arsenic[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.74, 5.68, 5.79], [112.17], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 143*10**9
        # Shear modulus
        G = 117*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -136450   # J/mol
        thermodynamics["Enthalpy"] = -144370       # J/mol
        thermodynamics["Entropy"] = 68.500         # J/(mol K)
        thermodynamics["Heat Capacity"] = 68.45    # J/(mol K)
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 2))
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            results = {}
            results["mineral"] = name
            results["M"] = molar_mass
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            try:
                results["p"] = round(p, 4)
            except:
                results["p"] = p
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            #
            return results
    #
    def create_acanthite(self):
        #
        name = "Ach"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        silver = PeriodicSystem(name="Ag").get_data()
        majors_name = ["S", "Ag"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Ag", silver[1], 2, silver[2]]], dtype=object)
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
                minors = ["Se"]
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
        molar_mass_pure = 2*silver[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.229, 6.931, 7.862], [99.61], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 62.48*10**9
        # Shear modulus
        G = 22.59*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -39448   # J/mol
        thermodynamics["Enthalpy"] = -31589       # J/mol
        thermodynamics["Entropy"] = 143.511       # J/(mol K)
        thermodynamics["Heat Capacity"] = 76.12   # J/(mol K)
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 2))
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            results = {}
            results["mineral"] = name
            results["M"] = molar_mass
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            try:
                results["p"] = round(p, 4)
            except:
                results["p"] = p
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            #
            return results
    #
    def create_argentite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        silver = PeriodicSystem(name="Ag").get_data()
        majors_name = ["S", "Ag"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Ag", silver[1], 2, silver[2]]], dtype=object)
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
                minors = ["Se"]
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
        mineral = "Argt"
        #
        # Molar mass
        molar_mass_pure = 2*silver[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.89], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 63.37*10**9
        # Shear modulus
        G = 22.61*10**9
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
        data.append(mineral)
        data.append(round(molar_mass, 3))
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_alabandite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["S", "Mn"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
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
                minors = ["Fe", "Mg", "Co"]
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
        name = "Ab"
        #
        # Molar mass
        molar_mass_pure = manganese[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.214], [99.61], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 56*10**9
        # Shear modulus
        G = 39*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -218313   # J/mol
        thermodynamics["Enthalpy"] = -213462       # J/mol
        thermodynamics["Entropy"] = 80.333         # J/(mol K)
        thermodynamics["Heat Capacity"] = 49.94    # J/(mol K)
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 2))
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            results = {}
            results["mineral"] = name
            results["M"] = molar_mass
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            try:
                results["p"] = round(p, 4)
            except:
                results["p"] = p
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            #
            return results
    #
    def create_berthierite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        antimony = PeriodicSystem(name="Sb").get_data()
        majors_name = ["S", "Fe", "Sb"]
        majors_data = np.array([["S", sulfur[1], 4, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["Sb", antimony[1], 2, antimony[2]]], dtype=object)
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
                minors = ["Mn", "Cu", "Pb", "Ag"]
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
        name = "Brt"
        #
        # Molar mass
        molar_mass_pure = iron[2] + 2*antimony[2] + 4*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.44, 14.12, 3.76], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 60.34*10**9
        # Shear modulus
        G = 24.92*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -255960   # J/mol
        thermodynamics["Enthalpy"] = -256426       # J/mol
        thermodynamics["Entropy"] = 245.015         # J/(mol K)
        thermodynamics["Heat Capacity"] = 177.09    # J/(mol K)
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 2))
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            results = {}
            results["mineral"] = name
            results["M"] = molar_mass
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            try:
                results["p"] = round(p, 4)
            except:
                results["p"] = p
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            #
            return results
    #
    def create_pyrrhotite_group(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["S", "Fe"]
        n_Fe = [1, 7, 11]
        n_S = [1, 8, 12]
        i = rd.randint(0, 2)
        x = round(rd.uniform(0, 1), 4)
        n_Fe = 1 + 10*x
        n_S = 1 + 11*x
        majors_data = np.array([["S", sulfur[1], n_S, sulfur[2]], ["Fe", iron[1], n_Fe, iron[2]]], dtype=object)
        majors_data_low = np.array([["S", sulfur[1], 1, sulfur[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
        majors_data_med = np.array([["S", sulfur[1], 8, sulfur[2]], ["Fe", iron[1], 7, iron[2]]], dtype=object)
        majors_data_high = np.array([["S", sulfur[1], 12, sulfur[2]], ["Fe", iron[1], 11, iron[2]]], dtype=object)
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
                minors = ["Ni", "Co", "Cu"]
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
        mineral = "Po"
        #
        # Molar mass
        molar_mass_pure = n_Fe*iron[2] + n_S*sulfur[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        molar_mass_pure_low = 1*iron[2] + 1*sulfur[2]
        molar_mass_low, amounts_low = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_low, majors=majors_data_low).calculate_molar_mass()
        element_low = [PeriodicSystem(name=amounts_low[i][0]).get_data() for i in range(len(amounts_low))]
        #
        molar_mass_pure_med = 7*iron[2] + 8*sulfur[2]
        molar_mass_med, amounts_med = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_med, majors=majors_data_med).calculate_molar_mass()
        element_med = [PeriodicSystem(name=amounts_med[i][0]).get_data() for i in range(len(amounts_med))]
        #
        molar_mass_pure_high = 11*iron[2] + 12*sulfur[2]
        molar_mass_high, amounts_high = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_high, majors=majors_data_high).calculate_molar_mass()
        element_high = [PeriodicSystem(name=amounts_high[i][0]).get_data() for i in range(len(amounts_high))]
        # Density
        dataV_low = CrystalPhysics([[3.407, 5.727], [], "hexagonal"])
        V_low = dataV_low.calculate_volume()
        dataRho_low = CrystalPhysics([molar_mass_pure_low, 2, V_low*10**(6)])
        rho_low = dataRho_low.calculate_bulk_density()
        rho_e_low = wg(amounts=amounts_low, elements=element_low, rho_b=rho_low).calculate_electron_density()
        #
        dataV_med = CrystalPhysics([[5.740, 6.151, 11.172], [83.295], "monoclinic"])
        V_med = dataV_med.calculate_volume()
        dataRho_med = CrystalPhysics([molar_mass_pure_med, 2, V_med])
        rho_med = dataRho_med.calculate_bulk_density()
        rho_e_med = wg(amounts=amounts_med, elements=element_med, rho_b=rho_med).calculate_electron_density()
        #
        dataV_high = CrystalPhysics([[6.624, 6.579, 16.231], [57.314], "monoclinic"])
        V_high = dataV_high.calculate_volume()
        dataRho_high = CrystalPhysics([molar_mass_pure_high, 2, V_high])
        rho_high = dataRho_high.calculate_bulk_density()
        rho_e_high = wg(amounts=amounts_high, elements=element_high, rho_b=rho_high).calculate_electron_density()
        #
        V = -65.31*(2*x + 1)**2 + 530.08*(2*x + 1) - 407.21
        rho = -166.29*(2*x + 1)**2 + 916.23*(2*x + 1) + 4320.97
        rho_e = -166.67*(2*x + 1)**2 + 910.38*(2*x + 1) + 4101.98
        # Bulk modulus
        K_l = 118.93*10**9
        K_m = 141.25*10**9
        K_h = 144.00*10**9
        K = (-9.785*(2*x + 1)**2 + 51.675*(2*x + 1) + 77.040)*10**9
        # Shear modulus
        G_l = 48.91*10**9
        G_m = 58.33*10**9
        G_h = 59.29*10**9
        G = (-4.230*(2*x + 1)**2 + 22.110*(2*x + 1) + 31.030)*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -99289     # J/mol
        thermodynamics["Enthalpy"] = -96291         # J/mol
        thermodynamics["Entropy"] = 69.429         # J/(mol K)
        thermodynamics["Heat Capacity"] = 50.50    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_pyrrhotite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["S", "Fe"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
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
                minors = ["Ni", "Co", "Cu"]
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
        mineral = "Po"
        #
        # Molar mass
        x = round(rd.uniform(0, 0.17), 4)
        molar_mass_pure = (1-x)*iron[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[12.811, 6.87, 11.885], [117, 3], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 26, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 143.18*10**9
        # Shear modulus
        G = 59.31*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -99289     # J/mol
        thermodynamics["Enthalpy"] = -96291         # J/mol
        thermodynamics["Entropy"] = 69.429         # J/(mol K)
        thermodynamics["Heat Capacity"] = 50.50    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    # def create_cobaltite(self):
    #     # Major elements
    #     sulfur = PeriodicSystem(name="S").get_data()
    #     cobalt = PeriodicSystem(name="Co").get_data()
    #     arsenic = PeriodicSystem(name="As").get_data()
    #     majors_name = ["S", "Co", "As"]
    #     majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Co", cobalt[1], 1, cobalt[2]],
    #                             ["As", arsenic[1], 1, arsenic[2]]], dtype=object)
    #     # Minor elements
    #     traces_data = []
    #     if len(self.traces_list) > 0:
    #         self.impurity = "impure"
    #     if self.impurity == "pure":
    #         var_state = "fixed"
    #     else:
    #         var_state = "variable"
    #         if self.impurity == "random":
    #             self.traces_list = []
    #             minors = ["Cu", "Pb", "Sb", "Fe", "Ni"]
    #             n = rd.randint(1, len(minors))
    #             while len(self.traces_list) < n:
    #                 selection = rd.choice(minors)
    #                 if selection not in self.traces_list and selection not in majors_name:
    #                     self.traces_list.append(selection)
    #                 else:
    #                     continue
    #         traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
    #         x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
    #         for i in range(len(self.traces_list)):
    #             traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
    #         if len(traces_data) > 0:
    #             traces_data = np.array(traces_data, dtype=object)
    #             traces_data = traces_data[traces_data[:, 1].argsort()]
    #     #
    #     data = []
    #     mineral = "Cob"
    #     #
    #     # Molar mass
    #     molar_mass_pure = cobalt[2] + arsenic[2] + sulfur[2]
    #     molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
    #                                   majors=majors_data).calculate_molar_mass()
    #     element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
    #     # Density
    #     dataV = CrystalPhysics([[5.57, 5.582, 5.582], [], "orthorhombic"])
    #     V = dataV.calculate_volume()
    #     dataRho = CrystalPhysics([molar_mass_pure, 4, V])
    #     rho = dataRho.calculate_bulk_density()
    #     rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
    #     # Bulk modulus
    #     K = 139*10**9
    #     # Shear modulus
    #     G = 104*10**9
    #     # Young's modulus
    #     E = (9*K*G)/(3*K + G)
    #     # Poisson's ratio
    #     nu = (3*K - 2*G)/(2*(3*K + G))
    #     # vP/vS
    #     vPvS = ((K + 4/3*G)/G)**0.5
    #     # P-wave velocity
    #     vP = ((K + 4/3*G)/rho)**0.5
    #     # S-wave velocity
    #     vS = (G/rho)**0.5
    #     # Gamma ray
    #     gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
    #     # Photoelectricity
    #     pe = wg(amounts=amounts, elements=element).calculate_pe()
    #     U = pe*rho_e*10**(-3)
    #     # Electrical resistivity
    #     p = None
    #     #
    #     data.append(mineral)
    #     data.append(round(molar_mass, 3))
    #     data.append(round(rho, 2))
    #     data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
    #     data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
    #     data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
    #     data.append(amounts)
    #     #
    #     return data
    #
    def create_carrollite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        cobalt = PeriodicSystem(name="Co").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["S", "Co", "Cu"]
        majors_data = np.array([["S", sulfur[1], 4, sulfur[2]], ["Co", cobalt[1], 2, cobalt[2]],
                                ["Cu", copper[1], 1, copper[2]]], dtype=object)
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
        data = []
        mineral = "Carr"
        #
        # Molar mass
        molar_mass_pure = copper[2] + 2*cobalt[2] + 4*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.477], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 8, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 106.49*10**9
        # Shear modulus
        G = 45.71*10**9
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
        data.append(mineral)
        data.append(round(molar_mass, 3))
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_chalcocite(self):
        #
        name = "Cc"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["S", "Cu"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Cu", copper[1], 2, copper[2]]], dtype=object)
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
                minors = ["Fe"]
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
        molar_mass_pure = 2*copper[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.881, 27.323, 13.491], [116.35], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 96, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 101.98*10**9
        # Shear modulus
        G = 39.12*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -84205    # J/mol
        thermodynamics["Enthalpy"] = -78884        # J/mol
        thermodynamics["Entropy"] = 116.200         # J/(mol K)
        thermodynamics["Heat Capacity"] = 76.32    # J/(mol K)
        #
        if self.data_type == False:
            data = []
            data.append(name)
            data.append(round(molar_mass, 2))
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            results = {}
            results["mineral"] = name
            results["M"] = molar_mass
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            try:
                results["p"] = round(p, 4)
            except:
                results["p"] = p
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            #
            return results
    #
    def create_digenite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["S", "Cu"]
        majors_data = np.array([["S", sulfur[1], 5, sulfur[2]], ["Cu", copper[1], 9, copper[2]]], dtype=object)
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
                minors = ["Fe"]
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
        mineral = "Dg"
        #
        # Molar mass
        molar_mass_pure = 9*copper[2] + 5*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.919, 48], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 3, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 101.13*10**9
        # Shear modulus
        G = 38.28*10**9
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
    def create_tennantite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["S", "Cu", "As"]
        majors_data = np.array([["S", sulfur[1], 13, sulfur[2]], ["Cu", copper[1], 12, copper[2]],
                                ["As", arsenic[1], 4, arsenic[2]]], dtype=object)
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
                minors = ["Ag", "Fe", "Zn", "Pb"]
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
        mineral = "Tn"
        #
        # Molar mass
        molar_mass_pure = 12*copper[2] + 4*arsenic[2] + 13*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.186], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 70.29*10**9
        # Shear modulus
        G = 28.87*10**9
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
        data.append(mineral)
        data.append(round(molar_mass, 3))
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_tetrahedrite(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        antimony = PeriodicSystem(name="Sb").get_data()
        majors_name = ["S", "Cu", "Sb"]
        majors_data = np.array([["S", sulfur[1], 13, sulfur[2]], ["Cu", copper[1], 12, copper[2]],
                                ["Sb", antimony[1], 4, antimony[2]]], dtype=object)
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
                minors = ["Ag", "Fe", "Zn", "Pb"]
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
        mineral = "Td"
        #
        # Molar mass
        molar_mass_pure = 12*copper[2] + 4*antimony[2] + 13*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.33], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 59*10**9
        # Shear modulus
        G = 13*10**9
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
        data.append(mineral)
        data.append(round(molar_mass, 3))
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_tennantite_group(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["S", "Fe", "Cu", "Zn", "As"]
        majors_data = np.array([["S", sulfur[1], 13, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["Cu", copper[1], 10, copper[2]], ["Zn", zinc[1], 1, zinc[2]],
                                ["As", arsenic[1], 4, arsenic[2]]], dtype=object)
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
                minors = ["Ag", "Pb"]
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
        mineral = "Tn"
        #
        # Molar mass
        x = round(rd.uniform(0, 1), 4)
        molar_mass_pure = 10*copper[2] + 2*(x*iron[2] + (1-x)*zinc[2]) + 4*arsenic[2] + 13*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.186], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 70.29*10**9
        # Shear modulus
        G = 28.87*10**9
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
        data.append(mineral)
        data.append([round(molar_mass, 3), x])
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_tetrahedrite_group(self):
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        antimony = PeriodicSystem(name="Sb").get_data()
        majors_name = ["S", "Fe", "Cu", "Zn", "Sb"]
        majors_data = np.array([["S", sulfur[1], 13, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["Cu", copper[1], 10, copper[2]], ["Zn", zinc[1], 1, zinc[2]],
                                ["Sb", antimony[1], 4, antimony[2]]], dtype=object)
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
                minors = ["Ag", "Pb"]
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
        mineral = "Td"
        #
        # Molar mass
        x = round(rd.uniform(0, 1), 4)
        molar_mass_pure = 10*copper[2] + 2*(x*iron[2] + (1-x)*zinc[2]) + 4*antimony[2] + 13*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.33], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 59*10**9
        # Shear modulus
        G = 13*10**9
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
        data.append(mineral)
        data.append([round(molar_mass, 3), x])
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_fahlore(self):
        name = "Fh"
        # Major elements
        majors_name = ["S", "Fe", "Cu", "As", "Sb"]
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
                minors = ["Ag", "Pb", "Zn"]
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
        # Molar mass
        x = round(rd.uniform(0.75, 0.925), 4)
        y = round(rd.uniform(0, 1), 4)
        majors_data = np.array([
            ["S", self.sulfur[1], 13, self.sulfur[2]], ["Fe", self.iron[1], 12*(1-x), self.iron[2]],
            ["Cu", self.copper[1], 12*x, self.copper[2]], ["As", self.arsenic[1], 4*y, self.arsenic[2]],
            ["Sb", self.antimony[1], 4*(1-y), self.antimony[2]]], dtype=object)
        molar_mass_pure = (12*(x*self.copper[2] + (1-x)*self.iron[2]) + 4*(y*self.arsenic[2] + (1-y)*self.antimony[2]) +
                           13*self.sulfur[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

        majors_data_tn = np.array([
            ["S", self.sulfur[1], 13, self.sulfur[2]], ["Fe", self.iron[1], 0, self.iron[2]],
            ["Cu", self.copper[1], 12, self.copper[2]], ["As", self.arsenic[1], 4, self.arsenic[2]],
            ["Sb", self.antimony[1], 0, self.antimony[2]]], dtype=object)
        molar_mass_pure_tn = 12*self.copper[2] + 4*self.arsenic[2] + 13*self.sulfur[2]
        molar_mass_tn, amounts_tn = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_tn, majors=majors_data_tn).calculate_molar_mass()
        element_tn = [PeriodicSystem(name=amounts_tn[i][0]).get_data() for i in range(len(amounts_tn))]

        majors_data_td = np.array([
            ["S", self.sulfur[1], 13, self.sulfur[2]], ["Fe", self.iron[1], 0, self.iron[2]],
            ["Cu", self.copper[1], 12, self.copper[2]], ["As", self.arsenic[1], 0, self.arsenic[2]],
            ["Sb", self.antimony[1], 4, self.antimony[2]]], dtype=object)
        molar_mass_pure_td = 12*self.copper[2] + 4*self.antimony[2] + 13*self.sulfur[2]
        molar_mass_td, amounts_td = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_td, majors=majors_data_td).calculate_molar_mass()
        element_td = [PeriodicSystem(name=amounts_td[i][0]).get_data() for i in range(len(amounts_td))]
        # Density
        dataV_tn = CrystalPhysics([[10.186], [], "cubic"])
        V_tn = dataV_tn.calculate_volume()
        Z_tn = 2
        V_m_tn = MineralChemistry().calculate_molar_volume(volume_cell=V_tn, z=Z_tn)
        dataRho_tn = CrystalPhysics([molar_mass_tn, Z_tn, V_tn])
        rho_tn = dataRho_tn.calculate_bulk_density()
        rho_e_tn = wg(amounts=amounts_tn, elements=element_tn, rho_b=rho_tn).calculate_electron_density()

        dataV_td = CrystalPhysics([[10.33], [], "cubic"])
        V_td = dataV_td.calculate_volume()
        Z_td = 2
        V_m_td = MineralChemistry().calculate_molar_volume(volume_cell=V_td, z=Z_td)
        dataRho_td = CrystalPhysics([molar_mass_td, Z_td, V_td])
        rho_td = dataRho_td.calculate_bulk_density()
        rho_e_td = wg(amounts=amounts_td, elements=element_td, rho_b=rho_td).calculate_electron_density()

        V_m = y*V_m_tn + (1 - y)*V_m_td
        rho = y*rho_tn + (1 - y)*rho_td
        rho_e = y*rho_e_tn + (1 - y)*rho_e_td
        # Bulk modulus
        K_tn = 70.29*10**9
        K_td = 59*10**9
        K = y*K_tn + (1-y)*K_td
        # Shear modulus
        G_tn = 28.87*10**9
        G_td = 13*10**9
        G = y*G_tn + (1-y)*G_td
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
        # Results
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

        return results

    def create_covellite(self):
        #
        name = "Cv"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["S", "Cu"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Cu", copper[1], 1, copper[2]]], dtype=object)
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
                minors = ["Fe", "Se", "Ag", "Pb"]
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
        molar_mass_pure = copper[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.792, 16.344], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass_pure, 6, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 74*10**9
        # Shear modulus
        G = 24*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -48636     # J/mol
        thermodynamics["Enthalpy"] = -47982         # J/mol
        thermodynamics["Entropy"] = 67.400          # J/(mol K)
        thermodynamics["Heat Capacity"] = 47.54     # J/(mol K)
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
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_millerite(self):
        #
        name = "Mlr"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["S", "Ni"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["Ni", nickel[1], 1, nickel[2]]], dtype=object)
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
                minors = ["Fe", "Co", "Cu"]
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
        molar_mass_pure = nickel[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.62, 3.149], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 9, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 129.85*10**9
        # Shear modulus
        G = 53.78*10**9
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
    def create_pentlandite(self):
        #
        name = "Pn"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["S", "Fe", "Ni"]
        x = round(rd.uniform(0.4, 0.6), 4)
        majors_data = np.array([["S", sulfur[1], 8, sulfur[2]], ["Fe", iron[1], x*9, iron[2]],
                                ["Ni", nickel[1], (1-x)*9, nickel[2]]], dtype=object)
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
                minors = ["Co", "Ag", "Cu"]
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
        molar_mass_pure = 9*(x*iron[2] + (1-x)*nickel[2]) + 8*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.04], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 113.02*10**9
        # Shear modulus
        G = 48.54*10**9
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
    def create_realgar(self):
        #
        name = "Rl"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["S", "As"]
        majors_data = np.array([["S", sulfur[1], 1, sulfur[2]], ["As", arsenic[1], 1, arsenic[2]]], dtype=object)
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
        molar_mass_pure = arsenic[2] + sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.29, 13.53, 6.57], [106.883], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 16, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 25.6*10**9
        # Shear modulus
        G = 13.64*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -30519     # J/mol
        thermodynamics["Enthalpy"] = -32417         # J/mol
        thermodynamics["Entropy"] = 61.379         # J/(mol K)
        thermodynamics["Heat Capacity"] = 60.70    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_orpiment(self):
        #
        name = "Orp"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["S", "As"]
        majors_data = np.array([["S", sulfur[1], 3, sulfur[2]], ["As", arsenic[1], 2, arsenic[2]]], dtype=object)
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
                minors = ["Hg", "Ge", "Sb"]
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
        molar_mass_pure = 2*arsenic[2] + 3*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.49, 9.59, 4.25], [90.45], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 27.89*10**9
        # Shear modulus
        G = 14.22*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -85981     # J/mol
        thermodynamics["Enthalpy"] = -87158         # J/mol
        thermodynamics["Entropy"] = 163.594         # J/(mol K)
        thermodynamics["Heat Capacity"] = 115.19    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_marcasite(self):
        #
        name = "Mr"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["S", "Fe"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
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
                minors = ["Cu", "As"]
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
        molar_mass_pure = iron[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.445, 5.425, 3.388], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 150*10**9
        # Shear modulus
        G = 125*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = -155411     # J/mol
        thermodynamics["Enthalpy"] = -166600         # J/mol
        thermodynamics["Entropy"] = 53.900          # J/(mol K)
        thermodynamics["Heat Capacity"] = 62.43     # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_chalcopyrite_group(self):
        #
        name = "Ccp"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        gallium = PeriodicSystem(name="Ga").get_data()
        silver = PeriodicSystem(name="Ag").get_data()
        indium = PeriodicSystem(name="In").get_data()
        #
        majors_name = ["S", "Fe", "Cu", "Ga", "Ag", "In"]
        #
        a = rd.uniform(0.0, 1.0)
        a1 = round(a, 4)
        a2 = round(1 - a, 4)
        b = rd.uniform(0.0, 1.0)
        c = rd.uniform(0.0, 1.0-b)
        b1 = round(b, 4)
        b2 = round(c, 4)
        b3 = round(1 - b1 - b2, 4)
        c = rd.uniform(0.0, 1.0)
        c1 = round(c, 4)
        c2 = round(1 - c, 4)
        #
        majors_data = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], b1, iron[2]], ["Cu", copper[1], a1, copper[2]],
             ["Ga", gallium[1], b3, gallium[2]], ["Ag", silver[1], a2, silver[2]],
             ["In", indium[1], b2, indium[2]]], dtype=object)
        majors_data = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], a1*b1 + a2*c1, iron[2]], ["Cu", copper[1], a1, copper[2]],
             ["Ga", gallium[1], b3, gallium[2]], ["Ag", silver[1], a2, silver[2]],
             ["In", indium[1], a1*b2 + a2*c2, indium[2]]], dtype=object)
        majors_data_Ccp = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 1, iron[2]], ["Cu", copper[1], 1, copper[2]],
             ["Ga", gallium[1], 0, gallium[2]], ["Ag", silver[1], 0, silver[2]],
             ["In", indium[1], 0, indium[2]]], dtype=object)
        majors_data_Rqt = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 0, iron[2]], ["Cu", copper[1], 1, copper[2]],
             ["Ga", gallium[1], 0, gallium[2]], ["Ag", silver[1], 0, silver[2]],
             ["In", indium[1], 1, indium[2]]], dtype=object)
        majors_data_Glt = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 0, iron[2]], ["Cu", copper[1], 1, copper[2]],
             ["Ga", gallium[1], 1, gallium[2]], ["Ag", silver[1], 0, silver[2]],
             ["In", indium[1], 0, indium[2]]], dtype=object)
        majors_data_Lnt = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 1, iron[2]], ["Cu", copper[1], 0, copper[2]],
             ["Ga", gallium[1], 0, gallium[2]], ["Ag", silver[1], 1, silver[2]],
             ["In", indium[1], 0, indium[2]]], dtype=object)
        majors_data_Lft = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 0, iron[2]], ["Cu", copper[1], 0, copper[2]],
             ["Ga", gallium[1], 0, gallium[2]], ["Ag", silver[1], 1, silver[2]],
             ["In", indium[1], 1, indium[2]]], dtype=object)
        #
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
                minors = ["Ni"]
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
        molar_mass_pure = a1*copper[2] + a2*silver[2] + (b1 + c1)*iron[2] + (b2 + c2)*indium[2] + b3*gallium[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        molar_mass_pure_Ccp = copper[2] + iron[2] + 2*sulfur[2]
        molar_mass_Ccp, amounts_Ccp = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Ccp, majors=majors_data_Ccp).calculate_molar_mass()
        element_Ccp = [PeriodicSystem(name=amounts_Ccp[i][0]).get_data() for i in range(len(amounts_Ccp))]
        #
        molar_mass_pure_Rqt = copper[2] + indium[2] + 2*sulfur[2]
        molar_mass_Rqt, amounts_Rqt = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Rqt, majors=majors_data_Rqt).calculate_molar_mass()
        element_Rqt = [PeriodicSystem(name=amounts_Rqt[i][0]).get_data() for i in range(len(amounts_Rqt))]
        #
        molar_mass_pure_Glt = copper[2] + gallium[2] + 2*sulfur[2]
        molar_mass_Glt, amounts_Glt = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Glt, majors=majors_data_Glt).calculate_molar_mass()
        element_Glt = [PeriodicSystem(name=amounts_Glt[i][0]).get_data() for i in range(len(amounts_Glt))]
        #
        molar_mass_pure_Lnt = silver[2] + iron[2] + 2*sulfur[2]
        molar_mass_Lnt, amounts_Lnt = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Lnt, majors=majors_data_Lnt).calculate_molar_mass()
        element_Lnt = [PeriodicSystem(name=amounts_Lnt[i][0]).get_data() for i in range(len(amounts_Lnt))]
        #
        molar_mass_pure_Lft = silver[2] + indium[2] + 2*sulfur[2]
        molar_mass_Lft, amounts_Lft = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Lft, majors=majors_data_Lft).calculate_molar_mass()
        element_Lft = [PeriodicSystem(name=amounts_Lft[i][0]).get_data() for i in range(len(amounts_Lft))]
        #
        # Density
        dataV_Ccp = CrystalPhysics([[5.28, 10.41], [], "tetragonal"])
        V_Ccp = dataV_Ccp.calculate_volume()
        dataRho_Ccp = CrystalPhysics([molar_mass_Ccp, 4, V_Ccp])
        rho_Ccp = dataRho_Ccp.calculate_bulk_density()
        rho_e_Ccp = wg(amounts=amounts_Ccp, elements=element_Ccp, rho_b=rho_Ccp).calculate_electron_density()
        #
        dataV_Glt = CrystalPhysics([[5.35, 10.48], [], "tetragonal"])
        V_Glt = dataV_Glt.calculate_volume()
        dataRho_Glt = CrystalPhysics([molar_mass_Glt, 4, V_Glt])
        rho_Glt = dataRho_Glt.calculate_bulk_density()
        rho_e_Glt = wg(amounts=amounts_Glt, elements=element_Glt, rho_b=rho_Glt).calculate_electron_density()
        #
        dataV_Rqt = CrystalPhysics([[5.51, 11.05], [], "tetragonal"])
        V_Rqt = dataV_Rqt.calculate_volume()
        dataRho_Rqt = CrystalPhysics([molar_mass_Rqt, 4, V_Rqt])
        rho_Rqt = dataRho_Rqt.calculate_bulk_density()
        rho_e_Rqt = wg(amounts=amounts_Rqt, elements=element_Rqt, rho_b=rho_Rqt).calculate_electron_density()
        #
        dataV_Lft = CrystalPhysics([[5.88, 11.21], [], "tetragonal"])
        V_Lft = dataV_Lft.calculate_volume()
        dataRho_Lft = CrystalPhysics([molar_mass_Lft, 4, V_Lft])
        rho_Lft = dataRho_Lft.calculate_bulk_density()
        rho_e_Lft = wg(amounts=amounts_Lft, elements=element_Lft, rho_b=rho_Lft).calculate_electron_density()
        #
        dataV_Lnt = CrystalPhysics([[5.4371, 10.8479], [], "tetragonal"])
        V_Lnt = dataV_Lnt.calculate_volume()
        dataRho_Lnt = CrystalPhysics([molar_mass_Lnt, 4, V_Lnt])
        rho_Lnt = dataRho_Lnt.calculate_bulk_density()
        rho_e_Lnt = wg(amounts=amounts_Lnt, elements=element_Lnt, rho_b=rho_Lnt).calculate_electron_density()
        #
        V = (a*b)*V_Ccp + ((1-a)*b)*V_Lnt + (a*c)*V_Rqt + (a*(1-b-c))*V_Glt + ((1-a)*c)*V_Lft
        rho = a1*(b1*rho_Ccp + b2*rho_Rqt + b3*rho_Glt) + a2*(c1*rho_Lnt + c2*rho_Lft)
        rho_e = (a*b)*rho_e_Ccp + ((1-a)*b)*rho_e_Lnt + (a*c)*rho_e_Rqt + (a*(1-b-c))*rho_e_Glt + ((1-a)*c)*rho_e_Lft
        #
        # Bulk modulus
        K_Ccp = 56*10**9
        K_Lnt = 58*10**9
        K_Rqt = 64*10**9
        K_Glt = 75*10**9
        K_Lft = 52*10**9
        K = a1*(b1*K_Ccp + b2*K_Rqt + b3*K_Glt) + a2*(c1*K_Lnt + c2*K_Lft)
        # Shear modulus
        G_Ccp = 19*10**9
        G_Lnt = 15*10**9
        G_Rqt = 25*10**9
        G_Glt = 36*10**9
        G_Lft = 16*10**9
        G = a1*(b1*G_Ccp + b2*G_Rqt + b3*G_Glt) + a2*(c1*G_Lnt + c2*G_Lft)
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None       # J/mol
        thermodynamics["Enthalpy"] = None           # J/mol
        thermodynamics["Entropy"] = None            # J/(mol K)
        thermodynamics["Heat Capacity"] = None      # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_gallite(self):
        #
        name = "Glt"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        gallium = PeriodicSystem(name="Ga").get_data()
        majors_name = ["S", "Cu", "Ga"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Cu", copper[1], 1, copper[2]],
                                ["Ga", gallium[1], 1, gallium[2]]], dtype=object)
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
                minors = ["Pb", "Zn", "Fe", "Ge"]
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
        molar_mass_pure = copper[2] + gallium[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.35, 10.48], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 75*10**9
        # Shear modulus
        G = 36*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_roquesite(self):
        #
        name = "Rqt"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        indium = PeriodicSystem(name="In").get_data()
        majors_name = ["S", "Cu", "In"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Cu", copper[1], 1, copper[2]],
                                ["In", indium[1], 1, indium[2]]], dtype=object)
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
                minors = ["Fe"]
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
        molar_mass_pure = copper[2] + indium[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.51, 11.05], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 64*10**9
        # Shear modulus
        G = 25*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_lenaite(self):
        #
        name = "Lnt"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        silver = PeriodicSystem(name="Ag").get_data()
        majors_name = ["S", "Fe", "Ag"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], 1, iron[2]],
                                ["Ag", silver[1], 1, silver[2]]], dtype=object)
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
        molar_mass_pure = silver[2] + iron[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.4371, 10.8479], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 58*10**9
        # Shear modulus
        G = 15*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_laforetite(self):
        #
        name = "Lft"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        silver = PeriodicSystem(name="Ag").get_data()
        indium = PeriodicSystem(name="In").get_data()
        majors_name = ["S", "Ag", "In"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Ag", silver[1], 1, silver[2]],
                                ["In", indium[1], 1, indium[2]]], dtype=object)
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
                minors = ["Cu", "Fe", "Se"]
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
        molar_mass_pure = silver[2] + indium[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.88, 11.21], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 52*10**9
        # Shear modulus
        G = 16*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
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
            results["M"] = molar_mass
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V, 4)
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
            results["thermodynamics"] = thermodynamics
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_vaesite(self):   # Ni S2
        #
        name = "Vst"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["S", "Ni"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Ni", nickel[1], 1, nickel[2]]], dtype=object)
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
                minors = ["Co", "Fe"]
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
        molar_mass_pure = nickel[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.6793], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 101.63*10**9
        # Shear modulus
        G = 42.05*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
        #
        # Output
        results = {}
        results["mineral"] = name
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
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
        results["thermodynamics"] = thermodynamics
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_cattierite(self):   # Co S2
        #
        name = "Cat"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        cobalt = PeriodicSystem(name="Co").get_data()
        majors_name = ["S", "Co"]
        majors_data = np.array([["S", sulfur[1], 2, sulfur[2]], ["Co", cobalt[1], 1, cobalt[2]]], dtype=object)
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
                minors = ["Fe", "Ni"]
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
        molar_mass_pure = cobalt[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.535], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 4, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 124*10**9
        # Shear modulus
        G = 73*10**9
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
        #
        # Output
        results = {}
        results["mineral"] = name
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
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
        results["thermodynamics"] = thermodynamics
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_pyrite_group(self):   # (Fe,Ni,Co) S2
        #
        name = "Py"
        #
        # Major elements
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        cobalt = PeriodicSystem(name="Co").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["S", "Co"]
        #
        a = rd.uniform(0.5, 1.0)
        b = rd.uniform(0.0, (1 - a))
        c = 1 - a - b
        #
        majors_data = np.array(
            [["S", sulfur[1], 2, sulfur[2]], ["Fe", iron[1], a, iron[2]], ["Co", cobalt[1], c, cobalt[2]],
             ["Ni", nickel[1], b, nickel[2]]], dtype=object)
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
                minors = ["Fe", "Ni", "Co"]
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
        molar_mass_pure = a*iron[2] + b*nickel[2] + c*cobalt[2] + 2*sulfur[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[5.417], [], "cubic"])
        V_Fe = dataV_Fe.calculate_volume()
        dataRho_Fe = CrystalPhysics([molar_mass, 4, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Ni = CrystalPhysics([[5.535], [], "cubic"])
        V_Ni = dataV_Ni.calculate_volume()
        dataRho_Ni = CrystalPhysics([molar_mass, 4, V_Ni])
        rho_Ni = dataRho_Ni.calculate_bulk_density()
        rho_e_Ni = wg(amounts=amounts, elements=element, rho_b=rho_Ni).calculate_electron_density()
        #
        dataV_Co = CrystalPhysics([[5.535], [], "cubic"])
        V_Co = dataV_Co.calculate_volume()
        dataRho_Co = CrystalPhysics([molar_mass, 4, V_Co])
        rho_Co = dataRho_Co.calculate_bulk_density()
        rho_e_Co = wg(amounts=amounts, elements=element, rho_b=rho_Co).calculate_electron_density()
        #
        V = a*V_Fe + b*V_Ni + c*V_Co
        rho = a*rho_Fe + b*rho_Ni + c*rho_Co
        rho_e = a*rho_e_Fe + b*rho_e_Ni + c*rho_e_Co
        #
        # Bulk modulus
        K_Fe = 146*10**9
        K_Ni = 124*10**9
        K_Co = 124*10**9
        K = a*K_Fe + b*K_Ni + c*K_Co
        # Shear modulus
        G_Fe = 135*10**9
        G_Ni = 73*10**9
        G_Co = 73*10**9
        G = a*G_Fe + b*G_Ni + c*G_Co
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
        #
        # Output
        results = {}
        results["mineral"] = name
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
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
        results["thermodynamics"] = thermodynamics
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_cobaltite(self):  # Co As S
        ## General Information
        name = "Cbt"
        elements_list = ["S", "Co", "As"]
        #
        sulfur = PeriodicSystem(name="S").get_data()
        cobalt = PeriodicSystem(name="Co").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        #
        molar_mass_ideal = cobalt[2] + arsenic[2] + sulfur[2]
        amounts_elements = {
            "S": round(sulfur[2]/molar_mass_ideal, 6),
            "Co": round(cobalt[2]/molar_mass_ideal, 6),
            "As": round(arsenic[2]/molar_mass_ideal, 6)}
        #
        ## Trace elements
        composition_sulfide = {}
        for element in elements_list:
            composition_sulfide[element] = int(amounts_elements[element]*10**6)
        #
        element_traces = {
            "3+": ["Sb", "Fe"],
            "2+": ["Cu", "Pb", "Ni"],
            "All": ["Cu", "Pb", "Sb", "Fe", "Ni"]}
        #
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
            #
            for trace_element, value in self.traces_list.items():
                if trace_element in ["Sb", "Fe"]:
                    elements_list.append(trace_element)
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0])
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_sulfide[trace_element] = amount_ppm
                    composition_sulfide["As"] -= amount_ppm
                    #
                elif trace_element in ["Cu", "Pb", "Ni"]:
                    elements_list.append(trace_element)
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0])
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_sulfide[trace_element] = amount_ppm
                    composition_sulfide["Co"] -= amount_ppm
            #
        else:
            self.impurity == "pure"
            var_state = "fixed"
        #
        compositon_data = TraceElements(
            tracer=self.traces_list).calculate_composition_sulfides(
            var_elements=elements_list, var_composition=composition_sulfide, var_mineral="Cobaltite")
        #
        ## Molar mass
        molar_mass_pure = molar_mass_ideal
        molar_mass = 0
        amounts = []
        #
        for element in compositon_data:
            chem_data = PeriodicSystem(name=element).get_data()
            molar_mass += compositon_data[element]["x"] * chem_data[2]
            amounts.append([chem_data[0], chem_data[1], compositon_data[element]["w"]])
        #
        magic_factor = molar_mass / molar_mass_pure
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        ## Density
        dataV = CrystalPhysics([[5.57, 5.582, 5.582], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)*magic_factor
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()*magic_factor
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        ## Bulk modulus
        K = 139*10**9*magic_factor
        ## Shear modulus
        G = 104*10**9*magic_factor
        ## Young's modulus
        E = (9*K*G)/(3*K + G)
        ## Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        ## vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        ## P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        ## S-wave velocity
        vS = (G/rho)**0.5
        ## Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        ## Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        ## Electrical resistivity
        p = None
        #
        ## Data Export
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_marmatite(self):   # (Zn,Fe) S
        ## General Information
        name = "Sp"
        elements_list = ["S", "Fe", "Zn"]
        #
        sulfur = PeriodicSystem(name="S").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        #
        majors_data_fe = np.array(
            [["S", sulfur[1], 1, sulfur[2]], ["Fe", iron[1], 1, iron[2]], ["Zn", zinc[1], 0, zinc[2]]], dtype=object)
        majors_data_zn = np.array(
            [["S", sulfur[1], 1, sulfur[2]], ["Fe", iron[1], 0, iron[2]], ["Zn", zinc[1], 1, zinc[2]]], dtype=object)
        #
        x = round(rd.uniform(0.95, 1.0), 4)
        #
        molar_mass_x = x*zinc[2] + (1 - x)*iron[2] + sulfur[2]
        amounts_elements = {
            "S": round(sulfur[2]/molar_mass_x, 6),
            "Fe": round((1 - x)*iron[2]/molar_mass_x, 6),
            "Zn": round(x*zinc[2]/molar_mass_x, 6)}
        #
        ## Trace elements
        composition_sulfide = {}
        for element in elements_list:
            composition_sulfide[element] = int(amounts_elements[element]*10**6)
        #
        element_traces = {
            "4+": ["Ge"],
            "3+": ["In", "Ga", "Sb"],
            "2+": ["Mn", "Cd", "Hg", "Sn", "Pb", "Co"],
            "1+": ["Tl", "Ag"],
            "All": ["Mn", "Cd", "Hg", "In", "Tl", "Ga", "Ge", "Sb", "Sn", "Pb", "Ag", "Co"]}
        #
        if len(self.traces_list) > 0:
            var_state = "variable"
            #
            for trace_element, value in self.traces_list.items():
                if trace_element in ["Mn", "Cd", "Hg", "In", "Tl", "Ga", "Ge", "Sb", "Sn", "Pb", "Ag", "Co"]:
                    elements_list.append(trace_element)
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0])
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    if composition_sulfide["Zn"] >= composition_sulfide["Fe"]:
                        composition_sulfide[trace_element] = amount_ppm
                        composition_sulfide["Zn"] -= amount_ppm
                    else:
                        composition_sulfide[trace_element] = amount_ppm
                        composition_sulfide["Fe"] -= amount_ppm
            #
        else:
            var_state = "fixed"
        #
        compositon_data = TraceElements(
            tracer=self.traces_list).calculate_composition_sulfides(
            var_elements=elements_list, var_composition=composition_sulfide, var_mineral="Marmatite", var_x=x)
        #
        ## Molar mass
        molar_mass_pure = molar_mass_x
        molar_mass = 0
        amounts = []
        #
        for element in compositon_data:
            chem_data = PeriodicSystem(name=element).get_data()
            molar_mass += compositon_data[element]["x"] * chem_data[2]
            amounts.append([chem_data[0], chem_data[1], compositon_data[element]["w"]])
        #
        magic_factor = molar_mass/molar_mass_pure
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        molar_mass_pure_fe = iron[2] + sulfur[2]
        molar_mass_fe, amounts_fe = MineralChemistry(
            w_traces=[], molar_mass_pure=molar_mass_pure_fe, majors=majors_data_fe).calculate_molar_mass()
        element_fe = [PeriodicSystem(name=amounts_fe[i][0]).get_data() for i in range(len(amounts_fe))]
        #
        molar_mass_pure_zn = zinc[2] + sulfur[2]
        molar_mass_zn, amounts_zn = MineralChemistry(
            w_traces=[], molar_mass_pure=molar_mass_pure_zn, majors=majors_data_zn).calculate_molar_mass()
        element_zn = [PeriodicSystem(name=amounts_zn[i][0]).get_data() for i in range(len(amounts_zn))]
        #
        # Density
        dataV_Fe = CrystalPhysics([[3.60, 5.45], [], "tetragonal"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 4
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass_fe, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts_fe, elements=element_fe, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Zn = CrystalPhysics([[5.406], [], "cubic"])
        V_Zn = dataV_Zn.calculate_volume()
        Z_Zn = 4
        V_m_Zn = MineralChemistry().calculate_molar_volume(volume_cell=V_Zn, z=Z_Zn)
        dataRho_Zn = CrystalPhysics([molar_mass_zn, Z_Zn, V_Zn])
        rho_Zn = dataRho_Zn.calculate_bulk_density()
        rho_e_Zn = wg(amounts=amounts_zn, elements=element_zn, rho_b=rho_Zn).calculate_electron_density()
        #
        V_m = (x*V_m_Zn + (1 - x)*V_m_Fe)*magic_factor
        rho = (x*rho_Zn + (1 - x)*rho_Fe)*magic_factor
        rho_e = (x*rho_e_Zn + (1 - x)*rho_e_Fe)*magic_factor
        #
        # Bulk modulus
        K_Fe = 30*10**9
        K_Zn = 68*10**9
        K = (x*K_Zn + (1 - x)*K_Fe)*magic_factor
        # Shear modulus
        G_Fe = 19*10**9
        G_Zn = 33*10**9
        G = (x*G_Zn + (1 - x)*G_Fe)*magic_factor
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
        # Thermodynamics
        thermodynamics = {}
        thermodynamics["Gibbs Energy"] = None     # J/mol
        thermodynamics["Enthalpy"] = None         # J/mol
        thermodynamics["Entropy"] = None          # J/(mol K)
        thermodynamics["Heat Capacity"] = None    # J/(mol K)
        #
        # Output
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
        results["thermodynamics"] = thermodynamics
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
