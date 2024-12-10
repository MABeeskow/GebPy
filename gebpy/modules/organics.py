#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		organics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		21.09.2022

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry

# Halogenes
class Organics:
    """ Class that generates geophysical and geochemical data of halogene minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1):
        if self.mineral == "Organic Matter":
            if number > 1:
                data = [self.create_organic_matter() for n in range(number)]
            else:
                data = self.create_organic_matter()
        else:
            data = "Nothing found!"
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Organic Matter":
                data_mineral = self.create_organic_matter()
            elif self.mineral == "Lignin":
                data_mineral = self.create_lignin()
            elif self.mineral == "Lipid":
                data_mineral = self.create_lipid()
            elif self.mineral == "Carbohydrate":
                data_mineral = self.create_carbohydrates()
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

        # if self.mineral == "Organic Matter":
        #     dataset = [self.create_organics_matter() for n in range(number)]
        #
        # return dataset
    #
    def create_organic_matter(self):
        # Major elements
        carbohydrates = Organics(data_type=True).create_carbohydrates()
        lignin = Organics(data_type=True).create_lignin()
        lipid = Organics(data_type=True).create_lipid()
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
        data = []
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
        #
        majors_data = np.array([["H", hydrogen[1], w_ch*carbohydrates["chemistry"]["H"] + w_lg*lignin["chemistry"]["H"] + w_lp*lipid["chemistry"]["H"], hydrogen[2]],
                                ["C", carbon[1], w_ch*carbohydrates["chemistry"]["C"] + w_lg*lignin["chemistry"]["C"] + w_lp*lipid["chemistry"]["C"], carbon[2]],
                                ["N", nitrogen[1], w_lg*lignin["chemistry"]["N"], nitrogen[2]],
                                ["O", oxygen[1], w_ch*carbohydrates["chemistry"]["O"] + w_lg*lignin["chemistry"]["O"] + w_lp*lipid["chemistry"]["O"], oxygen[2]],
                                ["S", sulfur[1], w_lg*lignin["chemistry"]["S"], sulfur[2]]], dtype=object)
        #
        molar_mass_pure = w_ch*carbohydrates["M"] + w_lg*lignin["M"] + w_lp*lipid["M"]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        amounts2 = []
        w_H = round((w_ch*carbohydrates["chemistry"]["H"] + w_lg*lignin["chemistry"]["H"] + w_lp*lipid["chemistry"]["H"])*hydrogen[2]/molar_mass_pure, 4)
        w_C = round((w_ch*carbohydrates["chemistry"]["C"] + w_lg*lignin["chemistry"]["C"] + w_lp*lipid["chemistry"]["C"])*carbon[2]/molar_mass_pure, 4)
        w_N = round((w_lg*lignin["chemistry"]["N"])*nitrogen[2]/molar_mass_pure, 4)
        w_S = round((w_lg*lignin["chemistry"]["S"])*sulfur[2]/molar_mass_pure, 4)
        w_O = round(1 - w_H - w_C - w_N - w_S, 4)
        data_H = ["H", 1, w_H]
        data_C = ["C", 6, w_C]
        data_N = ["N", 7, w_N]
        data_O = ["O", 8, w_O]
        data_S = ["S", 16, w_S]
        amounts2.extend([data_H, data_C, data_N, data_O, data_S])
        amounts = amounts2

        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = w_ch*carbohydrates["rho"] + w_lg*lignin["rho"] + w_lp*lipid["rho"]
        V_m = molar_mass/rho
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = w_ch*carbohydrates["K"] + w_lg*lignin["K"] + w_lp*lipid["K"]*10**9
        # Shear modulus
        G = w_ch*carbohydrates["G"] + w_lg*lignin["G"] + w_lp*lipid["G"]*10**9
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
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
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
    def create_carbohydrates(self):
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        majors_name = ["H", "C", "O"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        mineral = "carbohydrates"
        #
        # Molar mass
        molar_mass_pure = 0.06*hydrogen[2] + 0.44*carbon[2] + 0.50*oxygen[2]
        #
        majors_data = np.array([["H", hydrogen[1], 0.06, hydrogen[2]], ["C", carbon[1], 0.44, carbon[2]],
                                ["O", oxygen[1], 0.50, oxygen[2]]], dtype=object)
        #
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = 1586
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 10*10**9
        # Shear modulus
        G = 5*10**9
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
    def create_lignin(self):
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
        mineral = "lignin"
        #
        # Molar mass
        molar_mass_pure = 0.06*hydrogen[2] + 0.626*carbon[2] + 0.003*nitrogen[2] + 0.31*oxygen[2] + 0.001*sulfur[2]
        #
        majors_data = np.array([["H", hydrogen[1], 0.06, hydrogen[2]], ["C", carbon[1], 0.626, carbon[2]],
                                ["N", nitrogen[1], 0.003, nitrogen[2]], ["O", oxygen[1], 0.31, oxygen[2]],
                                ["S", sulfur[1], 0.001, sulfur[2]]], dtype=object)
        #
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = 680
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 4*10**9
        # Shear modulus
        G = 2*10**9
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
    def create_lipid(self):
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        majors_name = ["H", "C", "O"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        mineral = "lipid"
        #
        # Molar mass
        molar_mass_pure = 0.10*hydrogen[2] + 0.80*carbon[2] + 0.10*oxygen[2]
        #
        majors_data = np.array([["H", hydrogen[1], 0.10, hydrogen[2]], ["C", carbon[1], 0.80, carbon[2]], ["O", oxygen[1], 0.10, oxygen[2]]], dtype=object)
        #
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = 850
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 4*10**9
        # Shear modulus
        G = 2*10**9
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