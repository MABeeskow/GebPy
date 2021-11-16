#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		pyllosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		15.11.2021

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from scipy import stats
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry

# Halogenes
class Pyllosilicates:
    """ Class that generates geophysical and geochemical data of phyllosilicate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", dict=False):
        self.traces_list = traces_list
        self.impurity = impurity
        self.dict = dict
    #
    def create_illite(self): # (K,H3O) (Al,Mg,Fe)2 (Si,Al)4 O10 [(OH)2,(H2O)]
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Mg", "Al", "Si", "K", "Fe"]
        #
        a = round(rd.uniform(0.6, 0.8), 4)
        b = round(rd.uniform(0.75, 1.0), 4)
        b2 = round(rd.uniform(0.0, float(1-b)), 4)
        c = round(rd.uniform(0.975, 1.0), 4)
        d = round(rd.uniform(0.75, 1.0), 4)
        #
        majors_data = np.array([["H", hydrogen[1], 3*(1-a)+4*d, hydrogen[2]], ["O", oxygen[1], (1-a)+10+3*d, oxygen[2]],
                                ["Mg", magnesium[1], 2*b2, magnesium[2]], ["Al", aluminium[1], 2*b+4*(1-c), aluminium[2]],
                                ["Si", silicon[1], 4*c, silicon[2]], ["K", potassium[1], a, potassium[2]],
                                ["Fe", iron[1], 2*(1-b-b2), iron[2]]], dtype=object)
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
        mineral = "Ilt"
        #
        # Molar mass
        molar_mass_pure = a*potassium[2] + (1-a)*(3*hydrogen[2] + oxygen[2]) \
                          + 2*(b*aluminium[2] + b2*magnesium[2] + (1-b-b2)*iron[2]) \
                          + 4*(c*silicon[2] + (1-c)*aluminium[2]) + 10*oxygen[2] \
                          + d*(2*(oxygen[2]+hydrogen[2]) + (2*hydrogen[2]+oxygen[2]))
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.18, 8.98, 10.32], [101.83], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = (35.72 + (62.21-35.72)/(2.706-2.546)*(rho/1000-2.546))*10**9
        # Shear modulus
        G = (17.80 + (25.70-17.80)/(2.706-2.546)*(rho/1000-2.546))*10**9
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
        if self.dict == False:
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
            results["M"] = round(molar_mass, 3)
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
    def create_kaolinite(self): # Al2(OH)4Si2O5
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 4, hydrogen[2]], ["O", oxygen[1], 9, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 2, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Mg", "Na", "K", "Ti", "Ca"]
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
        mineral = "Kln"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + 4*(oxygen[2]+hydrogen[2]) + 2*silicon[2] + 5*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.13, 8.89, 7.25], [90.0, 104.5, 89.8], "triclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 2, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 122.54*10**9
        # Shear modulus
        G = 66.63*10**9
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
        if self.dict == False:
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
            results["M"] = round(molar_mass, 3)
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
    def create_montmorillonite(self): # (Na,Ca)0.3(Al,Mg)2Si4O10(OH2)(H2O)10
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "O", "Na", "Mg", "Al", "Si", "Ca"]
        #
        x = round(rd.uniform(0.6, 0.7), 2)
        y = round(rd.uniform(0.9, 1), 2)
        n = rd.randint(10, 12)
        #
        majors_data = np.array([["H", hydrogen[1], 2+2*n, hydrogen[2]], ["O", oxygen[1], 10+2+n, oxygen[2]],
                                ["Na", sodium[1], 0.3*x, sodium[2]], ["Mg", magnesium[1], 2*(1-y), magnesium[2]],
                                ["Al", aluminium[1], 2*y, aluminium[2]], ["Si", silicon[1], 4, silicon[2]],
                                ["Ca", calcium[1], 0.3*(1-x), calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "K"]
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
        mineral = "Mnt"
        #
        # Molar mass
        molar_mass_pure = 0.3*(x*sodium[2]+(1-x)*calcium[2]) + 2*(y*aluminium[2]+(1-y)*magnesium[2]) + 4*silicon[2] \
                          + 10*oxygen[2] + 2*(hydrogen[2]+oxygen[2]) + n*(2*hydrogen[2]+oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.17, 8.94, 9.95], [99.9], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 1, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2738, 2788, 3250, 2504, 3182]
        y_K = [37.30, 35.31, 49.46, 29.71, 66.59]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        y_G = [17.00, 20.19, 24.70, 16.30, 27.00]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
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
        if self.dict == False:
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
            results["M"] = round(molar_mass, 3)
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