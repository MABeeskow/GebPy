#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		phospides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		20.09.2022

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem
from modules.minerals import CrystalPhysics
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry, TraceElements

# OXIDES
class Phospides:
    """ Class that generates geophysical and geochemical data of phospide minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1):
        if self.mineral == "Allabogdanite":
            data = self.create_allabogdanite()
        else:
            data = "Nothing found!"
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Allabogdanite":
                data_mineral = self.create_allabogdanite()
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
    def create_allabogdanite(self): # (Fe,Ni)2 P
        #
        name = "Abgd"
        #
        # Major elements
        phosphorus = PeriodicSystem(name="P").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["P", "Fe", "Ni"]
        w_Fe = round(rd.uniform(0.25, 0.75), 4)
        majors_data = np.array([["P", phosphorus[1], 1, phosphorus[2]], ["Fe", iron[1], round(2*w_Fe, 4), iron[2]],
                                ["Ni", nickel[1], round(2*(1-w_Fe), 4), nickel[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = 2*(w_Fe*iron[2] + (1-w_Fe)*nickel[2]) + phosphorus[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[5.811, 3.431], [], "hexagonal"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 3
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe*10**(6)])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Ni = CrystalPhysics([[5.873, 3.349], [], "hexagonal"])
        V_Ni = dataV_Ni.calculate_volume()
        Z_Ni = 3
        V_m_Ni = MineralChemistry().calculate_molar_volume(volume_cell=V_Ni, z=Z_Ni)
        dataRho_Ni = CrystalPhysics([molar_mass, Z_Ni, V_Ni*10**(6)])
        rho_Ni = dataRho_Ni.calculate_bulk_density()
        rho_e_Ni = wg(amounts=amounts, elements=element, rho_b=rho_Ni).calculate_electron_density()
        #
        V_m = w_Fe*V_m_Fe + (1-w_Fe)*V_m_Ni
        rho = w_Fe*rho_Fe + (1-w_Fe)*rho_Ni
        rho_e = w_Fe*rho_e_Fe + (1-w_Fe)*rho_e_Ni
        # Bulk modulus
        K_Fe = 216*10**9
        K_Ni = 194*10**9
        K = w_Fe*K_Fe + (1-w_Fe)*K_Ni
        # Shear modulus
        G_Fe = 89*10**9
        G_Ni = 42*10**9
        G = w_Fe*G_Fe + (1-w_Fe)*G_Ni
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