#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		oxides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		04.06.2021

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.elements import elements
from modules.chemistry import PeriodicSystem
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg

# OXIDES
class Quartz(): # SiO2
    """ Class that generates geophysical and geochemical data of quartz"""
    #
    def __init__(self, traces_list=[], impurity="pure"):
        self.traces_list = traces_list
        self.impurity = impurity
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            minors = ["H", "Al", "Li", "Fe", "Ti", "Na", "Mg", "Ge"]
            n = rd.randint(1, len(minors))
            while len(self.traces_list) < n:
                selection = rd.choice(minors)
                if selection not in self.traces_list:
                    self.traces_list.append(selection)
                else:
                    continue
    #
    def create_quartz(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        # Minor elements
        aluminium = elements.Al(self)
        titanium = elements.Ti(self)
        lithium = elements.Li(self)
        #
        data = []
        #
        mineral = "Qz"
        #
        # Molar mass
        majors = [oxygen, silicon]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        x_hydrogen = 0
        for i in range(len(self.traces_list)):
            if self.traces_list[i] == "H":
                x_hydrogen = x_traces[i]
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            else:
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [(2.0-x_hydrogen), (1.0-(np.sum(x_traces)-x_hydrogen))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[4.9135, 5.4050], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 29*10**9
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
        GR = 0
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 2*10**14
        #
        data.append(mineral)
        data.append(round(M, 3))
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data

class Uraninite():  # UO2
    """ Class that generates geophysical and geochemical data of uraninite"""
    #
    def __init__(self, traces_list=[], impurity="pure"):
        self.traces_list = traces_list
        self.impurity = impurity
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            minors = ["Th", "Zr", "Pb", "Ra", "Ac", "Po", "Ce", "Y", "Er", "La"]
            n = rd.randint(1, len(minors))
            while len(self.traces_list) < n:
                selection = rd.choice(minors)
                if selection not in self.traces_list:
                    self.traces_list.append(selection)
                else:
                    continue
    #
    def create_uraninite(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        uranium = PeriodicSystem(name="U").get_data()
        #
        data = []
        #
        mineral = "Urn"
        #
        # Molar mass
        majors = [oxygen, uranium]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        for i in range(len(self.traces_list)):
            M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
            weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [2.0, (1-np.sum(x_traces))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[5.4682], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 200*10**9
        # Shear modulus
        G = 89*10**9
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
        w_U = round(uranium[2]/M, 4)
        GR = w_U*1000000*8
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 28*10**(-8)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data

class Magnetite():  # Fe3O4
    """ Class that generates geophysical and geochemical data of magnetite"""
    #
    def __init__(self, traces_list=[], impurity="pure"):
        self.traces_list = traces_list
        self.impurity = impurity
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            minors = ["Mg", "Zn", "Mn", "Ni", "Cr", "Ti", "V", "Al"]
            n = rd.randint(1, len(minors))
            while len(self.traces_list) < n:
                selection = rd.choice(minors)
                if selection not in self.traces_list:
                    self.traces_list.append(selection)
                else:
                    continue
    #
    def create_magnetite(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        #
        data = []
        #
        mineral = "Mag"
        #
        # Molar mass
        majors = [oxygen, iron]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        for i in range(len(self.traces_list)):
            M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
            weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [4.0, (3-np.sum(x_traces))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[8.396], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 176*10**9
        # Shear modulus
        G = 64*10**9
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 2850
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data