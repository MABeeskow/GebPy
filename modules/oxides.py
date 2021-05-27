#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		oxides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.05.2021

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
    def __init__(self, traces_list=[]):
        self.traces_list = traces_list
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
        for i in range(len(self.traces_list)):
            M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
            weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            #weights.append([int(traces[i][1]), x_traces[i]])
        x_majors = [2.0, (1-np.sum(x_traces))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
            #weights.append([majors[i][1], x_majors[i]])
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