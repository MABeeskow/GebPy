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
    def __init__(self, traces_list=None, traces=None):
        self.traces = traces
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
        if self.traces == None and self.traces_list == None:
            M = round(silicon[2] + 2*oxygen[2], 3)
            w_Si = round(silicon[2]/M, 4)
            w_O = round(2*oxygen[2]/M, 4)
            element = [oxygen, silicon]
            composition = [w_O, w_Si]
            amounts = [1, 2]
        elif self.traces_list != None:
            x_list = []
            for i in range(len(self.traces_list)):
                condition = False
                while condition == False:
                    x = rd.uniform(0., 0.001)
                    M = round((1-x)*silicon[2] + x*aluminium[2] + 2*oxygen[2], 3)
                    w_O = round(2*oxygen[2]/M, 6)
                    w_Al = round(x*aluminium[2]/M, 6)
                    w_Si = round((1-x)*silicon[2]/M, 6)
                    if 1*10**(-6) <= w_Al <= 500*10**(-6):
                        element = [oxygen, aluminium, silicon]
                        composition = [w_O, w_Al, w_Si]
                        amounts = [x, 1-x, 2]
                        condition = True
                    else:
                        continue
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
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 2*10**14
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 1), round(G*10**(-9), 1), round(E*10**(-9), 1), round(nu, 4)])
        data.append([round(vP, 1), round(vS, 1), round(vPvS, 2)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data