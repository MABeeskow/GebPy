#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		fluids.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		16.02.2020

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals
from modules.elements import elements
from modules.geophysics import BoreholeGeophysics as bg

class Gas:
    #
    def __init__(self):
        pass
    #
    def air(self):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        nitrogen = elements.N(self)
        oxygen = elements.O(self)
        argon = elements.Ar(self)
        element = [nitrogen, oxygen, argon]
        #
        data = []
        #
        fluid = "air"
        #
        # Molar mass
        w_N = 0.78
        w_O = 0.21
        w_Ar = 0.01
        M = round(w_N*nitrogen[2] + w_O*oxygen[2] + w_Ar*argon[2], 3)
        weights = [w_N, w_O, w_Ar]
        # Density
        rho = 1.2041
        # Bulk modulus
        K = round(131126.49, 3)
        # Shear modulus
        G = 0.0
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = np.inf
        # P-wave velocity
        vP = ((K + 4/3*G)/(rho))**(0.5)
        # S-wave velocity
        vS = ((G)/(rho))**(0.5)
        # Gamma ray
        GR = 0
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(fluid)
        data.append(round(M,2))
        data.append(round(rho,1))
        data.append([round(K*10**(-9),2), round(G*10**(-9),2), round(E*10**(-9),2), round(nu,2), round(vPvS,2)])
        data.append([round(vP,1), round(vS,1)])
        data.append([round(GR,2), round(PE,2), round(U,2)])
        data.append(weights)
        #
        return data
#
class Water:
    #
    def __init__(self):
        pass
    #
    def seawater(self):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemNa = elements.Na(self)
        chemCl = elements.Cl(self)
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        chemWater = Water.water("")
        chemHalite = minerals.halides.halite("")
        #
        chemSeawater = []
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        saltAmount = 0.035
        M = (1-saltAmount)*chemWater[1]+saltAmount*(chemHalite[1])                      # Molar mass
        x_Water = round((1-saltAmount)*chemWater[1]/M, 4)
        x_Salt = round(saltAmount*chemHalite[1]/M, 4)
        weights = [x_Water, x_Salt]
        chemSeawater.append(M)
        rho = round(np.random.normal(997, 0),1)                                      # Density
        chemSeawater.append(rho)
        K = round(chemWater[2],2)            # Bulk modulus
        chemSeawater.append(K)
        G = round(chemWater[3],2)             # Shear modulus
        chemSeawater.append(G)
        vP = ((K*10**9 + 4/3 * G*10**9)/(rho))**(0.5)   # P-wave velocity
        chemSeawater.append(round(vP,1))
        vS = ((G*10**9)/(rho))**(0.5)                   # S-wave velocity
        chemSeawater.append(round(vS,1))
        GR = 0                                          # Gamma ray
        chemSeawater.append(GR)
        #
        data = []
        #
        fluid = "H2O"
        #
        # Molar mass
        M = (1-saltAmount)*chemWater[1]+saltAmount*(chemHalite[1])
        x_Water = round((1-saltAmount)*chemWater[1]/M, 4)
        x_Salt = round(saltAmount*chemHalite[1]/M, 4)
        weights = [x_Water, x_Salt]
        # Density
        dataV = minerals.CrystalPhysics([[4.51, 7.35], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = minerals.CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 2.1*10**9
        # Shear modulus
        G = 0.0
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = np.inf
        # P-wave velocity
        vP = ((K + 4/3*G)/(rho))**(0.5)
        # S-wave velocity
        vS = ((G)/(rho))**(0.5)
        # Gamma ray
        GR = 0
        # Photoelectricity
        PE = x_Water*chemWater[5][1] + x_Salt*chemHalite[5][1]
        U = PE*rho*10**(-3)
        #
        data.append(fluid)
        data.append(round(M,2))
        data.append(round(rho,1))
        data.append([round(K*10**(-9),2), round(G*10**(-9),2), round(E*10**(-9),2), round(nu,2), round(vPvS,2)])
        data.append([round(vP,1), round(vS,1)])
        data.append([round(GR,2), round(PE,2), round(U,2)])
        #
        return data
    #
    def water(self):   # H2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_H = elements.H(self)
        chem_O = elements.O(self)
        element = [chem_H, chem_O]
        #
        data = []
        #
        mineral = "H2O"
        #
        # Molar mass
        M = round(2*chem_H[2] + chem_O[2],3)
        x_H = round(2*chem_H[2]/M, 4)
        x_O = round(chem_O[2]/M, 4)
        weights = [x_H, x_O]
        # Density
        dataV = minerals.CrystalPhysics([[4.51, 7.35], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = minerals.CrystalPhysics([M, 4, V])
        rho_calc = dataRho.calculate_bulk_density()
        rho = rho_calc*1000/rho_calc
        # Bulk modulus
        K = 2.1*10**9
        # Shear modulus
        G = 0.0
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = np.inf
        # P-wave velocity
        vP = ((K + 4/3*G)/(rho))**(0.5)
        # S-wave velocity
        vS = ((G)/(rho))**(0.5)
        # Gamma ray
        GR = 0
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
    #
class Hydrocarbons:
    #
    def __init__(self):
        pass
    #
    def oil(self):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        carbon = elements.C(self)
        nitrogen = elements.N(self)
        oxygen = elements.O(self)
        sulfur = elements.S(self)
        element = [hydrogen, carbon, nitrogen, oxygen, sulfur]
        #
        data = []
        #
        fluid = "oil"
        #
        # Molar mass
        condition = False
        condition_2 = False
        while condition == False and condition_2 == False:
            w_C = round(rd.uniform(0.83, 0.85), 4)
            w_H = round(rd.uniform(0.10, 0.14), 4)
            w_N = round(rd.uniform(0.001, 0.02), 4)
            w_O = round(rd.uniform(0.0005, 0.015), 4)
            w_S = round(1-w_C-w_H-w_N-w_O, 4)
            w = w_H + w_C + w_N + w_O + w_S
            if w == 1.0:
                M = round(w_H*hydrogen[2] + w_C*carbon[2] + w_N*nitrogen[2] + w_O*oxygen[2] + w_S*sulfur[2], 3)
                weights = [w_H, w_C, w_N, w_O, w_S]
                # Density
                V_m = round(w_H*hydrogen[3] + w_C*carbon[3] + w_N*nitrogen[3] + w_O*oxygen[3] + w_S*sulfur[3], 3)
                rho = round(M/1000/V_m/0.0135, 4)
                if 850 <= rho <= 950:
                    condition = True
                    condition_2 = True
                else:
                    condition = False
                    condition_2 = False
            else:
                condition = False
        # Bulk modulus
        K = round(w_H*hydrogen[5] + w_C*carbon[5] + w_N*nitrogen[5] + w_O*oxygen[5] + w_S*sulfur[5], 3)
        # Shear modulus
        G = 0.0
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = np.inf
        # P-wave velocity
        vP = ((K + 4/3*G)/(rho))**(0.5)
        # S-wave velocity
        vS = ((G)/(rho))**(0.5)
        # Gamma ray
        GR = 0
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(fluid)
        data.append(round(M,2))
        data.append(round(rho,1))
        data.append([round(K*10**(-9),2), round(G*10**(-9),2), round(E*10**(-9),2), round(nu,2), round(vPvS,2)])
        data.append([round(vP,1), round(vS,1)])
        data.append([round(GR,2), round(PE,2), round(U,2)])
        data.append(weights)
        #
        return data
    #
    def natural_gas(self):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        carbon = elements.C(self)
        nitrogen = elements.N(self)
        oxygen = elements.O(self)
        sulfur = elements.S(self)
        element = [hydrogen, carbon, nitrogen, oxygen, sulfur]
        #
        data = []
        #
        fluid = "gas"
        #
        # Molar mass
        condition = False
        condition_2 = False
        while condition == False and condition_2 == False:
            w_C = round(rd.uniform(0.83, 0.85), 4)
            w_H = round(rd.uniform(0.10, 0.14), 4)
            w_N = round(rd.uniform(0.001, 0.02), 4)
            w_O = round(rd.uniform(0.0005, 0.015), 4)
            w_S = round(1-w_C-w_H-w_N-w_O, 4)
            w = w_H + w_C + w_N + w_O + w_S
            if w == 1.0:
                M = round(w_H*hydrogen[2] + w_C*carbon[2] + w_N*nitrogen[2] + w_O*oxygen[2] + w_S*sulfur[2], 3)
                weights = [w_H, w_C, w_N, w_O, w_S]
                # Density
                V_m = round(w_H*hydrogen[3] + w_C*carbon[3] + w_N*nitrogen[3] + w_O*oxygen[3] + w_S*sulfur[3], 3)
                rho = round(M/1000/V_m/0.0145, 4)
                if 750 <= rho <= 850:
                    condition = True
                    condition_2 = True
                else:
                    continue
            else:
                condition = False
        # Bulk modulus
        K = round(w_H*hydrogen[5] + w_C*carbon[5] + w_N*nitrogen[5] + w_O*oxygen[5] + w_S*sulfur[5], 3)
        # Shear modulus
        G = 0.0
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = np.inf
        # P-wave velocity
        vP = ((K + 4/3*G)/(rho))**(0.5)
        # S-wave velocity
        vS = ((G)/(rho))**(0.5)
        # Gamma ray
        GR = 0
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(fluid)
        data.append(round(M,2))
        data.append(round(rho,1))
        data.append([round(K*10**(-9),2), round(G*10**(-9),2), round(E*10**(-9),2), round(nu,2), round(vPvS,2)])
        data.append([round(vP,1), round(vS,1)])
        data.append([round(GR,2), round(PE,2), round(U,2)])
        data.append(weights)
        #
        return data