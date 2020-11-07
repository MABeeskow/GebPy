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
from modules import minerals
from modules.elements import elements
from modules.geophysics import BoreholeGeophysics as bg

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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M,2))
        data.append(round(rho,1))
        data.append([round(K*10**(-9),2), round(G*10**(-9),2), round(E*10**(-9),2), round(nu,2), round(vPvS,2)])
        data.append([round(vP,1), round(vS,1)])
        data.append([round(GR,2), round(PE,2), round(U,2)])
        #
        return data
    #