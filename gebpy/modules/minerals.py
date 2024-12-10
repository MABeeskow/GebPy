#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		minerals.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		28.05.2021

# -----------------------------------------------

## MODULES
import numpy as np
from numpy import round
import random as rd
from modules.elements import elements
import scipy.constants
from scipy import stats
from modules.geophysics import BoreholeGeophysics as bg



u = scipy.constants.physical_constants["atomic mass constant"]
angstrom = scipy.constants.physical_constants["Angstrom star"]
avogadro = scipy.constants.physical_constants["Avogadro constant"]

##################
# CRYSTALPHYSICS #
##################
#
class CrystalPhysics:
    #
    def __init__(self, properties):
        self.properties = properties
    #
    def calculate_bulk_density(self):
        # properties = [ molar mass, formula unit, unit cell volume ]
        M = self.properties[0]  # in g/mol
        Z = self.properties[1]
        V = self.properties[2]  # in cm^3
        #
        # density rho in kg/m^3
        rho = (Z*M)/(V*avogadro[0])*1000
        #
        return rho
    #
    def calculate_electron_density(self):
        # properties = [ elements, amounts, bulk density ]
        Z = [self.properties[1][i]*self.properties[0][i][1]
             for i in range(len(self.properties[0]))][0]/np.sum(self.properties[1])
        A = [self.properties[1][i]*self.properties[0][i][2]
             for i in range(len(self.properties[0]))][0]/np.sum(self.properties[1])
        rho_b = self.properties[2]
        #
        rho_e = 2*Z/A * rho_b
        #
        return rho_e
    #
    def calculate_volume(self):
        # properties = [ list of lattice lengths, list of lattice angles, crystal system ]
        lenghts = self.properties[0]    # in angstrom
        angles = self.properties[1]     # in degree
        crystalsystem = self.properties[2]
        #
        if crystalsystem == "cubic":
            a = lenghts[0]*10**(-8)
            V = a**3
            return V
        elif crystalsystem == "tetragonal":
            a = lenghts[0]*10**(-8)
            c = lenghts[1]*10**(-8)
            V = a**2 * c
            return V
        elif crystalsystem == "hexagonal":
            a = lenghts[0]*10**(-10)
            c = lenghts[1]*10**(-10)
            angle = 60
            V = (a**2 * c)*np.sin(angle*np.pi/180)
            return V
        elif crystalsystem == "trigonal":
            a = lenghts[0]*10**(-8)
            c = lenghts[1]*10**(-8)
            angle = 60
            V = (a**2 * c)*np.sin(angle*np.pi/180)
            return V
        elif crystalsystem == "orthorhombic":
            a = lenghts[0]*10**(-8)
            b = lenghts[1]*10**(-8)
            c = lenghts[2]*10**(-8)
            V = a * b * c
            return V
        elif crystalsystem == "monoclinic":
            a = lenghts[0]*10**(-8)
            b = lenghts[1]*10**(-8)
            c = lenghts[2]*10**(-8)
            beta = angles[0]
            V = (a * b * c)*np.sin(beta*np.pi/180)
            return V
        elif crystalsystem == "triclinic":
            a = lenghts[0]*10**(-8)
            b = lenghts[1]*10**(-8)
            c = lenghts[2]*10**(-8)
            alpha = angles[0]
            beta = angles[1]
            gamma = angles[2]
            V = (a * b * c)*(1 - np.cos(alpha*np.pi/180)**2 - np.cos(beta*np.pi/180)**2 - np.cos(gamma*np.pi/180)**2 + 2*(abs(np.cos(alpha*np.pi/180)*np.cos(beta*np.pi/180)*np.cos(gamma*np.pi/180))))**(0.5)
            return V
        #
#
#####################
# ORGANIC COMPOUNDS #
#####################
#
class Organics:
    #
    def __init__(self):
        pass
    #
    def carbohydrates(self):
        # CHEMISTRY
        hydrogen = elements.H(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        name = "carbohydrates"
        #
        # Molar mass
        w_H = 0.06
        w_C = 0.44
        w_O = 0.50
        M = round(w_H*hydrogen[2] + w_C*carbon[2] + w_O*oxygen[2], 3)
        composition = [w_H, w_C, w_O]
        # Density
        rho = 1586
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
        GR = 0
        # Photoelectricity
        PE = round(((0.06*hydrogen[1] + 0.44*carbon[1] + 0.5*oxygen[1])/10)**3.6, 3)
        #
        data.append(name)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        data.append(composition)
        #
        return data
    #
    def lignin(self):
        # CHEMISTRY
        hydrogen = elements.H(self)
        carbon = elements.C(self)
        nitrogen = elements.N(self)
        oxygen = elements.O(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        name = "lignin"
        #
        # Molar mass
        w_H = 0.06
        w_C = 0.63
        w_N = 0.003
        w_O = 0.31
        w_S = 0.001
        M = round(w_H*hydrogen[2] + w_C*carbon[2] + w_N*nitrogen[2] + w_O*oxygen[2] + w_S*sulfur[2], 3)
        composition = [w_H, w_C, w_N, w_O, w_S]
        # Density
        rho = 680
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
        GR = 0
        # Photoelectricity
        PE = round(((0.06*hydrogen[1] + 0.63*carbon[1] + 0.003*nitrogen[1] + 0.31*oxygen[1] + 0.001*sulfur[1])/10)**3.6, 3)
        #
        data.append(name)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        data.append(composition)
        #
        return data
    #
    def lipid(self):
        # CHEMISTRY
        hydrogen = elements.H(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        name = "lipid"
        #
        # Molar mass
        w_H = 0.10
        w_C = 0.80
        w_O = 0.10
        M = round(w_H*hydrogen[2] + w_C*carbon[2] + w_O*oxygen[2], 3)
        composition = [w_H, w_C, w_O]
        # Density
        rho = 850
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
        GR = 0
        # Photoelectricity
        PE = round(((0.10*hydrogen[1] + 0.80*carbon[1] + 0.10*oxygen[1])/10)**3.6, 3)
        #
        data.append(name)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        data.append(composition)
        #
        return data
    #
    def organic_matter(self):
        # CHEMISTRY
        carbohydrates = Organics.carbohydrates("")
        lignin = Organics.lignin("")
        lipid = Organics.lipid("")
        hydrogen = elements.H(self)
        carbon = elements.C(self)
        nitrogen = elements.N(self)
        oxygen = elements.O(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        name = "Org"
        #
        # Molar
        condition = False
        while condition == False:
            w_ch = round(rd.uniform(0.4, 0.6), 4)
            w_lg = round(rd.uniform(0.2, float(1-w_ch)), 4)
            w_lp = round(1 - w_ch - w_lg, 4)
            if w_ch+w_lg+w_lp == 1.0:
                condition = True
        M = round(w_ch*carbohydrates[1] + w_lg*lignin[1] + w_lp*lipid[1], 3)
        w_H = round((w_ch*carbohydrates[6][0]*hydrogen[2] + w_lg*lignin[6][0]*hydrogen[2] + w_lp*lipid[6][0]*hydrogen[2])/M, 4)
        w_C = round((w_ch*carbohydrates[6][1]*carbon[2] + w_lg*lignin[6][1]*carbon[2] + w_lp*lipid[6][1]*carbon[2])/M, 4)
        w_N = round((w_lg*lignin[6][2]*nitrogen[2])/M, 4)
        w_O = round((w_ch*carbohydrates[6][2]*oxygen[2] + w_lg*lignin[6][3]*oxygen[2] + w_lp*lipid[6][2]*oxygen[2])/M, 4)
        w_S = round(1-w_H-w_C-w_N-w_O, 4)
        composition = [w_H, w_C, w_N, w_O, w_S]
        # Density
        rho = w_ch*carbohydrates[2] + w_lg*lignin[2] + w_lp*lipid[2]
        # Bulk modulus
        K = (w_ch*carbohydrates[3][0] + w_lg*lignin[3][0] + w_lp*lipid[3][0])*10**9
        # Shear modulus
        G = (w_ch*carbohydrates[3][1] + w_lg*lignin[3][1] + w_lp*lipid[3][1])*10**9
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
        PE = round(w_ch*carbohydrates[5][1] + w_lg*lignin[5][1] + w_lp*lipid[5][1], 3)
        #
        data.append(name)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        data.append(composition)
        #
        return data
#
############
# MINERALS #
############
#
#
class natives:
    #
    def __init__(self):
        pass
    #
    def organic_matter(self):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        carbon = elements.C(self)
        #
        data = []
        #
        mineral = "Org"
        #
        # Molar mass
        M = round(carbon[2], 3)
        w_C = carbon[2]/M
        composition = [w_C]
        # Density
        dataV = CrystalPhysics([[2.467, 6.71], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 119*10**9
        # Shear modulus
        G = 93*10**9
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
        PE = round((carbon[1]/10)**3.6, 3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        data.append(composition)
        #
        return data
    #
    def sulfur(self):   # S
        # CHEMISTRY
        sulfur = elements.S(self)
        #
        data = []
        #
        name = "S"
        #
        # Molar mass
        w_S = 1.0
        M = round(8*w_S*sulfur[2], 3)
        composition = [w_S]
        # Density
        dataV = CrystalPhysics([[10.45, 12.845, 24.46], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 16, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 16.87*10**9
        # Shear modulus
        G = 9.14*10**9
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
        PE = round(((sulfur[1])/10)**3.6, 3)
        #
        data.append(name)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        data.append(composition)
        #
        return data
    #
    def arsenic(self):   # As
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        arsenic = elements.As(self)
        #
        data = []
        #
        mineral = "As"
        #
        # Molar mass
        M = round(arsenic[2], 3)
        # Density
        dataV = CrystalPhysics([[3.768, 10.574], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 22*10**9
        # Shear modulus
        G = 2.78*10**9
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
        PE = round((arsenic[1]/10)**3.6, 3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        #
        return data
    #
    def bismuth(self):  # Bi
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Bi = elements.Bi(self)
        #
        data = []
        #
        mineral = "Bi"
        #
        # Molar mass
        M = round(chem_Bi[2], 3)
        # Density
        dataV = CrystalPhysics([[4.537, 11.838], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 29*10**9
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
        GR = 0
        # Photoelectricity
        PE = round((chem_Bi[1]/10)**3.6, 3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        #
        return data
    #
    def antimony(self):   # Sb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        antimony = elements.Sb(self)
        #
        data = []
        #
        mineral = "Sb"
        #
        # Molar mass
        M = round(antimony[2], 3)
        # Density
        dataV = CrystalPhysics([[4.299, 11.25], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 36*10**9    # estimated
        # Shear modulus
        G = 24*10**9    # estimated
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
        PE = round((antimony[1]/10)**3.6, 3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        #
        return data
#
class oxides:
    #
    def __init__(self):
        pass
    #
    def quartz(self, traces=None, tr_Al=False, tr_Ti=False, tr_Li=False):   # SiO2
        # Major elements
        oxygen = elements.O(self)
        silicon = elements.Si(self)
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
        if traces == None and tr_Al == False and tr_Ti == False and tr_Li == False:
            M = round(silicon[2] + 2*oxygen[2], 3)
            w_Si = round(silicon[2]/M, 4)
            w_O = round(2*oxygen[2]/M, 4)
            element = [oxygen, silicon]
            composition = [w_O, w_Si]
            amounts = [1, 2]
        elif tr_Al == True:
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
        elif tr_Ti == True:
            condition = False
            while condition == False:
                x = rd.uniform(0., 0.001)
                M = round((1-x)*silicon[2] + x*titanium[2] + 2*oxygen[2], 3)
                w_O = round(2*oxygen[2]/M, 6)
                w_Ti = round(x*titanium[2]/M, 6)
                w_Si = round((1-x)*silicon[2]/M, 6)
                if 1*10**(-6) <= w_Ti <= 500*10**(-6):
                    element = [oxygen, titanium, silicon]
                    composition = [w_O, w_Ti, w_Si]
                    amounts = [x, 1-x, 2]
                    condition = True
                else:
                    continue
        elif tr_Li == True:
            condition = False
            while condition == False:
                x = rd.uniform(0., 0.001)
                M = round((1-x)*silicon[2] + x*lithium[2] + 2*oxygen[2], 3)
                w_O = round(2*oxygen[2]/M, 6)
                w_Li = round(x*lithium[2]/M, 6)
                w_Si = round((1-x)*silicon[2]/M, 6)
                if 1*10**(-6) <= w_Li <= 500*10**(-6):
                    element = [oxygen, lithium, silicon]
                    composition = [w_O, w_Li, w_Si]
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
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def uraninite(self):   # UO2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        uranium = elements.U(self)
        oxygen = elements.O(self)
        element = [uranium, oxygen]
        #
        data = []
        #
        mineral = "Urn"
        #
        # Molar mass
        M = round(uranium[2] + 2*oxygen[2], 3)
        x_U = round(uranium[2]/M, 4)
        x_O = round(2*oxygen[2]/M, 4)
        amounts = [1, 2]
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
        w_U = round(uranium[2]/M, 4)
        w_O = round(2*oxygen[2]/M, 4)
        weights = [w_U, w_O]
        compositon = [w_O, w_U]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
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
        data.append(compositon)
        #
        return data
    #
    def magnetite(self):   # Fe3O4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        oxygen = elements.O(self)
        element = [iron, oxygen]
        #
        data = []
        #
        mineral = "Mag"
        #
        # Molar mass
        M = round(3*iron[2] + 4*oxygen[2], 3)
        x_Fe = round(3*iron[2]/M, 4)
        x_O = round(4*oxygen[2]/M, 4)
        weights = [x_Fe, x_O]
        amounts = [3, 4]
        # Density
        dataV = CrystalPhysics([[8.397], [], "cubic"])
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
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
        #
        return data
    #
    def hematite(self):   # Fe2O3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        oxygen = elements.O(self)
        element = [iron, oxygen]
        #
        data = []
        #
        mineral = "Hem"
        #
        # Molar mass
        M = round(2*iron[2] + 3*oxygen[2], 3)
        w_Fe = round(2*iron[2]/M, 4)
        w_O = round(3*oxygen[2]/M, 4)
        weights = [w_Fe, w_O]
        composition = [w_O, w_Fe]
        amounts = [2, 3]
        # Density
        dataV = CrystalPhysics([[5.0317, 13.737], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 143.58*10**9
        # Shear modulus
        G = 53.43*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 5*10**6
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def aluminiumoxide(self):   # Al2O3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        aluminium = elements.Al(self)
        oxygen = elements.O(self)
        element = [aluminium, oxygen]
        #
        data = []
        #
        mineral = "Al2O3"
        #
        # Molar mass
        M = round(2*aluminium[2] + 3*oxygen[2], 3)
        x_Al = round(2*aluminium[2]/M, 4)
        x_O = round(3*oxygen[2]/M, 4)
        weights = [x_Al, x_O]
        amounts = [2, 3]
        # Density
        dataV = CrystalPhysics([[4.75, 12.98], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 232*10**9
        # Shear modulus
        G = 147*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def wuestite(self):   # FeO
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        oxygen = elements.O(self)
        element = [iron, oxygen]
        #
        data = []
        #
        mineral = "FeO"
        #
        # Molar mass
        M = round(iron[2] + oxygen[2], 3)
        x_Fe = round(iron[2]/M, 4)
        x_O = round(oxygen[2]/M, 4)
        weights = [x_Fe, x_O]
        amounts = [1, 1]
        # Density
        dataV = CrystalPhysics([[4.31], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 147*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def potassiumoxide(self):   # K2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        oxygen = elements.O(self)
        element = [potassium, oxygen]
        #
        data = []
        #
        mineral = "K2O"
        #
        # Molar mass
        M = round(2*potassium[2] + oxygen[2], 3)
        x_K = round(2*potassium[2]/M, 4)
        x_O = round(oxygen[2]/M, 4)
        weights = [x_K, x_O]
        amounts = [2, 1]
        # Density
        dataV = CrystalPhysics([[6.436], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 27*10**9
        # Shear modulus
        G = 12*10**9
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
        GR = potassium[2]/M*100*16
        # Photoelectricity
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def water(self):   # H2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        element = [hydrogen, oxygen]
        #
        chemWater = []
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        M = round(hydrogen[2] + 2*oxygen[2], 3)              # Molar mass
        chemWater.append(M)
        dataV = CrystalPhysics([[7.143, 7.604], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()                # Density
        chemWater.append(rho)
        K = round(np.random.normal(2.2, 0), 2)           # Bulk modulus
        chemWater.append(K)
        G = round(np.random.normal(0.0, 0), 2)           # Shear modulus
        chemWater.append(G)
        vP = ((K*10**9 + 4/3 * G*10**9)/rho)**0.5   # P-wave velocity
        chemWater.append(round(vP, 1))
        vS = ((G * 10**9)/rho)**0.5                   # S-wave velocity
        chemWater.append(round(vS, 1))
        GR = 0                                          # Gamma ray
        chemWater.append(GR)
        #
        return chemWater
    #
    def ice(self):   # H2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        element = [hydrogen, oxygen]
        #
        data = []
        #
        mineral = "Ice"
        #
        # Molar mass
        M = round(2*hydrogen[2] + oxygen[2], 3)
        x_H = round(2*hydrogen[2]/M, 4)
        x_O = round(oxygen[2]/M, 4)
        weights = [x_H, x_O]
        amounts = [2, 1]
        # Density
        dataV = CrystalPhysics([[4.51, 7.35], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 58.51*10**9
        # Shear modulus
        G = 33.18*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def chromite(self):   # FeCr2O4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        chem_Cr = elements.Cr(self)
        oxygen = elements.O(self)
        element = [iron, chem_Cr, oxygen]
        #
        data = []
        #
        mineral = "Chr"
        #
        # Molar mass
        M = round(iron[2] + 2*chem_Cr[2] + 4*oxygen[2], 3)
        x_Fe = round(iron[2]/M, 4)
        x_Cr = round(2*chem_Cr[2]/M, 4)
        x_O = round(4*oxygen[2]/M, 4)
        weights = [x_Fe, x_Cr, x_O]
        amounts = [1, 1, 4]
        # Density
        dataV = CrystalPhysics([[8.36], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 147.21*10**9
        # Shear modulus
        G = 55.77*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def spinel(self):   # MgAl2O4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        oxygen = elements.O(self)
        element = [magnesium, aluminium, oxygen]
        #
        data = []
        #
        mineral = "Spl"
        #
        # Molar mass
        M = round(magnesium[2] + 2*aluminium[2] + 4*oxygen[2], 3)
        x_Mg = round(magnesium[2]/M, 4)
        x_Al = round(2*aluminium[2]/M, 4)
        x_O = round(4*oxygen[2]/M, 4)
        weights = [x_Mg, x_Al, x_O]
        amounts = [1, 2, 4]
        # Density
        dataV = CrystalPhysics([[8.08], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 180*10**9
        # Shear modulus
        G = 96*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def boehmite(self):   # AlO(OH)
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        aluminium = elements.Al(self)
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        element = [aluminium, oxygen, hydrogen]
        #
        data = []
        #
        mineral = "Bhm"
        #
        # Molar mass
        M = round(aluminium[2] + oxygen[2] + (oxygen[2]+hydrogen[2]), 3)
        x_Al = round(aluminium[2]/M, 4)
        x_O = round(2*oxygen[2]/M, 4)
        x_H = round(hydrogen[2]/M, 4)
        weights = [x_Al, x_O, x_H]
        amounts = [1, 2, 1]
        # Density
        dataV = CrystalPhysics([[2.868, 12.227, 3.7], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 114*10**9
        # Shear modulus
        G = 82*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def diaspore(self):   # AlO(OH)
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        aluminium = elements.Al(self)
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        element = [aluminium, oxygen, hydrogen]
        #
        data = []
        #
        mineral = "Dsp"
        #
        # Molar mass
        M = round(aluminium[2] + oxygen[2] + (oxygen[2]+hydrogen[2]), 3)
        x_Al = round(aluminium[2]/M, 4)
        x_O = round(2*oxygen[2]/M, 4)
        x_H = round(hydrogen[2]/M, 4)
        weights = [x_Al, x_O, x_H]
        amounts = [1, 2, 1]
        # Density
        dataV = CrystalPhysics([[4.397, 9.421, 2.8439], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 174.02*10**9
        # Shear modulus
        G = 96.08*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def gibbsite(self):   # Al(OH)3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        aluminium = elements.Al(self)
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        element = [aluminium, oxygen, hydrogen]
        #
        data = []
        #
        mineral = "Gbs"
        #
        # Molar mass
        M = round(aluminium[2] + 3*(oxygen[2]+hydrogen[2]), 3)
        x_Al = round(aluminium[2]/M, 4)
        x_O = round(3*oxygen[2]/M, 4)
        x_H = round(3*hydrogen[2]/M, 4)
        weights = [x_Al, x_O, x_H]
        amounts = [1, 3, 3]
        # Density
        dataV = CrystalPhysics([[8.641, 5.07, 9.719], [94.566], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 174.02*10**9
        # Shear modulus
        G = 96.08*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
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
    def cuprite(self):   # Cu2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        copper = elements.Cu(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Cup"
        #
        # Molar mass
        M = round(2*copper[2] + 1*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[4.2696], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 111*10**9
        # Shear modulus
        G = 8*10**9
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
        element = [copper, oxygen]
        x_Cu = round(2*copper[2]/M, 4)
        x_O = round(1*oxygen[2]/M, 4)
        weights = [x_Cu, x_O]
        amounts = [2, 1]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
##############
# CARBONATES #
##############
#
class carbonates:
    #
    def __init__(self):
        pass
    #
    def calcite(self):   # CaCO3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Cal"
        #
        # Molar mass
        M = round(calcium[2] + carbon[2] + 3*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[4.99, 17.06], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 76*10**9
        # Shear modulus
        G = 32*10**9
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
        element = [calcium, carbon, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_C = round(carbon[2]/M, 4)
        w_O = round(3*oxygen[2]/M, 4)
        weights = [w_Ca, w_C, w_O]
        composition = [w_C, w_O, w_Ca]
        amounts = [1, 1, 3]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 2*10**12
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def dolomite(self):   # CaMg(CO3)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        calcium = elements.Ca(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Dol"
        #
        # Molar mass
        M = round(calcium[2] + magnesium[2] + 2*(carbon[2] + 3*oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[4.81, 16.01], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 89*10**9
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
        element = [calcium, magnesium, carbon, oxygen]
        x_Ca = round(calcium[2]/M, 4)
        x_Mg = round(magnesium[2]/M, 4)
        x_C = round(2*carbon[2]/M, 4)
        x_O = round(2*3*oxygen[2]/M, 4)
        weights = [x_Ca, x_Mg, x_C, x_O]
        composition = [x_C, x_O, x_Mg, x_Ca]
        amounts = [1, 1, 2, 3]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 5*10**3
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def siderite(self):   # Fe(CO3)
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Sd"
        #
        # Molar mass
        M = round(iron[2] + (carbon[2] + 3*oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[4.69, 15.38], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 124*10**9
        # Shear modulus
        G = 51*10**9
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
        element = [iron, carbon, oxygen]
        x_Fe = round(iron[2]/M, 4)
        x_C = round(carbon[2]/M, 4)
        x_O = round(3*oxygen[2]/M, 4)
        weights = [x_Fe, x_C, x_O]
        composition = [x_C, x_O, x_Fe]
        amounts = [1, 1, 3]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 70
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def magnesite(self):   # Mg(CO3)
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Mgs"
        #
        # Molar mass
        M = round(magnesium[2] + (carbon[2] + 3*oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[4.63, 15.03], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 105*10**9
        # Shear modulus
        G = 63*10**9
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
        element = [magnesium, carbon, oxygen]
        x_Mg = round(magnesium[2]/M, 4)
        x_C = round(carbon[2]/M, 4)
        x_O = round(3*oxygen[2]/M, 4)
        weights = [x_Mg, x_C, x_O]
        composition = [x_C, x_O, x_Mg]
        amounts = [1, 1, 3]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def ankerite(self, keyword="None"):   # Ca(Fe,Mg,Mn)(CO3)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        iron = elements.Fe(self)
        magnesium = elements.Mg(self)
        manganese = elements.Mn(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Ank"
        #
        # Molar mass
        if keyword == "None":
            x_Fe = round(rd.uniform(0.0, 1.0), 2)
            x_Mg = round(rd.uniform(0.0, 1-x_Fe), 2)
            x_Mn = round(1 - x_Fe - x_Mg, 2)
            M = round(calcium[2] + (x_Fe*iron[2] + x_Mg*magnesium[2] + x_Mn*manganese[2]) + 2*(carbon[2]+3*oxygen[2]), 3)
        elif keyword == "Fe":
            x_Fe = round(rd.uniform(0.5, 1.0), 2)
            x_Mg = round(rd.uniform(0.0, 1-x_Fe), 2)
            x_Mn = round(1 - x_Fe - x_Mg, 2)
            M = round(calcium[2] + (x_Fe*iron[2] + x_Mg*magnesium[2] + x_Mn*manganese[2]) + 2*(carbon[2]+3*oxygen[2]), 3)
        elif keyword == "Mg":
            x_Mg = round(rd.uniform(0.5, 1.0), 2)
            x_Fe = round(rd.uniform(0.0, 1-x_Mg), 2)
            x_Mn = round(1 - x_Fe - x_Mg, 2)
            M = round(calcium[2] + (x_Fe*iron[2] + x_Mg*magnesium[2] + x_Mn*manganese[2]) + 2*(carbon[2]+3*oxygen[2]), 3)
        elif keyword == "Mn":
            x_Mn = round(rd.uniform(0.5, 1.0), 2)
            x_Fe = round(rd.uniform(0.0, 1-x_Mn), 2)
            x_Mg = round(1 - x_Fe - x_Mn, 2)
            M = round(calcium[2] + (x_Fe*iron[2] + x_Mg*magnesium[2] + x_Mn*manganese[2]) + 2*(carbon[2]+3*oxygen[2]), 3)
        elif keyword == "Pure" or keyword == "pure":
            x_Mn = round(0.0, 2)
            x_Fe = round(1.0, 2)
            x_Mg = round(0.0, 2)
            M = round(calcium[2] + (x_Fe*iron[2] + x_Mg*magnesium[2] + x_Mn*manganese[2]) + 2*(carbon[2]+3*oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[4.83, 16.167], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 89*10**9
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
        element = [calcium, iron, magnesium, manganese, carbon, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_Fe = round(x_Fe*iron[2]/M, 4)
        w_Mg = round(x_Mg*magnesium[2]/M, 4)
        w_Mn = round(x_Mn*manganese[2]/M, 4)
        w_C = round(2*carbon[2]/M, 4)
        w_O = round(2*3*oxygen[2]/M, 4)
        weights = [w_Ca, w_Fe, w_Mg, w_Mn, w_C, w_O]
        composition = [w_C, w_O, w_Mg, w_Ca, w_Mn, w_Fe]
        amounts = [1, x_Fe, x_Mg, x_Mn, 2, 6]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x_Fe, x_Mg, x_Mn])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def malachite(self):   # Cu2(CO3)(OH)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        copper = elements.Cu(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Mal"
        #
        # Molar mass
        M = round(2*copper[2] + (carbon[2]+3*oxygen[2]) + 2*(oxygen[2]+hydrogen[2]), 3)
        # Density
        dataV = CrystalPhysics([[9.502, 11.974, 3.24], [98.75], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 118.77*10**9
        # Shear modulus
        G = 50.06*10**9
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
        element = [copper, carbon, oxygen, hydrogen]
        x_Cu = round(2*copper[2]/M, 4)
        x_C = round(1*carbon[2]/M, 4)
        x_O = round(5*oxygen[2]/M, 4)
        x_H = round(2*oxygen[2]/M, 4)
        weights = [x_Cu, x_C, x_O, x_H]
        amounts = [2, 1, 5, 2]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
    def aragonite(self):   # CaCO3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        carbon = elements.C(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Arg"
        #
        # Molar mass
        M = round(1*calcium[2] + 1*carbon[2] + 3*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[4.959, 7.968, 5.741], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 69*10**9
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
        GR = 0
        # Photoelectricity
        element = [calcium, carbon, oxygen]
        x_Ca = round(1*calcium[2]/M, 4)
        x_C = round(1*carbon[2]/M, 4)
        x_O = round(3*oxygen[2]/M, 4)
        weights = [x_Ca, x_C, x_O]
        composition = [x_C, x_O, x_Ca]
        amounts = [1, 1, 3]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
#
###################
# PHYLLOSILICATES #
###################
#
class phyllosilicates:
    #
    def __init__(self):
        pass
    #
    def illite(self): # (K,H3O) (Al,Mg,Fe)2 (Si,Al)4 O10 [(OH)2,(H2O)]
        # CHEMISTRY
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        potassium = elements.K(self)
        iron = elements.Fe(self)
        #
        data = []
        #
        mineral = "Ilt"
        #
        # Molar mass
        a = round(rd.uniform(0.8, 0.9), 4)
        b = round(rd.uniform(0.6, 0.7), 4)
        b2 = round(rd.uniform(0.0, float(1-b)), 4)
        c = round(rd.uniform(0.8, 1.0), 4)
        d = round(rd.uniform(0.0, 1.0), 4)
        M = round(a*potassium[2] + (1-a)*(3*hydrogen[2] + oxygen[2]) + 2*(b*aluminium[2] + b2*magnesium[2] + (1-b-b2)*iron[2]) + 4*(c*silicon[2] + (1-c)*aluminium[2]) + 10*oxygen[2] + d*(2*(oxygen[2]+hydrogen[2]) + (2*hydrogen[2]+oxygen[2])), 3)
        # Density
        dataV = CrystalPhysics([[5.18, 8.98, 10.32], [101.83], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = potassium[2]/M*100*16
        # Photoelectricity
        element = [hydrogen, oxygen, magnesium, aluminium, silicon, potassium, iron]
        w_H = round(((1-a)*3*hydrogen[2] + d*2*hydrogen[2] + d*2*hydrogen[2])/M, 4)
        w_O = round(((1-a)*oxygen[2] + 10*oxygen[2] + d*2*oxygen[2] + d*oxygen[2])/M, 4)
        w_Mg = round(2*b2*magnesium[2]/M, 4)
        w_Al = round((2*b*aluminium[2] + 4*(1-c)*aluminium[2])/M, 4)
        w_Si = round(4*c*silicon[2]/M, 4)
        w_K = round(a*potassium[2]/M, 4)
        w_Fe = round(2*(1-b-b2)*iron[2]/M, 4)
        composition = [w_H, w_O, w_Mg, w_Al, w_Si, w_K, w_Fe]
        data_rho_e = CrystalPhysics([element, composition, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 52.5
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def chamosite(self):   # (Fe,Mg)5Al(Si3Al)O10(OH,O)8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        iron = elements.Fe(self)
        #
        data = []
        #
        mineral = "Chl"
        #
        # Molar mass
        x = rd.uniform(0, 1)
        y = rd.uniform(0, 1)
        M = round(5*(x*iron[2] + (1-x)*magnesium[2]) + aluminium[2] + (3*silicon[2] + aluminium[2]) + 10*oxygen[2] + 8*(y*(oxygen[2]+hydrogen[2]) + (1-y)*oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[5.373, 9.306, 14.222], [97.88], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 165.02*10**9
        # Shear modulus
        G = 52.10*10**9
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
        element = [iron, magnesium, oxygen, aluminium, hydrogen, silicon]
        w_Fe = round(5*x*iron[2]/M, 4)
        w_Mg = round(5*(1-x)*magnesium[2]/M, 4)
        w_O = round((10+8*y+8*(1-y))*oxygen[2]/M, 4)
        w_Al = round((1+1)*aluminium[2]/M, 4)
        w_H = round(8*y*hydrogen[2]/M, 4)
        w_Si = round(3*silicon[2]/M, 4)
        weights = [w_Fe, w_Mg, w_O, w_Al, w_H, w_Si]
        composition = [w_H, w_O, w_Mg, w_Al, w_Si, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def vermiculite(self): # (Mg,Fe,Al)3 (Al,Si)4 O10 (OH)2 * 4(H2O)
        # CHEMISTRY
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        iron = elements.Fe(self)
        #
        data = []
        #
        mineral = "Vrm"
        #
        # Molar mass
        a = round(rd.uniform(0.55, 0.65), 4)
        a2 = round(rd.uniform(0.25, float(1-a)), 4)
        b = round(rd.uniform(0.7, 0.8), 4)
        M = round(3*(a*magnesium[2] + a2*iron[2] + (1-a-a2)*aluminium[2]) + 4*(b*aluminium[2] + (1-b)*silicon[2]) + 10*oxygen[2] + 2*(oxygen[2] + hydrogen[2]) + 4*(2*hydrogen[2]+oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[5.26, 9.23, 14.97], [96.82], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 175*10**9
        # Shear modulus
        G = 90*10**9
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
        element = [hydrogen, oxygen, magnesium, aluminium, silicon, iron]
        w_H = round((2*hydrogen[2] + 4*2*hydrogen[2])/M, 4)
        w_O = round((10*oxygen[2] + 2*oxygen[2] + 4*oxygen[2])/M, 4)
        w_Mg = round(3*a*magnesium[2]/M, 4)
        w_Al = round((3*(1-a-a2)*aluminium[2] + 4*(b*aluminium[2]))/M, 4)
        w_Si = round((4*(1-b)*silicon[2])/M, 4)
        w_Fe = round(3*a2*iron[2]/M, 4)
        composition = [w_H, w_O, w_Mg, w_Al, w_Si, w_Fe]
        data_rho_e = CrystalPhysics([element, composition, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 0
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def biotite(self):   # KMg2.5Fe0.5AlSi3O10(OH)1.75F0.25
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        magnesium = elements.Mg(self)
        iron = elements.Fe(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        flourine = elements.F(self)
        #
        data = []
        #
        #mineral = "Bt"
        #
        # Molar mass
        M = round(potassium[2]+2.5*magnesium[2]+0.5*iron[2]+aluminium[2]+3*silicon[2]+10*oxygen[2]+1.75*(oxygen[2]+hydrogen[2])+0.25*flourine[2], 3)
        # Density
        rho = round(np.random.normal(3050, 10), 1)
        # Bulk modulus
        K = round(np.random.normal(55.35, 2.5), 2)
        # Shear modulus
        G = round(np.random.normal(27.23, 5), 2)
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K*10**9 + 4/3 * G*10**9)/rho)**0.5
        # S-wave velocity
        vS = ((G * 10**9)/rho)**0.5
        # Gamma ray
        GR = potassium[2]/M*100*16
        #
        #data.append(mineral)
        data.append(M)
        data.append(rho)
        data.append(K)
        data.append(G)
        data.append(round(vP, 1))
        data.append(round(vS, 1))
        data.append(GR)
        #
        return data
    #
    def annite(self):   # KFe3AlSi3O10(OH,F)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        iron = elements.Fe(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        flourine = elements.F(self)
        #
        data = []
        #
        mineral = "Ann"
        #
        # Molar mass
        x = rd.uniform(0, 1)
        M = round(potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(x*(oxygen[2]+hydrogen[2]) + (1-x)*flourine[2]), 3)
        # Density
        dataV = CrystalPhysics([[5.3860, 9.3241, 10.2683], [100.63], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 55.35*10**9
        # Shear modulus
        G = 27.23*10**9
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
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append(round(GR, 2))
        #
        return data
    #
    def muscovite(self):   # KAl2(AlSi3O10)(OH)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        flourine = elements.F(self)
        #
        data = []
        #
        mineral = "Ms"
        #
        # Molar mass
        M = round(potassium[2]+3*aluminium[2]+3*silicon[2]+10*oxygen[2]+1.8*(oxygen[2]+hydrogen[2])+0.2*flourine[2], 3)
        # Density
        dataV = CrystalPhysics([[5.19, 9.03, 20.05], [95.5], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 112.64*10**9
        # Shear modulus
        G = 68.33*10**9
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
        element = [potassium, aluminium, silicon, oxygen, hydrogen, flourine]
        w_K = round(potassium[2]/M, 4)
        w_Al = round(3*aluminium[2]/M, 4)
        w_Si = round(3*silicon[2]/M, 4)
        w_O = round((10+1.8)*oxygen[2]/M, 4)
        w_H = round(1.8*hydrogen[2]/M, 4)
        w_F = round(0.2*flourine[2]/M, 4)
        weights = [w_K, w_Al, w_Si, w_O, w_H, w_F]
        composition = [w_H, w_O, w_F, w_Al, w_Si, w_K]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        # Electrical resistivity
        p = 5*10**13
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def glauconite(self):   # K0.6Na0.05Fe1.5Mg0.4Al0.3Si3.8O10(OH)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        sodium = elements.Na(self)
        magnesium = elements.Mg(self)
        iron = elements.Fe(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Glt"
        #
        # Molar mass
        M = round(0.6*potassium[2]+0.05*sodium[2]+1.5*iron[2]+0.4*magnesium[2]+0.3*aluminium[2]+3.8*silicon[2]+10*oxygen[2]+2*(oxygen[2]+hydrogen[2]), 3)
        # Density
        rho = round(np.random.normal(2640, 10), 1)
        # Bulk modulus
        K = 112.64*10**9
        # Shear modulus
        G = 68.33*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3 * G)/(rho))**(0.5)
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        GR = potassium[2]/M*100*16
        # Photoelectricity
        element = [potassium, sodium, magnesium, iron, aluminium, silicon, oxygen, hydrogen]
        w_K = round(0.6*potassium[2]/M, 4)
        w_Na = round(0.05*sodium[2]/M, 4)
        w_Mg = round(0.4*magnesium[2]/M, 4)
        w_Fe = round(1.5*iron[2]/M, 4)
        w_Al = round(0.3*aluminium[2]/M, 4)
        w_Si = round(3.8*silicon[2]/M, 4)
        w_O = round((10+2)*oxygen[2]/M, 4)
        w_H = round(2*hydrogen[2]/M, 4)
        weights = [w_K, w_Na, w_Mg, w_Fe, w_Al, w_Si, w_O, w_H]
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
    def kaolinite(self):   # Al2(OH)4Si2O5
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        aluminium = elements.Al(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        silicon = elements.Si(self)
        #
        data = []
        #
        mineral = "Kln"
        #
        # Molar mass
        M = round(2*aluminium[2] + 4*(oxygen[2]+hydrogen[2]) + 2*silicon[2] + 5*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[5.13, 8.89, 7.25], [90.0, 104.5, 89.8], "triclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [aluminium, oxygen, hydrogen, silicon]
        w_Al = round(2*aluminium[2]/M, 4)
        w_O = round((4+5)*oxygen[2]/M, 4)
        w_H = round(4*hydrogen[2]/M, 4)
        w_Si = round(2*silicon[2]/M, 4)
        weights = [w_Al, w_O, w_H, w_Si]
        composition = [w_H, w_O, w_Al, w_Si]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        # Electrical resistivity
        p = 52.5
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        data.append(composition)
        #
        return data
    #
    def clinochlore(self):   # (Mg,Fe)5Al(Si3Al)O10(OH)8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        iron = elements.Fe(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Chl"
        #
        # Molar mass
        condition = False
        while condition == False:
            x = round(rd.uniform(0, 1), 2)
            M = round(5*(x*magnesium[2]+(1-x)*iron[2]) + aluminium[2] + (3*silicon[2]+aluminium[2]) + 10*oxygen[2] + 8*(oxygen[2]+hydrogen[2]), 3)
            # Density
            dataV = CrystalPhysics([[5.3, 9.3, 14.3], [97], "monoclinic"])
            V = dataV.calculate_volume()
            dataRho = CrystalPhysics([M, 2, V])
            rho = dataRho.calculate_bulk_density()
            if 2136 < rho < 3099:
                condition = True
            else:
                condition = False
        # Bulk modulus
        if 2680 <= rho <= 2840:
            K = (127.54 + (165.02-127.54)/(2840-2680)*(rho-2680))*10**9
        elif rho > 2136:
            a = (165.02*10**9-127.54*10**9)/(2840-2680)
            b = 127.54*10**9 - a *2680
            K = (a *rho + b)
        # Shear modulus
        if 2680 <= rho <= 2840:
            G = (84.20 + (52.10-84.20)/(2840-2680)*(rho-2680))*10**9
        elif rho < 3099:
            a = (52.10*10**9-84.20*10**9)/(2840-2680)
            b = 84.20*10**9 - a *2680
            G = (a *rho + b)
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
        element = [magnesium, iron, aluminium, silicon, oxygen, hydrogen]
        w_Mg = round(5*x*magnesium[2]/M, 4)
        w_Fe = round(5*(1-x)*iron[2]/M, 4)
        w_Al = round((1+1)*aluminium[2]/M, 4)
        w_Si = round(3*silicon[2]/M, 4)
        w_O = round((10+8)*oxygen[2]/M, 4)
        w_H = round(8*hydrogen[2]/M, 4)
        weights = [w_Mg, w_Fe, w_Al, w_Si, w_O, w_H]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append([round(M, 2), x])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
    #
    def montmorillonite(self):   # (Na,Ca)0.3(Al,Mg)2Si4O10(OH2)(H2O)10
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        calcium = elements.Ca(self)
        aluminium = elements.Al(self)
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Mnt"
        #
        # Molar mass
        condition = False
        while condition == False:
            x = round(rd.uniform(0.6, 0.7), 2)
            y = round(rd.uniform(0.9, 1), 2)
            n = rd.randint(10,12)
            M = round(0.3*(x*sodium[2]+(1-x)*calcium[2]) + 2*(y*aluminium[2]+(1-y)*magnesium[2]) + 4*silicon[2] + 10*oxygen[2] + (oxygen[2]+2*hydrogen[2]) + n*(2*hydrogen[2]+oxygen[2]), 2)
            # Density
            dataV = CrystalPhysics([[5.17, 8.94, 9.95], [99.9], "monoclinic"])
            V = dataV.calculate_volume()
            dataRho = CrystalPhysics([M, 1, V])
            rho = dataRho.calculate_bulk_density()
            if rho > 1807:
                condition = True
            else:
                condition = False
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
        GR = 0
        # Photoelectricity
        element = [sodium, calcium, aluminium, magnesium, silicon, oxygen, hydrogen]
        w_Na = round(0.3*x*sodium[2]/M, 4)
        w_Ca = round(0.3*(1-x)*calcium[2]/M, 4)
        w_Al = round(2*y*aluminium[2]/M, 4)
        w_Mg = round(2*(1-y)*magnesium[2]/M, 4)
        w_Si = round(4*silicon[2]/M, 4)
        w_O = round((10+1+n)*oxygen[2]/M, 4)
        w_H = round((2+n*2)*hydrogen[2]/M, 4)
        weights = [w_Na, w_Ca, w_Al, w_Mg, w_Si, w_O, w_H]
        composition = [w_H, w_O, w_Na, w_Mg, w_Al, w_Si, w_Ca]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append([M, x, y, n])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def clay(self):  # Ilt + Mnt + Kln
        # Elements
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        sodium = elements.Na(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        potassium = elements.K(self)
        calcium = elements.Ca(self)
        iron = elements.Fe(self)
        # Minerals
        ilt = phyllosilicates.illite("")
        mnt = phyllosilicates.montmorillonite("")
        kln = phyllosilicates.kaolinite("")
        #
        data = []
        #
        mineral = "Clay"
        #
        # Molar mass
        x = rd.uniform(0.0, 1.0)
        y = rd.uniform(0.0, float(1-x))
        M = x*ilt[1] + y*mnt[1][0] + (1-x-y)*kln[1]
        # Density
        rho = x*ilt[2] + y*mnt[2] + (1-x-y)*kln[2]
        # Bulk modulus
        K = (x*ilt[3][0] + y*mnt[3][0] + (1-x-y)*kln[3][0])*10**9
        # Shear modulus
        G = (x*ilt[3][1] + y*mnt[3][1] + (1-x-y)*kln[3][1])*10**9
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
        GR = x*ilt[5][0] + y*mnt[5][0] + (1-x-y)*kln[5][0]
        # Photoelectricity
        element = [hydrogen, oxygen, sodium, magnesium, aluminium, silicon, potassium, calcium, iron]
        w_H = round(x*ilt[6][0] + y*mnt[6][0] + (1-x-y)*kln[6][0], 4)
        w_O = round(x*ilt[6][1] + y*mnt[6][1] + (1-x-y)*kln[6][1], 4)
        w_Na = round(y*mnt[6][2], 4)
        w_Mg = round(x*ilt[6][2] + y*mnt[6][3], 4)
        w_Al = round(x*ilt[6][3] + y*mnt[6][4] + (1-x-y)*kln[6][2], 4)
        w_Si = round(x*ilt[6][4] + y*mnt[6][5] + (1-x-y)*kln[6][3], 4)
        w_K = round(x*ilt[6][5], 4)
        w_Ca = round(y*mnt[6][6], 4)
        w_Fe = round(x*ilt[6][6], 4)
        composition = [w_H, w_O, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe]
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
#
#################
# BIOTITE GROUP #
#################
#
class Biotites:
    #
    def __init__(self):
        pass
    #
    def biotite_group(self, keyword="None"):  # K [Mg(a)Fe(1-a)]3 Al(b)Fe(1-b) [Al(c)Si(1-c)]3 O10 [OH(d)F(1-d)]2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        flourine = elements.F(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        potassium = elements.K(self)
        iron = elements.Fe(self)
        element = [hydrogen, oxygen, flourine, magnesium, aluminium, silicon, potassium, iron]
        #
        self.keyword = keyword
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        if self.keyword == "None":
            x = round(rd.uniform(0, 1), 2)
            y = round(rd.uniform(0, 1), 2)
            z = round(rd.uniform(0, 1), 2)
            mineral = "Bt"
            #
            # Molar mass
            M = round(potassium[2] + (3*(1-y)+2*y)*(x*magnesium[2]+(1-x)*iron[2]) + aluminium[2] + (y*aluminium[2]+(3*(1-y)+2*y)*silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            # Density
            M_Phl = round(potassium[2] + 3*magnesium[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Ann = round(potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Sdp = round(potassium[2] + 2*iron[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Eas = round(potassium[2] + 2*magnesium[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            dataV_Phl = CrystalPhysics([[5.3078, 9.1901, 10.1547], [100.08], "monoclinic"])
            V_Phl = dataV_Phl.calculate_volume()
            dataRho_Phl = CrystalPhysics([M_Phl, 2, V_Phl])
            rho_Phl = dataRho_Phl.calculate_bulk_density()
            dataV_Ann = CrystalPhysics([[5.3860, 9.3241, 10.2683], [100.63], "monoclinic"])
            V_Ann = dataV_Ann.calculate_volume()
            dataRho_Ann = CrystalPhysics([M_Ann, 2, V_Ann])
            rho_Ann = dataRho_Ann.calculate_bulk_density()
            dataV_Sdp = CrystalPhysics([[5.348, 9.261, 10.263], [100.19], "monoclinic"])
            V_Sdp = dataV_Sdp.calculate_volume()
            dataRho_Sdp = CrystalPhysics([M_Sdp, 2, V_Sdp])
            rho_Sdp = dataRho_Sdp.calculate_bulk_density()
            dataV_Eas = CrystalPhysics([[5.27, 9.13, 10.15], [99.54], "monoclinic"]) # estimated!
            V_Eas = dataV_Eas.calculate_volume()
            dataRho_Eas = CrystalPhysics([M_Eas, 2, V_Eas])
            rho_Eas = dataRho_Eas.calculate_bulk_density()
            rho = x*(y*rho_Phl+(1-y)*rho_Eas) + (1-x)*(y*rho_Ann+(1-y)*rho_Sdp)
            # Bulk modulus
            K_Phl = 103.78*10**9
            K_Sdp = 110.91*10**9
            K_Ann = 114.72*10**9 # estimated!
            K_Eas = 100.44*10**9 # estimated!
            K = x*(y*K_Phl+(1-y)*K_Eas) + (1-x)*(y*K_Ann+(1-y)*K_Sdp)
            # Shear modulus
            G_Phl = 61.69*10**9
            G_Sdp = 59.61*10**9
            G_Ann = 58.61*10**9
            G_Eas = 62.77*10**9
            G = x*(y*G_Phl+(1-y)*G_Eas) + (1-x)*(y*G_Ann+(1-y)*G_Sdp)
            # Young's modulus
            E_Phl = (9*K_Phl*G_Phl)/(3*K_Phl + G_Phl)
            E_Sdp = (9*K_Sdp*G_Sdp)/(3*K_Sdp + G_Sdp)
            E_Ann = (9*K_Ann*G_Ann)/(3*K_Ann + G_Ann) # estimated!
            E_Eas = (9*K_Eas*G_Eas)/(3*K_Eas + G_Eas) # estimated!
            E = (9*K*G)/(3*K + G)
            # Poisson's ratio
            nu_Phl = (3*K_Phl - 2*G_Phl)/(2*(3*K_Phl + G_Phl))
            nu_Sdp = (3*K_Sdp - 2*G_Sdp)/(2*(3*K_Sdp + G_Sdp))
            nu_Ann = (3*K_Ann - 2*G_Ann)/(2*(3*K_Ann + G_Ann)) # estimated!
            nu_Eas = (3*K_Eas - 2*G_Eas)/(2*(3*K_Eas + G_Eas)) # estimated!
            nu = (3*K - 2*G)/(2*(3*K + G))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP_Phl = ((K_Phl + 4/3*G_Phl)/(rho_Phl))**(0.5)
            vP_Sdp = ((K_Sdp + 4/3*G_Sdp)/(rho_Sdp))**(0.5)
            vP_Ann = ((K_Ann + 4/3*G_Ann)/(rho_Ann))**(0.5)
            vP_Eas = ((K_Eas + 4/3*G_Eas)/(rho_Eas))**(0.5)
            vP = ((K + 4/3*G)/rho)**0.5
            # S-wave velocity
            vS_Phl = (G_Phl/rho_Phl)**0.5
            vS_Sdp = (G_Sdp/rho_Sdp)**0.5
            vS_Ann = (G_Ann/rho_Ann)**0.5
            vS_Eas = (G_Eas/rho_Eas)**0.5
            vS = (G/rho)**0.5
            # Gamma ray
            GR_Phl = (1-x)*potassium[2]/M_Phl*100*16
            GR_Sdp = (1-x)*potassium[2]/M_Sdp*100*16
            GR_Ann = (1-x)*potassium[2]/M_Ann*100*16
            GR_Eas = (1-x)*potassium[2]/M_Eas*100*16
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Phl = 2.30
            #PE_Sdp = 8.66
            #PE_Ann = 11.32
            #PE_Eas = 2.27
            #PE = x*(y*PE_Phl+(1-y)*PE_Eas) + (1-x)*(y*PE_Ann+(1-y)*PE_Sdp)
            #element = [hydrogen, oxygen, flourine, magnesium, aluminium, silicon, potassium, iron]
            w_H = round(2*z*hydrogen[2]/M, 4)
            w_O = round((10+2*z)*oxygen[2]/M, 4)
            w_F = round(2*(1-z)*flourine[2]/M, 4)
            w_Mg = round((3*(1-y)+2*y)*x*magnesium[2]/M, 4)
            w_Al = round((1+y)*aluminium[2]/M, 4)
            w_Si = round((3*(1-y)+2*y)*silicon[2]/M, 4)
            w_K = round(potassium[2]/M, 4)
            w_Fe = round((3*(1-y)+2*y)*(1-x)*iron[2]/M, 4)
            weights = [w_H, w_O, w_F, w_Mg, w_Al, w_Si, w_K, w_Fe]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 500100
            #
            data.append(mineral)
            data.append([round(M, 2), round(x, 2), round(y, 2), round(z, 2)])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(weights)
        elif self.keyword == "Mg":
            x = round(rd.uniform(0.68, 1), 2)
            y = round(rd.uniform(0, 1), 2)
            z = round(rd.uniform(0, 1), 2)
            mineral = "Bt"
            #
            # Molar mass
            M = round(potassium[2] + (3*(1-y)+2*y)*(x*magnesium[2]+(1-x)*iron[2]) + aluminium[2] + (y*aluminium[2]+(3*(1-y)+2*y)*silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            # Density
            M_Phl = round(potassium[2] + 3*magnesium[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Ann = round(potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Sdp = round(potassium[2] + 2*iron[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Eas = round(potassium[2] + 2*magnesium[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            dataV_Phl = CrystalPhysics([[5.3078, 9.1901, 10.1547], [100.08], "monoclinic"])
            V_Phl = dataV_Phl.calculate_volume()
            dataRho_Phl = CrystalPhysics([M_Phl, 2, V_Phl])
            rho_Phl = dataRho_Phl.calculate_bulk_density()
            dataV_Ann = CrystalPhysics([[5.3860, 9.3241, 10.2683], [100.63], "monoclinic"])
            V_Ann = dataV_Ann.calculate_volume()
            dataRho_Ann = CrystalPhysics([M_Ann, 2, V_Ann])
            rho_Ann = dataRho_Ann.calculate_bulk_density()
            dataV_Sdp = CrystalPhysics([[5.348, 9.261, 10.263], [100.19], "monoclinic"])
            V_Sdp = dataV_Sdp.calculate_volume()
            dataRho_Sdp = CrystalPhysics([M_Sdp, 2, V_Sdp])
            rho_Sdp = dataRho_Sdp.calculate_bulk_density()
            dataV_Eas = CrystalPhysics([[5.27, 9.13, 10.15], [99.54], "monoclinic"]) # estimated!
            V_Eas = dataV_Eas.calculate_volume()
            dataRho_Eas = CrystalPhysics([M_Eas, 2, V_Eas])
            rho_Eas = dataRho_Eas.calculate_bulk_density()
            rho = x*(y*rho_Phl+(1-y)*rho_Eas) + (1-x)*(y*rho_Ann+(1-y)*rho_Sdp)
            # Bulk modulus
            K_Phl = 103.78*10**9
            K_Sdp = 110.91*10**9
            K_Ann = 114.72*10**9 # estimated!
            K_Eas = 100.44*10**9 # estimated!
            K = x*(y*K_Phl+(1-y)*K_Eas) + (1-x)*(y*K_Ann+(1-y)*K_Sdp)
            # Shear modulus
            G_Phl = 61.69*10**9
            G_Sdp = 59.61*10**9
            G_Ann = 58.61*10**9
            G_Eas = 62.77*10**9
            G = x*(y*G_Phl+(1-y)*G_Eas) + (1-x)*(y*G_Ann+(1-y)*G_Sdp)
            # Young's modulus
            E_Phl = (9*K_Phl*G_Phl)/(3*K_Phl + G_Phl)
            E_Sdp = (9*K_Sdp*G_Sdp)/(3*K_Sdp + G_Sdp)
            E_Ann = (9*K_Ann*G_Ann)/(3*K_Ann + G_Ann) # estimated!
            E_Eas = (9*K_Eas*G_Eas)/(3*K_Eas + G_Eas) # estimated!
            E = (9*K*G)/(3*K + G)
            # Poisson's ratio
            nu_Phl = (3*K_Phl - 2*G_Phl)/(2*(3*K_Phl + G_Phl))
            nu_Sdp = (3*K_Sdp - 2*G_Sdp)/(2*(3*K_Sdp + G_Sdp))
            nu_Ann = (3*K_Ann - 2*G_Ann)/(2*(3*K_Ann + G_Ann)) # estimated!
            nu_Eas = (3*K_Eas - 2*G_Eas)/(2*(3*K_Eas + G_Eas)) # estimated!
            nu = (3*K - 2*G)/(2*(3*K + G))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP_Phl = ((K_Phl + 4/3*G_Phl)/(rho_Phl))**(0.5)
            vP_Sdp = ((K_Sdp + 4/3*G_Sdp)/(rho_Sdp))**(0.5)
            vP_Ann = ((K_Ann + 4/3*G_Ann)/(rho_Ann))**(0.5)
            vP_Eas = ((K_Eas + 4/3*G_Eas)/(rho_Eas))**(0.5)
            vP = ((K + 4/3*G)/rho)**0.5
            # S-wave velocity
            vS_Phl = (G_Phl/rho_Phl)**0.5
            vS_Sdp = (G_Sdp/rho_Sdp)**0.5
            vS_Ann = (G_Ann/rho_Ann)**0.5
            vS_Eas = (G_Eas/rho_Eas)**0.5
            vS = (G/rho)**0.5
            # Gamma ray
            GR_Phl = (1-x)*potassium[2]/M_Phl*100*16
            GR_Sdp = (1-x)*potassium[2]/M_Sdp*100*16
            GR_Ann = (1-x)*potassium[2]/M_Ann*100*16
            GR_Eas = (1-x)*potassium[2]/M_Eas*100*16
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Phl = 2.30
            #PE_Sdp = 8.66
            #PE_Ann = 11.32
            #PE_Eas = 2.27
            #PE = x*(y*PE_Phl+(1-y)*PE_Eas) + (1-x)*(y*PE_Ann+(1-y)*PE_Sdp)
            #element = [hydrogen, oxygen, flourine, magnesium, aluminium, silicon, potassium, iron]
            w_H = round(2*z*hydrogen[2]/M, 4)
            w_O = round((10+2*z)*oxygen[2]/M, 4)
            w_F = round(2*(1-z)*flourine[2]/M, 4)
            w_Mg = round((3*(1-y)+2*y)*x*magnesium[2]/M, 4)
            w_Al = round((1+y)*aluminium[2]/M, 4)
            w_Si = round((3*(1-y)+2*y)*silicon[2]/M, 4)
            w_K = round(potassium[2]/M, 4)
            w_Fe = round((3*(1-y)+2*y)*(1-x)*iron[2]/M, 4)
            weights = [w_H, w_O, w_F, w_Mg, w_Al, w_Si, w_K, w_Fe]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 500100
            #
            data.append(mineral)
            data.append([round(M, 2), round(x, 2), round(y, 2), round(z, 2)])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(weights)
        elif self.keyword == "Fe":
            x = round(rd.uniform(0, 0.33), 2)
            y = round(rd.uniform(0, 1), 2)
            z = round(rd.uniform(0, 1), 2)
            mineral = "Bt"
            #
            # Molar mass
            M = round(potassium[2] + (3*(1-y)+2*y)*(x*magnesium[2]+(1-x)*iron[2]) + aluminium[2] + (y*aluminium[2]+(3*(1-y)+2*y)*silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            # Density
            M_Phl = round(potassium[2] + 3*magnesium[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Ann = round(potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Sdp = round(potassium[2] + 2*iron[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Eas = round(potassium[2] + 2*magnesium[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            dataV_Phl = CrystalPhysics([[5.3078, 9.1901, 10.1547], [100.08], "monoclinic"])
            V_Phl = dataV_Phl.calculate_volume()
            dataRho_Phl = CrystalPhysics([M_Phl, 2, V_Phl])
            rho_Phl = dataRho_Phl.calculate_bulk_density()
            dataV_Ann = CrystalPhysics([[5.3860, 9.3241, 10.2683], [100.63], "monoclinic"])
            V_Ann = dataV_Ann.calculate_volume()
            dataRho_Ann = CrystalPhysics([M_Ann, 2, V_Ann])
            rho_Ann = dataRho_Ann.calculate_bulk_density()
            dataV_Sdp = CrystalPhysics([[5.348, 9.261, 10.263], [100.19], "monoclinic"])
            V_Sdp = dataV_Sdp.calculate_volume()
            dataRho_Sdp = CrystalPhysics([M_Sdp, 2, V_Sdp])
            rho_Sdp = dataRho_Sdp.calculate_bulk_density()
            dataV_Eas = CrystalPhysics([[5.27, 9.13, 10.15], [99.54], "monoclinic"]) # estimated!
            V_Eas = dataV_Eas.calculate_volume()
            dataRho_Eas = CrystalPhysics([M_Eas, 2, V_Eas])
            rho_Eas = dataRho_Eas.calculate_bulk_density()
            rho = x*(y*rho_Phl+(1-y)*rho_Eas) + (1-x)*(y*rho_Ann+(1-y)*rho_Sdp)
            # Bulk modulus
            K_Phl = 103.78*10**9
            K_Sdp = 110.91*10**9
            K_Ann = 114.72*10**9 # estimated!
            K_Eas = 100.44*10**9 # estimated!
            K = x*(y*K_Phl+(1-y)*K_Eas) + (1-x)*(y*K_Ann+(1-y)*K_Sdp)
            # Shear modulus
            G_Phl = 61.69*10**9
            G_Sdp = 59.61*10**9
            G_Ann = 58.61*10**9
            G_Eas = 62.77*10**9
            G = x*(y*G_Phl+(1-y)*G_Eas) + (1-x)*(y*G_Ann+(1-y)*G_Sdp)
            # Young's modulus
            E_Phl = (9*K_Phl*G_Phl)/(3*K_Phl + G_Phl)
            E_Sdp = (9*K_Sdp*G_Sdp)/(3*K_Sdp + G_Sdp)
            E_Ann = (9*K_Ann*G_Ann)/(3*K_Ann + G_Ann) # estimated!
            E_Eas = (9*K_Eas*G_Eas)/(3*K_Eas + G_Eas) # estimated!
            E = (9*K*G)/(3*K + G)
            # Poisson's ratio
            nu_Phl = (3*K_Phl - 2*G_Phl)/(2*(3*K_Phl + G_Phl))
            nu_Sdp = (3*K_Sdp - 2*G_Sdp)/(2*(3*K_Sdp + G_Sdp))
            nu_Ann = (3*K_Ann - 2*G_Ann)/(2*(3*K_Ann + G_Ann)) # estimated!
            nu_Eas = (3*K_Eas - 2*G_Eas)/(2*(3*K_Eas + G_Eas)) # estimated!
            nu = (3*K - 2*G)/(2*(3*K + G))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP_Phl = ((K_Phl + 4/3*G_Phl)/(rho_Phl))**(0.5)
            vP_Sdp = ((K_Sdp + 4/3*G_Sdp)/(rho_Sdp))**(0.5)
            vP_Ann = ((K_Ann + 4/3*G_Ann)/(rho_Ann))**(0.5)
            vP_Eas = ((K_Eas + 4/3*G_Eas)/(rho_Eas))**(0.5)
            vP = ((K + 4/3*G)/rho)**0.5
            # S-wave velocity
            vS_Phl = (G_Phl/rho_Phl)**0.5
            vS_Sdp = (G_Sdp/rho_Sdp)**0.5
            vS_Ann = (G_Ann/rho_Ann)**0.5
            vS_Eas = (G_Eas/rho_Eas)**0.5
            vS = (G/rho)**0.5
            # Gamma ray
            GR_Phl = (1-x)*potassium[2]/M_Phl*100*16
            GR_Sdp = (1-x)*potassium[2]/M_Sdp*100*16
            GR_Ann = (1-x)*potassium[2]/M_Ann*100*16
            GR_Eas = (1-x)*potassium[2]/M_Eas*100*16
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Phl = 2.30
            #PE_Sdp = 8.66
            #PE_Ann = 11.32
            #PE_Eas = 2.27
            #PE = x*(y*PE_Phl+(1-y)*PE_Eas) + (1-x)*(y*PE_Ann+(1-y)*PE_Sdp)
            #element = [hydrogen, oxygen, flourine, magnesium, aluminium, silicon, potassium, iron]
            w_H = round(2*z*hydrogen[2]/M, 4)
            w_O = round((10+2*z)*oxygen[2]/M, 4)
            w_F = round(2*(1-z)*flourine[2]/M, 4)
            w_Mg = round((3*(1-y)+2*y)*x*magnesium[2]/M, 4)
            w_Al = round((1+y)*aluminium[2]/M, 4)
            w_Si = round((3*(1-y)+2*y)*silicon[2]/M, 4)
            w_K = round(potassium[2]/M, 4)
            w_Fe = round((3*(1-y)+2*y)*(1-x)*iron[2]/M, 4)
            weights = [w_H, w_O, w_F, w_Mg, w_Al, w_Si, w_K, w_Fe]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 500100
            #
            data.append(mineral)
            data.append([round(M, 2), round(x, 2), round(y, 2), round(z, 2)])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(weights)
        elif self.keyword == "Biotite":
            x = round(rd.uniform(0.9, 1), 2)
            y = round(rd.uniform(0, 0.25), 2)
            z = round(rd.uniform(0, 0.5), 2)
            mineral = "Bt"
            #
            # Molar mass
            M = round(potassium[2] + (3*(1-y)+2*y)*(x*magnesium[2]+(1-x)*iron[2]) + aluminium[2] + (y*aluminium[2]+(3*(1-y)+2*y)*silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            # Density
            M_Phl = round(potassium[2] + 3*magnesium[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Ann = round(potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Sdp = round(potassium[2] + 2*iron[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            M_Eas = round(potassium[2] + 2*magnesium[2] + aluminium[2] + 2*(aluminium[2]+silicon[2]) + 10*oxygen[2] + 2*(z*(oxygen[2]+hydrogen[2])+(1-z)*flourine[2]), 3)
            dataV_Phl = CrystalPhysics([[5.3078, 9.1901, 10.1547], [100.08], "monoclinic"])
            V_Phl = dataV_Phl.calculate_volume()
            dataRho_Phl = CrystalPhysics([M_Phl, 2, V_Phl])
            rho_Phl = dataRho_Phl.calculate_bulk_density()
            dataV_Ann = CrystalPhysics([[5.3860, 9.3241, 10.2683], [100.63], "monoclinic"])
            V_Ann = dataV_Ann.calculate_volume()
            dataRho_Ann = CrystalPhysics([M_Ann, 2, V_Ann])
            rho_Ann = dataRho_Ann.calculate_bulk_density()
            dataV_Sdp = CrystalPhysics([[5.348, 9.261, 10.263], [100.19], "monoclinic"])
            V_Sdp = dataV_Sdp.calculate_volume()
            dataRho_Sdp = CrystalPhysics([M_Sdp, 2, V_Sdp])
            rho_Sdp = dataRho_Sdp.calculate_bulk_density()
            dataV_Eas = CrystalPhysics([[5.27, 9.13, 10.15], [99.54], "monoclinic"]) # estimated!
            V_Eas = dataV_Eas.calculate_volume()
            dataRho_Eas = CrystalPhysics([M_Eas, 2, V_Eas])
            rho_Eas = dataRho_Eas.calculate_bulk_density()
            rho = x*(y*rho_Phl+(1-y)*rho_Eas) + (1-x)*(y*rho_Ann+(1-y)*rho_Sdp)
            # Bulk modulus
            K_Phl = 103.78*10**9
            K_Sdp = 110.91*10**9
            K_Ann = 114.72*10**9 # estimated!
            K_Eas = 100.44*10**9 # estimated!
            K = x*(y*K_Phl+(1-y)*K_Eas) + (1-x)*(y*K_Ann+(1-y)*K_Sdp)
            # Shear modulus
            G_Phl = 61.69*10**9
            G_Sdp = 59.61*10**9
            G_Ann = 58.61*10**9
            G_Eas = 62.77*10**9
            G = x*(y*G_Phl+(1-y)*G_Eas) + (1-x)*(y*G_Ann+(1-y)*G_Sdp)
            # Young's modulus
            E_Phl = (9*K_Phl*G_Phl)/(3*K_Phl + G_Phl)
            E_Sdp = (9*K_Sdp*G_Sdp)/(3*K_Sdp + G_Sdp)
            E_Ann = (9*K_Ann*G_Ann)/(3*K_Ann + G_Ann) # estimated!
            E_Eas = (9*K_Eas*G_Eas)/(3*K_Eas + G_Eas) # estimated!
            E = (9*K*G)/(3*K + G)
            # Poisson's ratio
            nu_Phl = (3*K_Phl - 2*G_Phl)/(2*(3*K_Phl + G_Phl))
            nu_Sdp = (3*K_Sdp - 2*G_Sdp)/(2*(3*K_Sdp + G_Sdp))
            nu_Ann = (3*K_Ann - 2*G_Ann)/(2*(3*K_Ann + G_Ann)) # estimated!
            nu_Eas = (3*K_Eas - 2*G_Eas)/(2*(3*K_Eas + G_Eas)) # estimated!
            nu = (3*K - 2*G)/(2*(3*K + G))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP_Phl = ((K_Phl + 4/3*G_Phl)/(rho_Phl))**(0.5)
            vP_Sdp = ((K_Sdp + 4/3*G_Sdp)/(rho_Sdp))**(0.5)
            vP_Ann = ((K_Ann + 4/3*G_Ann)/(rho_Ann))**(0.5)
            vP_Eas = ((K_Eas + 4/3*G_Eas)/(rho_Eas))**(0.5)
            vP = ((K + 4/3*G)/rho)**0.5
            # S-wave velocity
            vS_Phl = (G_Phl/rho_Phl)**0.5
            vS_Sdp = (G_Sdp/rho_Sdp)**0.5
            vS_Ann = (G_Ann/rho_Ann)**0.5
            vS_Eas = (G_Eas/rho_Eas)**0.5
            vS = (G/rho)**0.5
            # Gamma ray
            GR_Phl = (1-x)*potassium[2]/M_Phl*100*16
            GR_Sdp = (1-x)*potassium[2]/M_Sdp*100*16
            GR_Ann = (1-x)*potassium[2]/M_Ann*100*16
            GR_Eas = (1-x)*potassium[2]/M_Eas*100*16
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Phl = 2.30
            #PE_Sdp = 8.66
            #PE_Ann = 11.32
            #PE_Eas = 2.27
            #PE = x*(y*PE_Phl+(1-y)*PE_Eas) + (1-x)*(y*PE_Ann+(1-y)*PE_Sdp)
            #element = [hydrogen, oxygen, flourine, magnesium, aluminium, silicon, potassium, iron]
            w_H = round(2*z*hydrogen[2]/M, 4)
            w_O = round((10+2*z)*oxygen[2]/M, 4)
            w_F = round(2*(1-z)*flourine[2]/M, 4)
            w_Mg = round((3*(1-y)+2*y)*x*magnesium[2]/M, 4)
            w_Al = round((1+y)*aluminium[2]/M, 4)
            w_Si = round((3*(1-y)+2*y)*silicon[2]/M, 4)
            w_K = round(potassium[2]/M, 4)
            w_Fe = round((3*(1-y)+2*y)*(1-x)*iron[2]/M, 4)
            weights = [w_H, w_O, w_F, w_Mg, w_Al, w_Si, w_K, w_Fe]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 500100
            #
            data.append(mineral)
            data.append([round(M, 2), round(x, 2), round(y, 2), round(z, 2)])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(weights)
        #
        elif self.keyword == "simple":
            x = round(rd.uniform(0.9, 1), 2)
            y = round(rd.uniform(0, 0.25), 2)
            mineral = "Bt"
            #
            # Molar mass
            M = round(potassium[2] + 3*(x*magnesium[2]+(1-x)*iron[2]) + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(y*(oxygen[2]+hydrogen[2])+(1-y)*flourine[2]), 3)
            # Density
            M_Phl = round(potassium[2] + 3*magnesium[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(y*(oxygen[2]+hydrogen[2])+(1-y)*flourine[2]), 3)
            M_Ann = round(potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(y*(oxygen[2]+hydrogen[2])+(1-y)*flourine[2]), 3)
            dataV_Phl = CrystalPhysics([[5.3078, 9.1901, 10.1547], [100.08], "monoclinic"])
            V_Phl = dataV_Phl.calculate_volume()
            dataRho_Phl = CrystalPhysics([M_Phl, 2, V_Phl])
            rho_Phl = dataRho_Phl.calculate_bulk_density()
            dataV_Ann = CrystalPhysics([[5.3860, 9.3241, 10.2683], [100.63], "monoclinic"])
            V_Ann = dataV_Ann.calculate_volume()
            dataRho_Ann = CrystalPhysics([M_Ann, 2, V_Ann])
            rho_Ann = dataRho_Ann.calculate_bulk_density()
            rho = x*((1-y)*rho_Phl) + (1-x)*(y*rho_Ann)
            # Bulk modulus
            K_Phl = 103.78*10**9
            K_Ann = 114.72*10**9 # estimated!
            K = x*((1-y)*K_Phl) + (1-x)*(y*K_Ann)
            # Shear modulus
            G_Phl = 61.69*10**9
            G_Ann = 58.61*10**9
            G = x*((1-y)*G_Phl) + (1-x)*(y*G_Ann)
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #element = [hydrogen, oxygen, flourine, magnesium, aluminium, silicon, potassium, iron]
            w_H = round(2*y*hydrogen[2]/M, 4)
            w_O = round((10+2*y)*oxygen[2]/M, 4)
            w_F = round(2*(1-y)*flourine[2]/M, 4)
            w_Mg = round(3*x*magnesium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_K = round(potassium[2]/M, 4)
            w_Fe = round(3*(1-x)*iron[2]/M, 4)
            weights = [w_H, w_O, w_F, w_Mg, w_Al, w_Si, w_K, w_Fe]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 500100
            #
            data.append(mineral)
            data.append([round(M, 2), round(x, 2), round(y, 2)])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(weights)
        #
        return data
#
#################
# NESOSILICATES #
#################
#
class nesosilicates:
    #
    def __init__(self):
        pass
    #
    def olivine(self, keyword="olivine"):  # (Mg,Fe,Mn)2SiO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        iron = elements.Fe(self)
        manganese = elements.Mn(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        element = [oxygen, magnesium, silicon, manganese, iron]
        self.keyword = keyword
        #
        if self.keyword == "Olivine" or self.keyword == "olivine":
            x_Mg = round(rd.uniform(0.5, 1), 2)
            x_Fe = round(rd.uniform(0, float(1-x_Mg)), 2)
            x_Mn = round(1-x_Mg-x_Fe, 2)
            mineral = "Ol"
        elif self.keyword == "Mg":
            x_Mg = round(rd.uniform(0.68, 1), 2)
            x_Fe = round(rd.uniform(0, float(1 - x_Mg)), 2)
            x_Mn = round(1-x_Mg-x_Fe, 2)
            mineral = "Ol"
        elif self.keyword == "Fe":
            x_Fe = round(rd.uniform(0.68, 1), 2)
            x_Mg = round(rd.uniform(0, float(1 - x_Fe)), 2)
            x_Mn = round(1-x_Mg-x_Fe, 2)
            mineral = "Ol"
        elif self.keyword == "Mn":
            x_Mn = round(rd.uniform(0.68, 1), 2)
            x_Mg = round(rd.uniform(0, float(1 - x_Mn)), 2)
            x_Fe = round(1-x_Mn-x_Mg, 2)
            mineral = "Ol"
        elif self.keyword in ["Fo", "Forsterite", "forsterite"]:
            x_Mg = 1.0
            x_Fe = 0.0
            x_Mn = 0.0
            mineral = "Fo"
        elif self.keyword in ["Fa", "Fayalite", "fayalite"]:
            x_Mg = 0.0
            x_Fe = 1.0
            x_Mn = 0.0
            mineral = "Fa"
        elif self.keyword in ["Tep", "Tephroite", "tephroite"]:
            x_Mg = 0.0
            x_Fe = 0.0
            x_Mn = 1.0
            mineral = "Tep"
        elif self.keyword == "None":
            x = round(rd.uniform(0, 1), 2)
            if 0 <= x < 0.5:
                x_Mg = round(rd.uniform(0.5, 1), 2)
                x_Fe = round(rd.uniform(0, float(1 - x_Mg)), 2)
                x_Mn = round(1-x_Mg-x_Fe, 2)
                mineral = "Fo"
            elif 0.5 <= x < 0.8:
                x_Fe = round(rd.uniform(0.5, 1), 2)
                x_Mg = round(rd.uniform(0, float(1 - x_Fe)), 2)
                x_Mn = round(1-x_Mg-x_Fe, 2)
                mineral = "Fa"
            elif 0.8 <= x <= 1.0:
                x_Mn = round(rd.uniform(0.5, 1), 2)
                x_Mg = round(rd.uniform(0, float(1 - x_Mn)), 2)
                x_Fe = round(1-x_Mn-x_Mg, 2)
                mineral = "Tep"
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], [GR, PE]]
        data = []
        #
        # Molar mass
        M = round(2*(x_Mg*magnesium[2]+x_Fe*iron[2]+x_Mn*manganese[2]) + silicon[2] + 4*oxygen[2], 3)
        # Density
        M_Fo = round(2*magnesium[2] + silicon[2] + 4*oxygen[2], 3)
        dataV_Fo = CrystalPhysics([[4.754, 10.1971, 5.9806], [], "orthorhombic"])   # source: https://www.mindat.org/min-1584.html
        V_Fo = dataV_Fo.calculate_volume()
        dataRho_Fo = CrystalPhysics([M_Fo, 4, V_Fo])
        rho_Fo = dataRho_Fo.calculate_bulk_density()
        M_Fa = round(2*iron[2] + silicon[2] + 4*oxygen[2], 3)
        dataV_Fa = CrystalPhysics([[4.79, 10.39, 6.06], [], "orthorhombic"])    # source: https://www.mindat.org/min-1458.html
        V_Fa = dataV_Fa.calculate_volume()
        dataRho_Fa = CrystalPhysics([M_Fa, 4, V_Fa])
        rho_Fa = dataRho_Fa.calculate_bulk_density()
        M_Tep = round(2*manganese[2] + silicon[2] + 4*oxygen[2], 3)
        dataV_Tep = CrystalPhysics([[4.88, 10.61, 6.24], [], "orthorhombic"])   # source: https://www.mindat.org/min-3913.html
        V_Tep = dataV_Tep.calculate_volume()
        dataRho_Tep = CrystalPhysics([M_Tep, 4, V_Tep])
        rho_Tep = dataRho_Tep.calculate_bulk_density()
        rho = x_Mg*rho_Fo + x_Fe*rho_Fa + x_Mn*rho_Tep
        # Bulk modulus
        K_Fo = 119*10**9
        K_Fa = 142.52*10**9
        K_Tep = 116*10**9
        K = x_Mg*K_Fo + x_Fe*K_Fa + x_Mn*K_Tep
        # Shear modulus
        G_Fo = 74*10**9
        G_Fa = 64.77*10**9
        G_Tep = 50*10**9
        G = x_Mg*G_Fo + x_Fe*G_Fa + x_Mn*G_Tep
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(6*K + 2*G)
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        GR = 0
        # Photoelectricity
        #PE_Fo = 1.53
        #PE_Fa = 17.09
        #PE_Tep = 14.67
        #PE = x_Mg*PE_Fo + x_Fe*PE_Fa + x_Mn*PE_Tep
        #element = [magnesium, iron, manganese, silicon, oxygen]
        w_Mg = round(2*x_Mg*magnesium[2]/M, 4)
        w_Fe = round(2*x_Fe*iron[2]/M, 4)
        w_Mn = round(2*x_Mn*manganese[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_O, w_Mg, w_Si, w_Mn, w_Fe]

        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x_Mg, x_Fe, x_Mn])
        data.append(round(rho, 1))
        data.append([round(K * 10 ** (-9), 2), round(G * 10 ** (-9), 2), round(E * 10 ** (-9), 2), round(nu, 4), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(weights)
        #
        return data
    #
    def forsterite(self):  # Mg2SiO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Fo"
        #
        # Molar mass
        M = round(2*magnesium[2] + silicon[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[4.8, 10.35, 6.06], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 119*10**9
        # Shear modulus
        G = 74*10**9
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
        element = [magnesium, silicon, oxygen]
        w_Mg = round(2*magnesium[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_Mg, w_Si, w_O]
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
    def fayalite(self):  # Fe2SiO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Fa"
        #
        # Molar mass
        M = round(2*iron[2] + silicon[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[4.82, 10.48, 6.09], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 142.52*10**9
        # Shear modulus
        G = 64.77*10**9
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
        element = [iron, silicon, oxygen]
        w_Fe = round(2*iron[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_Fe, w_Si, w_O]
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
    def tephroite(self):  # Mn2SiO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        manganese = elements.Mn(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Tep"
        #
        # Molar mass
        M = round(2*manganese[2] + silicon[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[4.9, 10.6, 6.26], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 116*10**9
        # Shear modulus
        G = 50*10**9
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
        element = [manganese, silicon, oxygen]
        w_Mn = round(2*manganese[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_Mn, w_Si, w_O]
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
    def larnite(self):  # Ca2SiO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Lrn"
        #
        # Molar mass
        M = round(2*calcium[2] + silicon[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[5.5, 6.74, 9.3], [94.6], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 96*10**9
        # Shear modulus
        G = 48*10**9
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
        element = [calcium, silicon, oxygen]
        w_Ca = round(2*calcium[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_Ca, w_Si, w_O]
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
    def garnet_pyralspite(self, keyword=None):  # (Mg,Fe,Mn)3 Al2 Si3 O12
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        manganese = elements.Mn(self)
        iron = elements.Fe(self)
        element = [oxygen, magnesium, aluminium, silicon, manganese, iron]
        self.keyword = keyword
        #
        if self.keyword == None:
            x_Mg = round(rd.uniform(0, 1), 2)
            x_Fe = round(rd.uniform(0, float(1-x_Mg)), 2)
            x_Mn = round(1-x_Mg-x_Fe, 2)
        elif self.keyword == "Mg":
            x_Mg = round(rd.uniform(0.5, 1), 2)
            x_Fe = round(rd.uniform(0, float(1-x_Mg)), 2)
            x_Mn = round(1-x_Mg-x_Fe, 2)
        elif self.keyword == "Fe":
            x_Fe = round(rd.uniform(0.5, 1), 2)
            x_Mg = round(rd.uniform(0, float(1-x_Fe)), 2)
            x_Mn = round(1-x_Mg-x_Fe, 2)
        elif self.keyword == "Mn":
            x_Mn = round(rd.uniform(0.5, 1), 2)
            x_Mg = round(rd.uniform(0, float(1-x_Mn)), 2)
            x_Fe = round(1-x_Mg-x_Mn, 2)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], [GR, PE]]
        data = []
        #
        mineral = "Grt"
        #
        # Molar mass
        M = round(3*(x_Mg*magnesium[2]+x_Fe*iron[2]+x_Mn*manganese[2]) + 2*aluminium[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        # Density
        M_Alm = round(3*iron[2] + 2*aluminium[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        dataV_Alm = CrystalPhysics([[11.526], [], "cubic"])
        V_Alm = dataV_Alm.calculate_volume()
        dataRho_Alm = CrystalPhysics([M_Alm, 8, V_Alm])
        rho_Alm = dataRho_Alm.calculate_bulk_density()
        M_Prp = round(3*magnesium[2] + 2*aluminium[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        dataV_Prp = CrystalPhysics([[11.459], [], "cubic"])
        V_Prp = dataV_Prp.calculate_volume()
        dataRho_Prp = CrystalPhysics([M_Prp, 8, V_Prp])
        rho_Prp = dataRho_Prp.calculate_bulk_density()
        M_Sps = round(3*manganese[2] + 2*aluminium[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        dataV_Sps = CrystalPhysics([[11.621], [], "cubic"])
        V_Sps = dataV_Sps.calculate_volume()
        dataRho_Sps = CrystalPhysics([M_Sps, 8, V_Sps])
        rho_Sps = dataRho_Sps.calculate_bulk_density()
        rho = x_Fe*rho_Alm + x_Mg*rho_Prp + x_Mn*rho_Sps
        # Bulk modulus
        K_Alm = 168.52*10**9
        K_Prp = 160*10**9
        K_Sps = 175.88*10**9
        K = x_Fe*K_Alm + x_Mg*K_Prp + x_Mn*K_Sps
        # Shear modulus
        G_Alm = 83.52*10**9
        G_Prp = 85*10**9
        G_Sps = 93.58*10**9
        G = x_Fe*G_Alm + x_Mg*G_Prp + x_Mn*G_Sps
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
        #element = [oxygen, magnesium, aluminium, silicon, manganese, iron]
        w_O = round(12*oxygen[2]/M, 4)
        w_Mg = round(3*x_Mg*magnesium[2]/M, 4)
        w_Al = round(2*aluminium[2]/M, 4)
        w_Si = round(3*silicon[2]/M, 4)
        w_Mn = round(3*x_Mn*manganese[2]/M, 4)
        w_Fe = round(3*x_Fe*iron[2]/M, 4)
        weights = [w_O, w_Mg, w_Al, w_Si, w_Mn, w_Fe]
        amounts = [12, 3*x_Mg, 2, 3, 3*x_Mn, 3*x_Fe]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x_Mg, x_Fe, x_Mn])
        data.append(round(rho, 1))
        data.append([round(K * 10 ** (-9), 2), round(G * 10 ** (-9), 2), round(E * 10 ** (-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
    #
    def garnet_ugrandite(self, keyword=None):  # Ca3 (Cr,Al,Fe)2 Si3 O12
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        oxygen = elements.O(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        calcium = elements.Ca(self)
        chromium = elements.Cr(self)
        iron = elements.Fe(self)
        element = [oxygen, aluminium, silicon, calcium, chromium, iron]
        self.keyword = keyword
        #
        if self.keyword == None:
            x_Al = round(rd.uniform(0, 1), 2)
            x_Cr = round(rd.uniform(0, float(1-x_Al)), 2)
            x_Fe = round(1-x_Al-x_Cr, 2)
        elif self.keyword == "Al":
            x_Al = round(rd.uniform(0.5, 1), 2)
            x_Cr = round(rd.uniform(0, float(1-x_Al)), 2)
            x_Fe = round(1-x_Al-x_Cr, 2)
        elif self.keyword == "Cr":
            x_Cr = round(rd.uniform(0.5, 1), 2)
            x_Al = round(rd.uniform(0, float(1-x_Cr)), 2)
            x_Fe = round(1-x_Al-x_Cr, 2)
        elif self.keyword == "Fe":
            x_Fe = round(rd.uniform(0.5, 1), 2)
            x_Al = round(rd.uniform(0, float(1-x_Fe)), 2)
            x_Cr = round(1-x_Fe-x_Al, 2)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], [GR, PE]]
        data = []
        #
        mineral = "Grt"
        #
        # Molar mass
        M = round(3*calcium[2] + 2*(x_Al*aluminium[2]+x_Cr*chromium[2]+x_Fe*iron[2]) + 3*(silicon[2]+4*oxygen[2]), 3)
        # Density
        M_Grs = round(3*calcium[2] + 2*aluminium[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        dataV_Grs = CrystalPhysics([[11.851], [], "cubic"])
        V_Grs = dataV_Grs.calculate_volume()
        dataRho_Grs = CrystalPhysics([M_Grs, 8, V_Grs])
        rho_Grs = dataRho_Grs.calculate_bulk_density()
        M_Uv = round(3*calcium[2] + 2*chromium[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        dataV_Uv = CrystalPhysics([[12.0], [], "cubic"])
        V_Uv = dataV_Uv.calculate_volume()
        dataRho_Uv = CrystalPhysics([M_Uv, 8, V_Uv])
        rho_Uv = dataRho_Uv.calculate_bulk_density()
        M_Adr = round(3*calcium[2] + 2*iron[2] + 3*(silicon[2]+4*oxygen[2]), 3)
        dataV_Adr = CrystalPhysics([[12.05], [], "cubic"])
        V_Adr = dataV_Adr.calculate_volume()
        dataRho_Adr = CrystalPhysics([M_Adr, 8, V_Adr])
        rho_Adr = dataRho_Adr.calculate_bulk_density()
        rho = x_Al*rho_Grs + x_Cr*rho_Uv + x_Fe*rho_Adr
        # Bulk modulus
        K_Grs = 154*10**9
        K_Uv = 136.17*10**9
        K_Adr = 128.57*10**9
        K = x_Al*K_Grs + x_Cr*K_Uv + x_Fe*K_Adr
        # Shear modulus
        G_Grs = 97*10**9
        G_Uv = 82.58*10**9
        G_Adr = 72.90*10**9
        G = x_Al*G_Grs + x_Cr*G_Uv + x_Fe*G_Adr
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
        #element = [oxygen, aluminium, silicon, calcium, chromium, iron]
        w_O = round(12*oxygen[2]/M, 4)
        w_Al = round(2*x_Al*aluminium[2]/M, 4)
        w_Si = round(3*silicon[2]/M, 4)
        w_Ca = round(3*calcium[2]/M, 4)
        w_Cr = round(2*x_Cr*chromium[2]/M, 4)
        w_Fe = round(2*x_Fe*iron[2]/M, 4)
        weights = [w_O, w_Al, w_Si, w_Ca, w_Cr, w_Fe]
        amounts = [12, 2*x_Al, 3, 3, 2*x_Cr, 2*x_Fe]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x_Al, x_Cr, x_Fe])
        data.append(round(rho, 1))
        data.append([round(K * 10 ** (-9), 2), round(G * 10 ** (-9), 2), round(E * 10 ** (-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
    #
    def topaz(self):  # Al2SiO4(F,OH)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        fluorine = elements.F(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Tpz"
        #
        # Molar mass
        x = round(rd.uniform(0.5, 0.6), 2)
        M = round(2*aluminium[2] + silicon[2] + 4*oxygen[2] + 2*(x*fluorine[2] + (1-x)*(oxygen[2]+hydrogen[2])), 3)
        # Density
        dataV = CrystalPhysics([[4.35, 8.8, 8.4], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 166*10**9
        # Shear modulus
        G = 114*10**9
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
        element = [hydrogen, oxygen, fluorine, aluminium, silicon]
        w_H = round(2*(1-x)*hydrogen[2]/M, 4)
        w_O = round((4*oxygen[2] + 2*(1-x)*oxygen[2])/M, 4)
        w_F = round(2*x*fluorine[2]/M, 4)
        w_Al = round(2*aluminium[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        weights = [w_H, w_O, w_F, w_Al, w_Si]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(weights)
        #
        return data
    #
################
# INOSILICATES #
################
#
class inosilicates:
    #
    def __init__(self):
        pass
    #
    def enstatite(self):  # Mg2Si2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "En"
        #
        # Molar mass
        M = round(2*magnesium[2] + 2*silicon[2] + 6*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[18.228, 8.805, 5.185], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 138.65*10**9
        # Shear modulus
        G = 82.51*10**9
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
        element = [magnesium, silicon, oxygen]
        w_Mg = round(2*magnesium[2]/M, 4)
        w_Si = round(2*silicon[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        weights = [w_Mg, w_Si, w_O]
        composition = [w_O, w_Mg, w_Si]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def wollastonite(self):  # CaSiO3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Wo"
        #
        # Molar mass
        M = round(calcium[2] + silicon[2] + 3*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[7.94, 7.32, 7.07], [90.033, 95.367, 103.433], "triclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 93*10**9
        # Shear modulus
        G = 47*10**9
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
        element = [calcium, silicon, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(3*oxygen[2]/M, 4)
        weights = [w_Ca, w_Si, w_O]
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
    def ferrosilite(self):  # FeSiO3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Fs"
        #
        # Molar mass
        M = round(iron[2] + silicon[2] + 3*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[18.418, 9.078, 5.237], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 141.76*10**9
        # Shear modulus
        G = 70.49*10**9
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
        element = [iron, silicon, oxygen]
        w_Fe = round(iron[2]/M, 4)
        w_Si = round(silicon[2]/M, 4)
        w_O = round(3*oxygen[2]/M, 4)
        weights = [w_Fe, w_Si, w_O]
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
    def tremolite(self):  # Ca2Mg5Si8O22(OH)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Tr"
        #
        # Molar mass
        M = round(2*calcium[2] + 5*magnesium[2] + 2*(4*silicon[2]+11*oxygen[2]+(oxygen[2]+hydrogen[2])), 3)
        # Density
        dataV = CrystalPhysics([[9.86, 18.05, 5.29], [104.8], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho_Or = CrystalPhysics([M, 2, V])
        rho = dataRho_Or.calculate_bulk_density()
        # Bulk modulus
        K = 124.07*10**9
        # Shear modulus
        G = 72.88*10**9
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
        element = [calcium, magnesium, silicon, oxygen, hydrogen]
        w_Ca = round(2*calcium[2]/M, 4)
        w_Mg = round(5*magnesium[2]/M, 4)
        w_Si = round(2*4*silicon[2]/M, 4)
        w_O = round(2*(11+1)*oxygen[2]/M, 4)
        w_H = round(2*hydrogen[2]/M, 4)
        weights = [w_Ca, w_Mg, w_Si, w_O, w_H]
        composition = [w_H, w_O, w_Mg, w_Si, w_Ca]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def actinolite(self):   # Ca2(Mg,Fe)5Si8O22(OH)2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        magnesium = elements.Mg(self)
        iron = elements.Fe(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Act"
        #
        # Molar mass
        x = round(rd.uniform(0, 1), 2)
        M = round(2*calcium[2] + 5*(x*magnesium[2]+(1-x)*iron[2]) + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2]+hydrogen[2]), 3)
        # Density
        dataV = CrystalPhysics([[9.891, 18.200, 5.305], [104.64], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho_Or = CrystalPhysics([M, 2, V])
        rho = dataRho_Or.calculate_bulk_density()
        # Bulk modulus
        K = 85*10**9
        # Shear modulus
        G = 47*10**9
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
        element = [calcium, magnesium, iron, silicon, oxygen, hydrogen]
        w_Ca = round(2*calcium[2]/M, 4)
        w_Mg = round(5*x*magnesium[2]/M, 4)
        w_Fe = round(5*(1-x)*iron[2]/M, 4)
        w_Si = round(8*silicon[2]/M, 4)
        w_O = round((22+2)*oxygen[2]/M, 4)
        w_H = round(2*hydrogen[2]/M, 4)
        weights = [w_Ca, w_Mg, w_Fe, w_Si, w_O, w_H]
        composition = [w_H, w_O, w_Mg, w_Si, w_Ca, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def augite(self):   # (Ca,Mg,Fe)(Mg,Fe)Si2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        magnesium = elements.Mg(self)
        iron = elements.Fe(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Aug"
        #
        # Molar mass
        x_Ca = round(rd.uniform(0.4, 0.9), 2)
        x_Mg1 = round(rd.uniform(0.0, 1-x_Ca), 2)
        x_Fe1 = 1-x_Ca-x_Mg1
        x_Mg2 = round(rd.uniform(0.0, 1.0), 2)
        x_Fe2 = 1-x_Mg2
        M = round((x_Ca*calcium[2]+x_Mg1*magnesium[2]+x_Fe1*iron[2]) + (x_Mg2*magnesium[2]+x_Fe2*iron[2]) + 2*silicon[2] + 6*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[9.699, 8.844, 5.272], [106.97], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho_Or = CrystalPhysics([M, 4, V])
        rho = dataRho_Or.calculate_bulk_density()
        # Bulk modulus
        K = (95.72 + (106.97-95.72)/(3420-3320)*(rho-3320))*10**9
        # Shear modulus
        G = (58.01 + (57.21-58.01)/(3420-3320)*(rho-3320))*10**9
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
        element = [calcium, magnesium, iron, silicon, oxygen]
        w_Ca = round(x_Ca*calcium[2]/M, 4)
        w_Mg = round((x_Mg1+x_Mg2)*magnesium[2]/M, 4)
        w_Fe = round((x_Fe1+x_Fe2)*iron[2]/M, 4)
        w_Si = round(2*silicon[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        weights = [w_Ca, w_Mg, w_Fe, w_Si, w_O]
        composition = [w_O, w_Mg, w_Si, w_Ca, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x_Ca, x_Mg1, x_Mg2])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def diopside(self):  # CaMgSi2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Di"
        #
        # Molar mass
        M = round(calcium[2] + magnesium[2] + 2*silicon[2] + 6*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[9.761, 8.926, 5.258], [105.8], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 102*10**9
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
        element = [calcium, magnesium, silicon, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_Mg = round(magnesium[2]/M, 4)
        w_Si = round(2*silicon[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        weights = [w_Ca, w_Mg, w_Si, w_O]
        composition = [w_O, w_Mg, w_Si, w_Ca]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def aegirine(self):  # NaFeSi2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        iron = elements.Fe(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Aeg"
        #
        # Molar mass
        M = round(sodium[2] + iron[2] + 2*silicon[2] + 6*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[9.65, 8.79, 5.29], [107.5], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 124.44*10**9
        # Shear modulus
        G = 74.16*10**9
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
        element = [sodium, iron, silicon, oxygen]
        w_Na = round(sodium[2]/M, 4)
        w_Fe = round(iron[2]/M, 4)
        w_Si = round(2*silicon[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        weights = [w_Na, w_Fe, w_Si, w_O]
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
    def hedenbergite(self):  # CaFeSi2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        iron = elements.Fe(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Hd"
        #
        # Molar mass
        M = round(calcium[2] + iron[2] + 2*silicon[2] + 6*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[9.827, 8.994, 5.261], [105.52], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 134.38*10**9
        # Shear modulus
        G = 77.75*10**9
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
        element = [calcium, iron, silicon, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_Fe = round(iron[2]/M, 4)
        w_Si = round(2*silicon[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        weights = [w_Ca, w_Fe, w_Si, w_O]
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
    def riebeckite(self):   # Na2Fe5(Si8O22)(OH)2
        # Chemistry
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        sodium = elements.Na(self)
        silicon = elements.Si(self)
        iron = elements.Fe(self)
        #
        data = []
        #
        mineral = "Rbk"
        #
        # Molar mass
        M = round(2*sodium[2] + 5*iron[2] + (8*silicon[2] + 22*oxygen[2]) + 2*(oxygen[2] + hydrogen[2]), 3)
        # Density
        data_V = CrystalPhysics([[9.769, 18.048, 5.335], [103.6], "monoclinic"])
        V = data_V.calculate_volume()
        Z = 2
        data_rho = CrystalPhysics([M, Z, V])
        rho = data_rho.calculate_bulk_density()
        # Bulk modulus
        K = 101.1*10**9
        # Shear modulus
        G = 43.7*10**9
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
        element = [hydrogen, oxygen, sodium, silicon, iron]
        w_H = round(2*hydrogen[2]/M, 4)
        w_O = round((22+2)*oxygen[2]/M, 4)
        w_Na = round(2*sodium[2]/M, 4)
        w_Si = round(8*silicon[2]/M, 4)
        w_Fe = round(5*iron[2]/M, 4)
        weights = [w_H, w_O, w_Na, w_Si, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(weights)
        #
        return data
    #
    def arfvedsonite(self):   # Na3Fe5(Si8O22)(OH)2
        # Chemistry
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        sodium = elements.Na(self)
        silicon = elements.Si(self)
        iron = elements.Fe(self)
        #
        data = []
        #
        mineral = "Arf"
        #
        # Molar mass
        M = round(3*sodium[2] + 5*iron[2] + (8*silicon[2] + 22*oxygen[2]) + 2*(oxygen[2] + hydrogen[2]), 3)
        # Density
        data_V = CrystalPhysics([[9.9, 18.0, 5.3], [104.0], "monoclinic"])
        V = data_V.calculate_volume()
        Z = 2
        data_rho = CrystalPhysics([M, Z, V])
        rho = data_rho.calculate_bulk_density()
        # Bulk modulus
        K = 101.1*10**9
        # Shear modulus
        G = 43.7*10**9
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
        element = [hydrogen, oxygen, sodium, silicon, iron]
        w_H = round(2*hydrogen[2]/M, 4)
        w_O = round((22+2)*oxygen[2]/M, 4)
        w_Na = round(3*sodium[2]/M, 4)
        w_Si = round(8*silicon[2]/M, 4)
        w_Fe = round(5*iron[2]/M, 4)
        weights = [w_H, w_O, w_Na, w_Si, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(weights)
        #
        return data
    #
    def amphibole_ca(self):  # Tr + Act
        # Elements
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        calcium = elements.Ca(self)
        iron = elements.Fe(self)
        # Minerals
        tr = inosilicates.tremolite("")
        act = inosilicates.actinolite("")
        #
        data = []
        #
        mineral = "Amph"
        #
        # Molar mass
        x = rd.uniform(0.1, 0.4)
        M = x*tr[1] + (1-x)*act[1][0]
        # Density
        rho = x*tr[2] + (1-x)*act[2]
        # Bulk modulus
        K = (x*tr[3][0] + (1-x)*act[3][0])*10**9
        # Shear modulus
        G = (x*tr[3][1] + (1-x)*act[3][1])*10**9
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
        GR = x*tr[5][0] + (1-x)*act[5][0]
        # Photoelectricity
        element = [hydrogen, oxygen, magnesium, silicon, calcium, iron]
        w_H = round(x*tr[6][0] + (1-x)*act[6][0], 4)
        w_O = round(x*tr[6][1] + (1-x)*act[6][1], 4)
        w_Mg = round(x*tr[6][2] + (1-x)*act[6][2], 4)
        w_Si = round(x*tr[6][3] + (1-x)*act[6][3], 4)
        w_Ca = round(x*tr[6][4] + (1-x)*act[6][4], 4)
        w_Fe = round((1-x)*act[6][5], 4)
        composition = [w_H, w_O, w_Mg, w_Si, w_Ca, w_Fe]
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def amphibole_na(self):  # Rbk + Arf
        # Elements
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        sodium = elements.Na(self)
        silicon = elements.Si(self)
        iron = elements.Fe(self)
        # Minerals
        riebeckite = inosilicates.riebeckite("")
        arfvedsonite = inosilicates.arfvedsonite("")
        #
        data = []
        #
        mineral = "Amph"
        #
        # Molar mass
        x = rd.uniform(0.0, 1.0)
        M = x*riebeckite[1] + (1-x)*arfvedsonite[1]
        # Density
        rho = x*riebeckite[2] + (1-x)*arfvedsonite[2]
        # Bulk modulus
        K = (x*riebeckite[3][0] + (1-x)*arfvedsonite[3][0])*10**9
        # Shear modulus
        G = (x*riebeckite[3][1] + (1-x)*arfvedsonite[3][1])*10**9
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
        GR = x*riebeckite[5][0] + (1-x)*arfvedsonite[5][0]
        # Photoelectricity
        element = [hydrogen, oxygen, sodium, silicon, iron]
        w_H = round(x*riebeckite[6][0] + (1-x)*arfvedsonite[6][0], 4)
        w_O = round(x*riebeckite[6][1] + (1-x)*arfvedsonite[6][1], 4)
        w_Na = round(x*riebeckite[6][2] + (1-x)*arfvedsonite[6][2], 4)
        w_Si = round(x*riebeckite[6][3] + (1-x)*arfvedsonite[6][3], 4)
        w_Fe = round(x*riebeckite[6][4] + (1-x)*arfvedsonite[6][4], 4)
        composition = [w_H, w_O, w_Na, w_Si, w_Fe]
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def pyroxene(self):  # Ens + Aug
        # Elements
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        calcium = elements.Ca(self)
        iron = elements.Fe(self)
        # Minerals
        en = inosilicates.enstatite("")
        aug = inosilicates.augite("")
        #
        data = []
        #
        mineral = "Pyx"
        #
        # Molar mass
        x = rd.uniform(0.0, 0.3)
        M = x*en[1] + (1-x)*aug[1][0]
        # Density
        rho = x*en[2] + (1-x)*aug[2]
        # Bulk modulus
        K = (x*en[3][0] + (1-x)*aug[3][0])*10**9
        # Shear modulus
        G = (x*en[3][1] + (1-x)*aug[3][1])*10**9
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
        GR = x*en[5][0] + (1-x)*aug[5][0]
        # Photoelectricity
        element = [oxygen, magnesium, silicon, calcium, iron]
        w_O = round(x*en[6][0] + (1-x)*aug[6][0], 4)
        w_Mg = round(x*en[6][1] + (1-x)*aug[6][1], 4)
        w_Si = round(x*en[6][2] + (1-x)*aug[6][2], 4)
        w_Ca = round((1-x)*aug[6][3], 4)
        w_Fe = round((1-x)*aug[6][4], 4)
        composition = [w_O, w_Mg, w_Si, w_Ca, w_Fe]
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def pyroxene_ca(self):  # Aug + Di
        # Elements
        oxygen = elements.O(self)
        magnesium = elements.Mg(self)
        silicon = elements.Si(self)
        calcium = elements.Ca(self)
        iron = elements.Fe(self)
        # Minerals
        augite = inosilicates.augite("")
        diopside = inosilicates.diopside("")
        #
        data = []
        #
        mineral = "Pyx"
        #
        # Molar mass
        x = rd.uniform(0.5, 1.0)
        M = x*augite[1][0] + (1-x)*diopside[1]
        # Density
        rho = x*augite[2] + (1-x)*diopside[2]
        # Bulk modulus
        K = (x*augite[3][0] + (1-x)*diopside[3][0])*10**9
        # Shear modulus
        G = (x*augite[3][1] + (1-x)*diopside[3][1])*10**9
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
        GR = x*augite[5][0] + (1-x)*diopside[5][0]
        # Photoelectricity
        element = [oxygen, magnesium, silicon, calcium, iron]
        w_O = round(x*augite[6][0] + (1-x)*diopside[6][0], 4)
        w_Mg = round(x*augite[6][1] + (1-x)*diopside[6][1], 4)
        w_Si = round(x*augite[6][2] + (1-x)*diopside[6][2], 4)
        w_Ca = round(x*augite[6][3] + (1-x)*diopside[6][3], 4)
        w_Fe = round(x*augite[6][4], 4)
        composition = [w_O, w_Mg, w_Si, w_Ca, w_Fe]
        PE = bg.calculate_pe(self, x_list=composition, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
##################
# TECTOSILICATES #
##################
#
class tectosilicates:
    #
    def __init__(self):
        pass
    #
    def orthoclase(self):   # KAlSi3O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Or"
        #
        # Molar mass
        M = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho_Or = CrystalPhysics([M, 4, V])
        rho = dataRho_Or.calculate_bulk_density()
        # Bulk modulus
        K = 89.78*10**9
        # Shear modulus
        G = 61.52*10**9
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
        GR = potassium[2]/M*100*16
        # Photoelectricity
        PE = 2.85
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2)])
        #
        return data
    #
    def albite(self):   # NaAlSi3O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        chemAlbite = []
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        M = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)     # Molar mass
        chemAlbite.append(M)
        rho = round(np.random.normal(2600, 10), 1)               # Density
        chemAlbite.append(rho)
        K = round(np.random.normal(75.6, 1), 2)                  # Bulk modulus
        chemAlbite.append(K)
        G = round(np.random.normal(25.6, 1), 2)                   # Shear modulus
        chemAlbite.append(G)
        vP = ((K*10**9 + 4/3 * G*10**9)/rho)**0.5           # P-wave velocity
        chemAlbite.append(round(vP, 1))
        vS = ((G * 10**9)/rho)**0.5                           # S-wave velocity
        chemAlbite.append(round(vS, 1))
        GR = 0                                                  # Gamma ray
        chemAlbite.append(GR)
        PE = 1.75
        chemAlbite.append(PE)
        #
        return chemAlbite
    #
    def anorthite(self):   # CaAl2Si2O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        #
        chemAnorthite = []
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        M = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)     # Molar mass
        chemAnorthite.append(M)
        rho = round(np.random.normal(2730, 10), 1)               # Density
        chemAnorthite.append(rho)
        K = round(np.random.normal(91.8, 1), 2)                  # Bulk modulus
        chemAnorthite.append(K)
        G = round(np.random.normal(40.76, 1), 2)                   # Shear modulus
        chemAnorthite.append(G)
        vP = ((K*10**9 + 4/3 * G*10**9)/rho)**0.5           # P-wave velocity
        chemAnorthite.append(round(vP, 1))
        vS = ((G * 10**9)/rho)**0.5                           # S-wave velocity
        chemAnorthite.append(round(vS, 1))
        GR = 0                                                  # Gamma ray
        chemAnorthite.append(GR)
        PE = 3.05
        chemAnorthite.append(PE)
        #
        return chemAnorthite
    #
    def scapolite(self, keyword):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        carbon = elements.C(self)
        oxygen = elements.O(self)
        sodium = elements.Na(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        chem_Cl = elements.Cl(self)
        calcium = elements.Ca(self)
        element = [carbon, oxygen, sodium, aluminium, silicon, chem_Cl, calcium]
        self.keyword = keyword
        #
        if self.keyword in ["Scapolite","scapolite", "Scp", "scp"]:
            x_Ca = round(rd.uniform(0, 1), 2)
            x_Na = round(1-x_Ca, 2)
        elif self.keyword in ["Ca", "Meionite", "meionite", "Mei", "mei"]:
            x_Ca = round(rd.uniform(0.5, 1), 2)
            x_Na = round(1-x_Ca, 2)
        elif self.keyword in ["Na", "Marialite", "marialite", "Mar", "mar"]:
            x_Ca = round(rd.uniform(0, 0.5), 2)
            x_Na = round(1-x_Ca, 2)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], [GR, PE]]
        data = []
        #
        mineral = "Scp"
        #
        # Molar mass
        M = round(x_Ca*(4*calcium[2]+6*aluminium[2]+6*silicon[2]+24*oxygen[2]+carbon[2]+3*oxygen[2]) + x_Na*(4*sodium[2]+3*aluminium[2]+9*silicon[2]+24*oxygen[2]+chem_Cl[2]), 3)
        # Density
        M_Mei = round(4*calcium[2]+6*aluminium[2]+6*silicon[2]+24*oxygen[2]+carbon[2]+3*oxygen[2], 3)
        dataV_Mei = CrystalPhysics([[12.26, 7.61], [], "tetragonal"])
        V_Mei = dataV_Mei.calculate_volume()
        dataRho_Mei = CrystalPhysics([M_Mei, 2, V_Mei])
        rho_Mei = dataRho_Mei.calculate_bulk_density()
        M_Mar = round(4*sodium[2]+3*aluminium[2]+9*silicon[2]+24*oxygen[2]+chem_Cl[2], 3)
        dataV_Mar = CrystalPhysics([[12.075, 7.516], [], "tetragonal"])
        V_Mar = dataV_Mar.calculate_volume()
        dataRho_Mar = CrystalPhysics([M_Mar, 2, V_Mar])
        rho_Mar = dataRho_Mar.calculate_bulk_density()
        rho = x_Ca*rho_Mei + x_Na*rho_Mar
        # Bulk modulus
        K_Mei = rd.uniform(0.95, 1)*94.58*10**9
        K_Mar = 94.58*10**9
        K = x_Ca*K_Mei + x_Na*K_Mar
        # Shear modulus
        G_Mei = rd.uniform(0.95, 1)*62.83*10**9
        G_Mar = 62.83*10**9
        G = x_Ca*G_Mei + x_Na*G_Mar
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
        #PE_Mei = 3.33
        #PE_Mar = 1.87
        #PE = x_Ca*PE_Mei + x_Na*PE_Mar
        #element = [carbon, oxygen, sodium, aluminium, silicon, chem_Cl, calcium]
        w_C = round(x_Ca*carbon[2]/M, 4)
        w_O = round((x_Ca*(24+3)+x_Na*24)*oxygen[2]/M, 4)
        w_Na = round(x_Na*4*sodium[2]/M, 4)
        w_Al = round((x_Ca*6+x_Na*3)*aluminium[2]/M, 4)
        w_Si = round((x_Ca*6+x_Na*9)*silicon[2]/M, 4)
        w_Cl = round(x_Na*chem_Cl[2]/M, 4)
        w_Ca = round(x_Ca*4*calcium[2]/M, 4)
        weights = [w_C, w_O, w_Na, w_Al, w_Si, w_Cl, w_Ca]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x_Ca, x_Na])
        data.append(round(rho, 1))
        data.append([round(K * 10 ** (-9), 2), round(G * 10 ** (-9), 2), round(E * 10 ** (-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
#
##################
# FELDSPAR GROUP #
##################
#
class feldspars:
    #
    def __init__(self):
        pass
    #
    def alkalifeldspar(self, keyword="None"):  # Na(x)K(1-x)AlSi3O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        potassium = elements.K(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        element = [sodium, potassium, aluminium, silicon, oxygen]
        self.keyword = keyword
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        if self.keyword == "None":
            x = round(rd.uniform(0, 1), 2)
            if x >= 0.9:
                mineral = "Ab"
            elif x < 0.9 and x >= 0.63:
                mineral = "Ano"
            elif x < 0.63 and x > 0.1:
                mineral = "Sa"
            elif x <= 0.1:
                magicnumber = rd.randint(0, 1)
                if magicnumber == 0:
                    mineral = "Or"
                if magicnumber == 1:
                    mineral = "Mc"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45*10**9
            K_Or = 89.78*10**9
            K = (x*K_Ab + (1-x)*K_Or)
            # Shear modulus
            G_Ab = 69.0*10**9
            G_Or = 61.52*10**9
            G = (x*G_Ab + (1-x)*G_Or)
            # Young's modulus
            E = (9*K*G)/(3*K + G)
            E_Ab = (9*K_Ab*G_Ab)/(3*K_Ab + G_Ab)
            E_Or = (9*K_Or*G_Or)/(3*K_Or + G_Or)
            # Poisson's ratio
            nu = (3*K - 2*G)/(2*(3*K + G))
            nu_Ab = (3*K_Ab - 2*G_Ab)/(2*(3*K_Ab + G_Ab))
            nu_Or = (3*K_Or - 2*G_Or)/(2*(3*K_Or + G_Or))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP = ((K + 4/3*G)/rho)**0.5
            vP_Ab = ((K_Ab + 4/3*G_Ab)/(rho_Ab))**(0.5)
            vP_Or = ((K_Or + 4/3*G_Or)/(rho_Or))**(0.5)
            # S-wave velocity
            vS = (G/rho)**0.5
            vS_Ab = ((G_Ab)/(rho_Ab))**(0.5)
            vS_Or = ((G_Or)/(rho_Or))**(0.5)
            # Gamma ray
            GR = (1-x)*potassium[2]/M*100*16
            GR_Ab = (1-1)*potassium[2]/M_Ab*100*16
            GR_Or = (1-0)*potassium[2]/M_Or*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif self.keyword == "Na":
            mineral = "Kfs"
            x = round(rd.uniform(0.75, 1), 2)
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.63:
            #    mineral = "Ano"
            #elif x < 0.63 and x > 0.1:
            #    mineral = "Sa"
            #elif x <= 0.1:
            #    magicnumber = rd.randint(0, 1)
            #    if magicnumber == 0:
            #        mineral = "Or"
            #    if magicnumber == 1:
            #        mineral = "Mc"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45
            K_Or = 89.78
            K = (x*K_Ab + (1-x)*K_Or)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_Or = 61.52
            G = (x*G_Ab + (1-x)*G_Or)*10**9
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif self.keyword == "K":
            mineral = "Kfs"
            x = round(rd.uniform(0, 0.25), 2)
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.63:
            #    mineral = "Ano"
            #elif x < 0.63 and x > 0.1:
            #    mineral = "Sa"
            #elif x <= 0.1:
            #    magicnumber = rd.randint(0, 1)
            #    if magicnumber == 0:
            #        mineral = "Or"
            #    if magicnumber == 1:
            #        mineral = "Mc"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45
            K_Or = 89.78
            K = (x*K_Ab + (1-x)*K_Or)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_Or = 61.52
            G = (x*G_Ab + (1-x)*G_Or)*10**9
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        if self.keyword == "Alkalifeldspar" or self.keyword == "Afs" or self.keyword == "Kfs":
            x = round(rd.uniform(0, 1), 2)
            mineral = "Kfs"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45
            K_Or = 89.78
            K = (x*K_Ab + (1-x)*K_Or)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_Or = 61.52
            G = (x*G_Ab + (1-x)*G_Or)*10**9
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        #
        return data
    #
    def plagioclase(self, keyword="None"):  # Na(x)Ca(1-x)Al(2-x)Si(2+x)O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        calcium = elements.Ca(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        oxygen = elements.O(self)
        element = [sodium, calcium, aluminium, silicon, oxygen]
        self.keyword = keyword
        #
        data = []
        #
        if self.keyword == "None":
            x = round(rd.uniform(0, 1), 2)
            if x >= 0.9:
                mineral = "Ab"
            elif x < 0.9 and x >= 0.7:
                mineral = "Olg"
            elif x < 0.7 and x >= 0.5:
                mineral = "Andes"
            elif x < 0.5 and x >= 0.3:
                mineral = "Lab"
            elif x < 0.3 and x > 0.1:
                mineral = "Byt"
            elif x <= 0.1:
                mineral = "An"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45*10**9
            K_An = 114.72*10**9
            K = (x*K_Ab + (1-x)*K_An)
            # Shear modulus
            G_Ab = 69.0*10**9
            G_An = 73.72*10**9
            G = (x*G_Ab + (1-x)*G_An)
            # Young's modulus
            E = (9*K*G)/(3*K + G)
            E_Ab = (9*K_Ab*G_Ab)/(3*K_Ab + G_Ab)
            E_An = (9*K_An*G_An)/(3*K_An + G_An)
            # Poisson's ratio
            nu = (3*K - 2*G)/(2*(3*K + G))
            nu_Ab = (3*K_Ab - 2*G_Ab)/(2*(3*K_Ab + G_Ab))
            nu_An = (3*K_An - 2*G_An)/(2*(3*K_An + G_An))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP = ((K + 4/3*G)/rho)**0.5
            vP_Ab = ((K_Ab + 4/3*G_Ab)/(rho_Ab))**(0.5)
            vP_An = ((K_An + 4/3*G_An)/(rho_An))**(0.5)
            # S-wave velocity
            vS = (G/rho)**0.5
            vS_Ab = ((G_Ab)/(rho_Ab))**(0.5)
            vS_An = ((G_An)/(rho_An))**(0.5)
            # Gamma ray
            GR = 0
            GR_Ab = 0
            GR_An = 0
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif self.keyword == "Na":
            x = round(rd.uniform(0.75, 1), 2)
            mineral = "Pl"
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.7:
            #    mineral = "Olg"
            #elif x < 0.7 and x >= 0.5:
            #    mineral = "Andes"
            #elif x < 0.5 and x >= 0.3:
            #    mineral = "Lab"
            #elif x < 0.3 and x > 0.1:
            #    mineral = "Byt"
            #elif x <= 0.1:
            #    mineral = "An"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45
            K_An = 114.72
            K = (x*K_Ab + (1-x)*K_An)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_An = 73.72
            G = (x*G_Ab + (1-x)*G_An)*10**9
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
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif self.keyword == "Ca":
            x = round(rd.uniform(0, 0.25), 2)
            mineral = "Pl"
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.7:
            #    mineral = "Olg"
            #elif x < 0.7 and x >= 0.5:
            #    mineral = "Andes"
            #elif x < 0.5 and x >= 0.3:
            #    mineral = "Lab"
            #elif x < 0.3 and x > 0.1:
            #    mineral = "Byt"
            #elif x <= 0.1:
            #    mineral = "An"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45
            K_An = 114.72
            K = (x*K_Ab + (1-x)*K_An)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_An = 73.72
            G = (x*G_Ab + (1-x)*G_An)*10**9
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
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif self.keyword == "Plagioclase" or self.keyword == "Pl":
            x = round(rd.uniform(0, 1), 2)
            mineral = "Pl"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45
            K_An = 114.72
            K = (x*K_Ab + (1-x)*K_An)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_An = 73.72
            G = (x*G_Ab + (1-x)*G_An)*10**9
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
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        #
        return data
#
class Feldspars:
    #
    def alkalifeldspar(keyword="None"):  # Na(x)K(1-x)AlSi3O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na("")
        potassium = elements.K("")
        aluminium = elements.Al("")
        silicon = elements.Si("")
        oxygen = elements.O("")
        element = [sodium, potassium, aluminium, silicon, oxygen]
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        if keyword == "None":
            x = round(rd.uniform(0, 1), 2)
            if x >= 0.9:
                mineral = "Ab"
            elif x < 0.9 and x >= 0.63:
                mineral = "Ano"
            elif x < 0.63 and x > 0.1:
                mineral = "Sa"
            elif x <= 0.1:
                magicnumber = rd.randint(0, 1)
                if magicnumber == 0:
                    mineral = "Or"
                if magicnumber == 1:
                    mineral = "Mc"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45*10**9
            K_Or = 89.78*10**9
            K = (x*K_Ab + (1-x)*K_Or)
            # Shear modulus
            G_Ab = 69.0*10**9
            G_Or = 61.52*10**9
            G = (x*G_Ab + (1-x)*G_Or)
            # Young's modulus
            E = (9*K*G)/(3*K + G)
            E_Ab = (9*K_Ab*G_Ab)/(3*K_Ab + G_Ab)
            E_Or = (9*K_Or*G_Or)/(3*K_Or + G_Or)
            # Poisson's ratio
            nu = (3*K - 2*G)/(2*(3*K + G))
            nu_Ab = (3*K_Ab - 2*G_Ab)/(2*(3*K_Ab + G_Ab))
            nu_Or = (3*K_Or - 2*G_Or)/(2*(3*K_Or + G_Or))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP = ((K + 4/3*G)/rho)**0.5
            vP_Ab = ((K_Ab + 4/3*G_Ab)/(rho_Ab))**(0.5)
            vP_Or = ((K_Or + 4/3*G_Or)/(rho_Or))**(0.5)
            # S-wave velocity
            vS = (G/rho)**0.5
            vS_Ab = ((G_Ab)/(rho_Ab))**(0.5)
            vS_Or = ((G_Or)/(rho_Or))**(0.5)
            # Gamma ray
            GR = (1-x)*potassium[2]/M*100*16
            GR_Ab = (1-1)*potassium[2]/M_Ab*100*16
            GR_Or = (1-0)*potassium[2]/M_Or*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif keyword == "Na":
            mineral = "Kfs"
            x = round(rd.uniform(0.6, 1), 2)
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.63:
            #    mineral = "Ano"
            #elif x < 0.63 and x > 0.1:
            #    mineral = "Sa"
            #elif x <= 0.1:
            #    magicnumber = rd.randint(0, 1)
            #    if magicnumber == 0:
            #        mineral = "Or"
            #    if magicnumber == 1:
            #        mineral = "Mc"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45
            K_Or = 89.78
            K = (x*K_Ab + (1-x)*K_Or)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_Or = 61.52
            G = (x*G_Ab + (1-x)*G_Or)*10**9
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif keyword == "K":
            mineral = "Kfs"
            x = round(rd.uniform(0, 0.4), 2)
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.63:
            #    mineral = "Ano"
            #elif x < 0.63 and x > 0.1:
            #    mineral = "Sa"
            #elif x <= 0.1:
            #    magicnumber = rd.randint(0, 1)
            #    if magicnumber == 0:
            #        mineral = "Or"
            #    if magicnumber == 1:
            #        mineral = "Mc"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45
            K_Or = 89.78
            K = (x*K_Ab + (1-x)*K_Or)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_Or = 61.52
            G = (x*G_Ab + (1-x)*G_Or)*10**9
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        if keyword == "Alkalifeldspar":
            x = round(rd.uniform(0, 1), 2)
            mineral = "Kfs"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
            V_Or = dataV_Or.calculate_volume()
            dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
            rho_Or = dataRho_Or.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_Or
            # Bulk modulus
            K_Ab = 103.45
            K_Or = 89.78
            K = (x*K_Ab + (1-x)*K_Or)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_Or = 61.52
            G = (x*G_Ab + (1-x)*G_Or)*10**9
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
            GR = (1-x)*potassium[2]/M*100*16
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_Or = 2.85
            #PE = x*PE_Ab + (1-x)*PE_Or
            #element = [sodium, potassium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_K = round((1-x)*potassium[2]/M, 4)
            w_Al = round(aluminium[2]/M, 4)
            w_Si = round(3*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_K, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_K]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        #
        return data
    #
    def plagioclase(keyword="None"):  # Na(x)Ca(1-x)Al(2-x)Si(2+x)O8
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na("")
        calcium = elements.Ca("")
        aluminium = elements.Al("")
        silicon = elements.Si("")
        oxygen = elements.O("")
        element = [sodium, calcium, aluminium, silicon, oxygen]
        keyword = keyword
        #
        data = []
        #
        if keyword == "None":
            x = round(rd.uniform(0, 1), 2)
            if x >= 0.9:
                mineral = "Ab"
            elif x < 0.9 and x >= 0.7:
                mineral = "Olg"
            elif x < 0.7 and x >= 0.5:
                mineral = "Andes"
            elif x < 0.5 and x >= 0.3:
                mineral = "Lab"
            elif x < 0.3 and x > 0.1:
                mineral = "Byt"
            elif x <= 0.1:
                mineral = "An"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45*10**9
            K_An = 114.72*10**9
            K = (x*K_Ab + (1-x)*K_An)
            # Shear modulus
            G_Ab = 69.0*10**9
            G_An = 73.72*10**9
            G = (x*G_Ab + (1-x)*G_An)
            # Young's modulus
            E = (9*K*G)/(3*K + G)
            E_Ab = (9*K_Ab*G_Ab)/(3*K_Ab + G_Ab)
            E_An = (9*K_An*G_An)/(3*K_An + G_An)
            # Poisson's ratio
            nu = (3*K - 2*G)/(2*(3*K + G))
            nu_Ab = (3*K_Ab - 2*G_Ab)/(2*(3*K_Ab + G_Ab))
            nu_An = (3*K_An - 2*G_An)/(2*(3*K_An + G_An))
            # vP/vS
            vPvS = ((K + 4/3*G)/G)**0.5
            # P-wave velocity
            vP = ((K + 4/3*G)/rho)**0.5
            vP_Ab = ((K_Ab + 4/3*G_Ab)/(rho_Ab))**(0.5)
            vP_An = ((K_An + 4/3*G_An)/(rho_An))**(0.5)
            # S-wave velocity
            vS = (G/rho)**0.5
            vS_Ab = ((G_Ab)/(rho_Ab))**(0.5)
            vS_An = ((G_An)/(rho_An))**(0.5)
            # Gamma ray
            GR = 0
            GR_Ab = 0
            GR_An = 0
            # Photoelectricity
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif keyword == "Na":
            x = round(rd.uniform(0.6, 1), 2)
            mineral = "Pl"
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.7:
            #    mineral = "Olg"
            #elif x < 0.7 and x >= 0.5:
            #    mineral = "Andes"
            #elif x < 0.5 and x >= 0.3:
            #    mineral = "Lab"
            #elif x < 0.3 and x > 0.1:
            #    mineral = "Byt"
            #elif x <= 0.1:
            #    mineral = "An"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45
            K_An = 114.72
            K = (x*K_Ab + (1-x)*K_An)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_An = 73.72
            G = (x*G_Ab + (1-x)*G_An)*10**9
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
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif keyword == "Ca":
            x = round(rd.uniform(0, 0.4), 2)
            mineral = "Pl"
            #if x >= 0.9:
            #    mineral = "Ab"
            #elif x < 0.9 and x >= 0.7:
            #    mineral = "Olg"
            #elif x < 0.7 and x >= 0.5:
            #    mineral = "Andes"
            #elif x < 0.5 and x >= 0.3:
            #    mineral = "Lab"
            #elif x < 0.3 and x > 0.1:
            #    mineral = "Byt"
            #elif x <= 0.1:
            #    mineral = "An"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45
            K_An = 114.72
            K = (x*K_Ab + (1-x)*K_An)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_An = 73.72
            G = (x*G_Ab + (1-x)*G_An)*10**9
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
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        elif keyword == "Plagioclase":
            x = round(rd.uniform(0, 1), 2)
            mineral = "Pl"
            #
            # Molar mass
            M = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
            # Density
            M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
            dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
            V_Ab = dataV_Ab.calculate_volume()
            dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
            rho_Ab = dataRho_Ab.calculate_bulk_density()
            M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
            dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
            V_An = dataV_An.calculate_volume()
            dataRho_An = CrystalPhysics([M_An, 8, V_An])
            rho_An = dataRho_An.calculate_bulk_density()
            rho = x*rho_Ab + (1-x)*rho_An
            # Bulk modulus
            K_Ab = 103.45
            K_An = 114.72
            K = (x*K_Ab + (1-x)*K_An)*10**9
            # Shear modulus
            G_Ab = 69.0
            G_An = 73.72
            G = (x*G_Ab + (1-x)*G_An)*10**9
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
            #PE_Ab = 1.75
            #PE_An = 3.05
            #PE = x*PE_Ab + (1-x)*PE_An
            #element = [sodium, calcium, aluminium, silicon, oxygen]
            w_Na = round(x*sodium[2]/M, 4)
            w_Ca = round((1-x)*calcium[2]/M, 4)
            w_Al = round((2-x)*aluminium[2]/M, 4)
            w_Si = round((2+x)*silicon[2]/M, 4)
            w_O = round(8*oxygen[2]/M, 4)
            weights = [w_Na, w_Ca, w_Al, w_Si, w_O]
            composition = [w_O, w_Na, w_Al, w_Si, w_Ca]
            PE = bg.calculate_pe("", x_list=weights, elements_list=element)
            U = PE*rho*10**(-3)
            # Electrical resistivity
            p = 5.5*10**11
            #
            data.append(mineral)
            data.append([round(M, 2), x])
            data.append(round(rho, 1))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
            data.append([round(vP, 1), round(vS, 1)])
            data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
            data.append(composition)
        #
        return data
#
###########
# HALIDES #
###########
#
class halides:
    #
    def __init__(self):
        pass
    #
    def halite(self):   # NaCl
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        chemCl = elements.Cl(self)
        element = [sodium, chemCl]
        #
        data = []
        #
        mineral = "Hl"
        #
        # Molar mass
        M = round(sodium[2] + chemCl[2], 3)
        x_Na = round(sodium[2]/M, 4)
        x_Cl = round(chemCl[2]/M, 4)
        weights = [x_Na, x_Cl]
        # Density
        dataV = CrystalPhysics([[5.6404], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 23*10**9
        # Shear modulus
        G = 14.33*10**9
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
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(weights)
        #
        return data
    #
    def sylvite(self):  # KCl
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        chemCl = elements.Cl(self)
        #
        data = []
        #
        mineral = "Syl"
        #
        # Molar mass
        M = round(potassium[2] + chemCl[2], 3)
        # Density
        dataV = CrystalPhysics([[6.2931], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = potassium[2]/M*100*16
        # Photoelectricity
        element = [potassium, chemCl]
        w_K = round(potassium[2]/M, 4)
        w_Cl = round(chemCl[2]/M, 4)
        weights = [w_K, w_Cl]
        composition = [w_Cl, w_K]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def bischofite(self):   # MgCl2*6H2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        magnesium = elements.Mg(self)
        chem_Cl = elements.Cl(self)
        hydrogen = elements.H(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Bis"
        #
        # Molar mass
        M = round(magnesium[2] + chem_Cl[2] + 6*(2*hydrogen[2]+oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[9.9, 7.15, 6.1], [93.7], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 61.98*10**9
        # Shear modulus
        G = 31.66*10**9
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
        element = [magnesium, chem_Cl, hydrogen, oxygen]
        w_Mg = round(magnesium[2]/M, 4)
        w_Cl = round(chem_Cl[2]/M, 4)
        w_H = round(6*2*hydrogen[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        weights = [w_Mg, w_Cl, w_H, w_O]
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
    def fluorite(self):   # CaF2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        flourine = elements.F(self)
        #
        data = []
        #
        mineral = "Fl"
        #
        # Molar mass
        M = round(calcium[2] + 2*flourine[2], 3)
        # Density
        dataV = CrystalPhysics([[5.4626], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 77*10**9
        # Shear modulus
        G = 38*10**9
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
        element = [calcium, flourine]
        w_Ca = round(calcium[2]/M, 4)
        w_F = round(2*flourine[2]/M, 4)
        weights = [w_Ca, w_F]
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
    def carobbiite(self):   # KF
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        potassium = elements.K(self)
        flourine = elements.F(self)
        #
        data = []
        #
        mineral = "Crb"
        #
        # Molar mass
        M = round(potassium[2] + flourine[2], 3)
        # Density
        dataV = CrystalPhysics([[5.34], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 29*10**9
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
        GR = 935.96
        # Photoelectricity
        element = [potassium, flourine]
        w_K = round(potassium[2]/M, 4)
        w_F = round(flourine[2]/M, 4)
        weights = [w_K, w_F]
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
############
# SULFIDES #
############
#
class sulfides:
    #
    def __init__(self):
        pass
    #
    def pyrite(self):   # FeS2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Py"
        #
        # Molar mass
        M = round(iron[2] + 2*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.417], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [iron, sulfur]
        w_Fe = round(iron[2]/M, 4)
        w_S = round(2*sulfur[2]/M, 4)
        weights = [w_Fe, w_S]
        composition = [w_S, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        # Electrical resistivity
        p = 3*10**(-1)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def bornite(self):   # Cu5FeS4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCu = elements.Cu(self)
        iron = elements.Fe(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Bn"
        #
        # Molar mass
        M = round(5*chemCu[2] + iron[2] + 4*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[10.95, 21.862, 10.95], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 16, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [chemCu, iron, sulfur]
        w_Cu = round(5*chemCu[2]/M, 4)
        w_Fe = round(iron[2]/M, 4)
        w_S = round(4*sulfur[2]/M, 4)
        weights = [w_Cu, w_Fe, w_S]
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
    def galena(self):   # PbS
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemPb = elements.Pb(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Gn"
        #
        # Molar mass
        M = round(chemPb[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.936], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [chemPb, sulfur]
        w_Pb = round(chemPb[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Pb, w_S]
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
    def chalcopyrite(self):   # CuFeS2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCu = elements.Cu(self)
        iron = elements.Fe(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Ccp"
        #
        # Molar mass
        M = round(chemCu[2] + iron[2] + 2*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.28, 10.41], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [chemCu, iron, sulfur]
        w_Cu = round(chemCu[2]/M, 4)
        w_Fe = round(iron[2]/M, 4)
        w_S = round(2*sulfur[2]/M, 4)
        weights = [w_Cu, w_Fe, w_S]
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
    def molybdenite(self):   # MoS2
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        molybdenum = elements.Mo(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Mol"
        #
        # Molar mass
        M = round(molybdenum[2] + 2*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[3.16, 12.3], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [molybdenum, sulfur]
        w_Mo = round(molybdenum[2]/M, 4)
        w_S = round(2*sulfur[2]/M, 4)
        weights = [w_Mo, w_S]
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
    def sphalerite(self):   # ZnS
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemZn = elements.Zn(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Sp"
        #
        # Molar mass
        M = round(chemZn[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.406], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [chemZn, sulfur]
        w_Zn = round(chemZn[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Zn, w_S]
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
    def stibnite(self):   # Sb2S3
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        antimony = elements.Sb(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Stbn"
        #
        # Molar mass
        M = round(2*antimony[2] + 3*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[11.229,11.31, 3.893], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [antimony, sulfur]
        w_Sb = round(2*antimony[2]/M, 4)
        w_S = round(3*sulfur[2]/M, 4)
        weights = [w_Sb, w_S]
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
    def arsenopyrite(self):   # FeAsS
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        arsenic = elements.As(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Apy"
        #
        # Molar mass
        M = round(iron[2] + arsenic[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.74, 5.68, 5.79], [112.17], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [iron, arsenic, sulfur]
        w_Fe = round(iron[2]/M, 4)
        w_As = round(arsenic[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Fe, w_As, w_S]
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
    def acanthite(self):   # Ag2S
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        silver = elements.Ag(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Ach"
        #
        # Molar mass
        M = round(2*silver[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[4.229, 6.931, 7.862], [99.61], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [silver, sulfur]
        w_Ag = round(2*silver[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Ag, w_S]
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
    def argentite(self):   # Ag2S
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        silver = elements.Ag(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Argt"
        #
        # Molar mass
        M = round(2*silver[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[4.89], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [silver, sulfur]
        w_Ag = round(2*silver[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Ag, w_S]
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
    def alabandite(self):   # MnS
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        manganese = elements.Mn(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Ab"
        #
        # Molar mass
        M = round(manganese[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.214], [99.61], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [manganese, sulfur]
        w_Mn = round(manganese[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Mn, w_S]
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
    def berthierite(self):   # FeSb2S4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        antimony = elements.Sb(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Brt"
        #
        # Molar mass
        M = round(iron[2] + 2*antimony[2] + 4*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[11.44, 14.12, 3.76], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [iron, antimony, sulfur]
        w_Fe = round(iron[2]/M, 4)
        w_Sb = round(2*antimony[2]/M, 4)
        w_S = round(4*sulfur[2]/M, 4)
        weights = [w_Fe, w_Sb, w_S]
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
    def pyrrhotite(self):   # Fe(1-x)S
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Po"
        #
        # Molar mass
        x = round(rd.uniform(0.0, 0.17), 3)
        M = round((1-x)*iron[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[12.811, 6.87,11.885], [117, 3], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 26, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [iron, sulfur]
        w_Fe = round((1-x)*iron[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Fe, w_S]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2),round(x, 2)])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
    #
    def cobaltite(self):   # CoAsS
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        cobalt = elements.Co(self)
        arsenic = elements.As(self)
        sulfur = elements.S(self)
        #
        data = []
        #
        mineral = "Cob"
        #
        # Molar mass
        M = round(cobalt[2] + arsenic[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[5.57, 5.582, 5.582], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 139*10**9
        # Shear modulus
        G = 104*10**9
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
        element = [cobalt, arsenic, sulfur]
        w_Co = round(cobalt[2]/M, 4)
        w_As = round(arsenic[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        weights = [w_Co, w_As, w_S]
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
    def carrollite(self):   # CuCo2S4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        copper = elements.Cu(self)
        cobalt = elements.Co(self)
        sulfur = elements.S(self)
        element = [copper, cobalt, sulfur]
        #
        data = []
        #
        mineral = "Carr"
        #
        # Molar mass
        M = round(copper[2] + 2*cobalt[2] + 4*sulfur[2], 3)
        w_Cu = round(copper[2]/M, 4)
        w_Co = round(2*cobalt[2]/M, 4)
        w_S = round(4*sulfur[2]/M, 4)
        amounts = [1, 2, 4]
        # Density
        dataV = CrystalPhysics([[9.477], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
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
        GR = 0
        # Photoelectricity
        weights = [w_Cu, w_Co, w_S]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 0
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
    def chalcocite(self):   # Cu2S
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        copper = elements.Cu(self)
        sulfur = elements.S(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Cc"
        #
        # Molar mass
        M = round(2*copper[2] + sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[11.881, 27.323, 13.491], [116.35], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 96, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [copper, sulfur]
        x_Cu = round(2*copper[2]/M, 4)
        x_S = round(1*sulfur[2]/M, 4)
        weights = [x_Cu, x_S]
        amounts = [2, 1]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
    def digenite(self):   # Cu9S5
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        copper = elements.Cu(self)
        sulfur = elements.S(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Dg"
        #
        # Molar mass
        M = round(9*copper[2] + 5*sulfur[2], 3)
        # Density
        dataV = CrystalPhysics([[3.919, 48], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
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
        GR = 0
        # Photoelectricity
        element = [copper, sulfur]
        x_Cu = round(9*copper[2]/M, 4)
        x_S = round(5*sulfur[2]/M, 4)
        weights = [x_Cu, x_S]
        amounts = [9, 5]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
############
# SULFATES #
############
#
class sulfates:
    #
    def __init__(self):
        pass
    #
    def anhydrite(self):   # CaSO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        sulfur = elements.S(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Anh"
        #
        # Molar mass
        M = round(calcium[2] + sulfur[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[6.245, 6.995, 6.993], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 93.14*10**9
        # Shear modulus
        G = 50.15*10**9
        # P-wave velocity
        vP = ((K + 4/3 * G)/(rho))**(0.5)
        # S-wave velocity
        vS = (G/rho)**0.5
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
        element = [calcium, sulfur, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_Ca, w_S, w_O]
        composition = [w_O, w_S, w_Ca]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def gypsum(self):   # CaSO4*H2O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        sulfur = elements.S(self)
        oxygen = elements.O(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Gp"
        #
        # Molar mass
        M = round(calcium[2] + sulfur[2] + 4*oxygen[2] + 2*(2*hydrogen[2] + oxygen[2]), 3)
        # Density
        dataV = CrystalPhysics([[5.68, 15.18, 6.29], [113.83], "monoclinic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = 0.875*dataRho.calculate_bulk_density()
        # Bulk modulus
        K = round(np.random.normal(44, 0.4), 2)*10**9
        # Shear modulus
        G = round(np.random.normal(21, 0.2), 2)*10**9
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
        element = [calcium, sulfur, oxygen, hydrogen]
        w_Ca = round(calcium[2]/M, 4)
        w_S = round(sulfur[2]/M, 4)
        w_O = round((4+2)*oxygen[2]/M, 4)
        w_H = round(2*2*hydrogen[2]/M, 4)
        weights = [w_Ca, w_S, w_O, w_H]
        composition = [w_H, w_O, w_S, w_Ca]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def scheelite(self):   # CaWO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        tungsten = elements.W(self)
        oxygen = elements.O(self)
        #
        data = []
        #
        mineral = "Sch"
        #
        # Molar mass
        M = round(calcium[2] + tungsten[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[5.242, 11.372], [], "tetragonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 117.44*10**9
        # Shear modulus
        G = 51.64*10**9
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
        element = [calcium, tungsten, oxygen]
        w_Ca = round(calcium[2]/M, 4)
        w_W = round(tungsten[2]/M, 4)
        w_O = round(4*oxygen[2]/M, 4)
        weights = [w_Ca, w_W, w_O]
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
    def barite(self):   # BaSO4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        barium = elements.Ba(self)
        sulfur = elements.S(self)
        oxygen = elements.O(self)
        #
        # data = [ mineral, molar mass, density, [K, G, E, nu, vP/vS], [vP, vS], GR ]
        data = []
        #
        mineral = "Bar"
        #
        # Molar mass
        M = round(1*barium[2] + 1*sulfur[2] + 4*oxygen[2], 3)
        # Density
        dataV = CrystalPhysics([[8.878, 5.45, 7.152], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 52*10**9
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
        GR = 0
        # Photoelectricity
        element = [barium, sulfur, oxygen]
        x_Ba = round(1*barium[2]/M, 4)
        x_S = round(1*sulfur[2]/M, 4)
        x_O = round(4*oxygen[2]/M, 4)
        weights = [x_Ba, x_S, x_O]
        amounts = [1, 1, 4]
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
###############
# ANTIMONIDES #
###############
#
class antimonides:
    #
    def __init__(self):
        pass
    #
    def allargentum(self):   # AgSb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        silver = elements.Ag(self)
        antimony = elements.Sb(self)
        #
        data = []
        #
        mineral = "Alrg"
        #
        # Molar mass
        x = round(rd.uniform(0.009, 0.16), 2)
        M = round((1-x)*silver[2] + x*antimony[2], 3)
        # Density
        dataV = CrystalPhysics([[2.945, 4.77], [], "hexagonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 2, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = np.random.normal(63, 1)*10**9    # estimated
        # Shear modulus
        G = np.random.normal(23, 1)*10**9    # estimated
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
        element = [silver, antimony]
        w_Ag = round((1-x)*silver[2]/M, 4)
        w_Sb = round(x*antimony[2]/M, 4)
        weights = [w_Ag, w_Sb]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
    #
    def dyscrasite(self):   # Ag3Sb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        silver = elements.Ag(self)
        antimony = elements.Sb(self)
        #
        data = []
        #
        mineral = "Dys"
        #
        # Molar mass
        M = round(3*silver[2] + antimony[2], 3)
        # Density
        dataV = CrystalPhysics([[2.996, 5.235, 4.83], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 1, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 65*10**9    # estimated
        # Shear modulus
        G = 17*10**9    # estimated
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
        element = [silver, antimony]
        w_Ag = round(3*silver[2]/M, 4)
        w_Sb = round(antimony[2]/M, 4)
        weights = [w_Ag, w_Sb]
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
    def antimony(self):   # Sb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        antimony = elements.Sb(self)
        #
        data = []
        #
        mineral = "Sb"
        #
        # Molar mass
        M = round(antimony[2], 3)
        # Density
        dataV = CrystalPhysics([[4.299, 11.25], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 36*10**9    # estimated
        # Shear modulus
        G = 24*10**9    # estimated
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
        element = [antimony]
        w_Sb = round(antimony[2]/M, 4)
        weights = [w_Sb]
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
#############
# ARSENIDES #
#############
#
class Arsenides:
    #
    def __init__(self):
        pass
    #
#
####################
# TOURMALINE GROUP #
####################
#
class Tourmalines:
    #
    def __init__(self):
        pass
    #
    def schorl(self):   # NaFe3Al6(BO3)3Si6O18(OH)4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        hydrogen = elements.H(self)
        boron = elements.B(self)
        oxygen = elements.O(self)
        sodium = elements.Na(self)
        aluminium = elements.Al(self)
        silicon = elements.Si(self)
        iron = elements.Fe(self)
        #
        data = []
        #
        mineral = "Tur"
        #
        # Molar mass
        M = round(sodium[2] + 3*iron[2] + 6*aluminium[2] + 3*(boron[2]+3*oxygen[2]) + 6*silicon[2] + 18*oxygen[2] + 4*(oxygen[2]+hydrogen[2]), 3)
        # Density
        dataV = CrystalPhysics([[15.99, 7.195], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
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
        element = [sodium, iron, aluminium, boron, oxygen, silicon, hydrogen]
        w_Na = round(sodium[2]/M, 4)
        w_Fe = round(3*iron[2]/M, 4)
        w_Al = round(6*aluminium[2]/M, 4)
        w_B = round(3*boron[2]/M, 4)
        w_O = round((3*3+18+4)*oxygen[2]/M, 4)
        w_Si = round(6*silicon[2]/M, 4)
        w_H = round(4*hydrogen[2]/M, 4)
        weights = [w_Na, w_Fe, w_Al, w_B, w_O, w_Si, w_H]
        composition = [w_H, w_B, w_O, w_Na, w_Al, w_Si, w_Fe]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        data.append(composition)
        #
        return data
    #
    def dravite(self):   # NaMg3Al6(BO3)3Si6O18(OH)4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        magnesium = elements.Mg(self)
        aluminium = elements.Al(self)
        boron = elements.B(self)
        oxygen = elements.O(self)
        silicon = elements.Si(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Tur"
        #
        # Molar mass
        M = round(sodium[2] + 3*magnesium[2] + 6*aluminium[2] + 3*(boron[2]+3*oxygen[2]) + 6*silicon[2] + 18*oxygen[2] + 4*(oxygen[2]+hydrogen[2]), 3)
        # Density
        dataV = CrystalPhysics([[15.941, 7.201], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
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
        element = [sodium, magnesium, aluminium, boron, oxygen, silicon, hydrogen]
        w_Na = round(sodium[2]/M, 4)
        w_Mg = round(3*magnesium[2]/M, 4)
        w_Al = round(6*aluminium[2]/M, 4)
        w_B = round(3*boron[2]/M, 4)
        w_O = round((3*3+18+4)*oxygen[2]/M, 4)
        w_Si = round(6*silicon[2]/M, 4)
        w_H = round(4*hydrogen[2]/M, 4)
        weights = [w_Na, w_Mg, w_Al, w_B, w_O, w_Si, w_H]
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
    def elbaite(self):   # NaLi2.5Al6.5(BO3)3Si6O18(OH)4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        sodium = elements.Na(self)
        lithium = elements.Li(self)
        aluminium = elements.Al(self)
        boron = elements.B(self)
        oxygen = elements.O(self)
        silicon = elements.Si(self)
        hydrogen = elements.H(self)
        #
        data = []
        #
        mineral = "Tur"
        #
        # Molar mass
        M = round(sodium[2] + 2.5*lithium[2] + 6.5*aluminium[2] + 3*(boron[2]+3*oxygen[2]) + 6*silicon[2] + 18*oxygen[2] + 4*(oxygen[2]+hydrogen[2]), 3)
        # Density
        dataV = CrystalPhysics([[15.838, 7.103], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
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
        element = [sodium, lithium, aluminium, boron, oxygen, silicon, hydrogen]
        w_Na = round(sodium[2]/M, 4)
        w_Li = round(2.5*lithium[2]/M, 4)
        w_Al = round(6.5*aluminium[2]/M, 4)
        w_B = round(3*boron[2]/M, 4)
        w_O = round((3*3+18+4)*oxygen[2]/M, 4)
        w_Si = round(6*silicon[2]/M, 4)
        w_H = round(4*hydrogen[2]/M, 4)
        weights = [w_Na, w_Li, w_Al, w_B, w_O, w_Si, w_H]
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
    def liddicoatite(self):   # Ca(Li,Al)3Al6(BO3)3Si6O18(O,OH,F)4
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        calcium = elements.Ca(self)
        lithium = elements.Li(self)
        aluminium = elements.Al(self)
        boron = elements.B(self)
        oxygen = elements.O(self)
        silicon = elements.Si(self)
        hydrogen = elements.H(self)
        flourine = elements.F(self)
        #
        data = []
        #
        mineral = "Tur"
        #
        # Molar mass
        x = round(rd.uniform(0.0, 1.0), 2)
        y = round(rd.uniform(0.0, 1.0), 2)
        M = round(calcium[2] + 3*(x*lithium[2]+(1-x)*aluminium[2]) + 6*aluminium[2] + 3*(boron[2]+3*oxygen[2]) + 6*silicon[2] + 18*oxygen[2] + 4*(y*(x*oxygen[2]+(1-x)*(oxygen[2]+hydrogen[2]))+(1-y)*flourine[2]), 3)
        # Density
        dataV = CrystalPhysics([[15.875, 7.126], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
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
        element = [calcium, lithium, aluminium, boron, oxygen, silicon, hydrogen, flourine]
        w_Ca = round(calcium[2]/M, 4)
        w_Li = round(3*x*lithium[2]/M, 4)
        w_Al = round((3*(1-x)+6)*aluminium[2]/M, 4)
        w_B = round(3*boron[2]/M, 4)
        w_O = round((3*3+18+4*y*(x+(1-x)))*oxygen[2]/M, 4)
        w_Si = round(6*silicon[2]/M, 4)
        w_H = round(4*y*(1-x)*hydrogen[2]/M, 4)
        w_F = round(4*(1-y)*flourine[2]/M, 4)
        weights = [w_Ca, w_Li, w_Al, w_B, w_O, w_Si, w_H, w_F]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho*10**(-3)
        #
        data.append(mineral)
        data.append([round(M, 2), x, y])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2)])
        #
        return data
#
################
# COLTAN GROUP #
################
#
class Coltans:
    #
    def __init__(self):
        pass
    #
    def columbite(self):   # (Fe,Mn)Nb2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        manganese = elements.Mn(self)
        niobium = elements.Nb(self)
        oxygen = elements.O(self)
        element = [iron, manganese, niobium, oxygen]
        element_Fe = [iron, niobium, oxygen]
        element_Mn = [manganese, niobium, oxygen]
        #
        data = []
        #
        mineral = "Col"
        #
        # Molar mass
        x_Fe = round(rd.uniform(0.0, 1.0), 2)
        x_Mn = round(1 - x_Fe,2)
        M = round((x_Fe*iron[2] + x_Mn*manganese[2]) + 2*niobium[2] + 6*oxygen[2], 3)
        M_Fe = round(iron[2] + 2*niobium[2] + 6*oxygen[2], 3)
        M_Mn = round(manganese[2] + 2*niobium[2] + 6*oxygen[2], 3)
        w_Fe = round(x_Fe*iron[2]/M, 4)
        w_Mn = round(x_Mn*manganese[2]/M, 4)
        w_Nb = round(2*niobium[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        amounts = [x_Fe, x_Mn, 2, 6]
        amounts_Fe = [1, 2, 6]
        amounts_Mn = [1, 2, 6]
        # Density
        dataV_Fe = CrystalPhysics([[5.746, 14.308, 5.075], [], "orthorhombic"])
        V_Fe = dataV_Fe.calculate_volume()
        dataRho_Fe = CrystalPhysics([M_Fe, 4, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        data_rho_e_Fe = CrystalPhysics([element_Fe, amounts_Fe, rho_Fe])
        rho_e_Fe = data_rho_e_Fe.calculate_electron_density()
        dataV_Mn = CrystalPhysics([[5.767, 14.434, 5.085], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        dataRho_Mn = CrystalPhysics([M_Mn, 4, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        data_rho_e_Mn = CrystalPhysics([element_Mn, amounts_Mn, rho_Mn])
        rho_e_Mn = data_rho_e_Mn.calculate_electron_density()
        dataV = CrystalPhysics([[x_Fe*5.746+x_Mn*5.767, x_Fe*14.308+x_Mn*14.434, x_Fe*5.075+x_Mn*5.085], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K_Fe = 166.16*10**9
        K_Mn = 162.81*10**9
        K = x_Fe*K_Fe + x_Mn*K_Mn
        # Shear modulus
        G_Fe = 75.01*10**9
        G_Mn = 73.84*10**9
        G = x_Fe*G_Fe + x_Mn*G_Mn
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
        weights = [w_Fe, w_Mn, w_Nb, w_O]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 0
        #
        data.append(mineral)
        data.append([round(M, 2), x_Fe, x_Mn])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
    def tantalite(self):   # (Fe,Mn)Ta2O6
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        iron = elements.Fe(self)
        manganese = elements.Mn(self)
        tantalum = elements.Ta(self)
        oxygen = elements.O(self)
        element = [iron, manganese, tantalum, oxygen]
        element_Fe = [iron, tantalum, oxygen]
        element_Mn = [manganese, tantalum, oxygen]
        #
        data = []
        #
        mineral = "Tant"
        #
        # Molar mass
        x_Fe = round(rd.uniform(0.0, 1.0), 2)
        x_Mn = round(1 - x_Fe,2)
        M = round((x_Fe*iron[2] + x_Mn*manganese[2]) + 2*tantalum[2] + 6*oxygen[2], 3)
        M_Fe = round(iron[2] + 2*tantalum[2] + 6*oxygen[2], 3)
        M_Mn = round(manganese[2] + 2*tantalum[2] + 6*oxygen[2], 3)
        w_Fe = round(x_Fe*iron[2]/M, 4)
        w_Mn = round(x_Mn*manganese[2]/M, 4)
        w_Ta = round(2*tantalum[2]/M, 4)
        w_O = round(6*oxygen[2]/M, 4)
        amounts = [x_Fe, x_Mn, 2, 6]
        amounts_Fe = [1, 2, 6]
        amounts_Mn = [1, 2, 6]
        # Density
        dataV_Fe = CrystalPhysics([[5.73, 14.24, 5.08], [], "orthorhombic"])
        V_Fe = dataV_Fe.calculate_volume()
        dataRho_Fe = CrystalPhysics([M_Fe, 4, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        data_rho_e_Fe = CrystalPhysics([element_Fe, amounts_Fe, rho_Fe])
        rho_e_Fe = data_rho_e_Fe.calculate_electron_density()
        dataV_Mn = CrystalPhysics([[5.772, 14.465, 5.097], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        dataRho_Mn = CrystalPhysics([M_Mn, 4, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        data_rho_e_Mn = CrystalPhysics([element_Mn, amounts_Mn, rho_Mn])
        rho_e_Mn = data_rho_e_Mn.calculate_electron_density()
        dataV = CrystalPhysics([[x_Fe*5.73+x_Mn*5.772, x_Fe*14.24+x_Mn*14.465, x_Fe*5.08+x_Mn*5.097], [], "orthorhombic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K_Fe = 199.31*10**9
        K_Mn = 201.87*10**9
        K = x_Fe*K_Fe + x_Mn*K_Mn
        # Shear modulus
        G_Fe = 88.85*10**9
        G_Mn = 94.55*10**9
        G = x_Fe*G_Fe + x_Mn*G_Mn
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
        weights = [w_Fe, w_Mn, w_Ta, w_O]
        PE = bg.calculate_pe(self, x_list=weights, elements_list=element)
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 0
        #
        data.append(mineral)
        data.append([round(M, 2), x_Fe, x_Mn])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data
    #
    def coltan(self):
        # Coltans
        chem_columbite = Coltans.columbite("")
        chem_tantalite = Coltans.tantalite("")
        #
        def lower_limit(x_Mn):
            x = x_Mn*100
            return 2.33*10**(-4)*x**3 - 1.23*10**(-2)*x**2 + 0.5*x + 53
        def upper_limit(x_Mn):
            x = x_Mn*100
            return 1.98*10**(-4)*x**3 - 2.12*10**(-2)*x**2 + 1.11*x + 76
        #
        condition = False
        x_Col = round(rd.uniform(0.0, 1.0), 2)
        x_Tant = round(1 - x_Col,2)
        #
        y_Col_low = lower_limit(chem_columbite[1][2])
        y_Col_high = upper_limit(chem_columbite[1][2])
        #
        while condition == False:
            if x_Col <= y_Col_low or y_Col_high <= x_Col <= 100 or lower_limit(chem_columbite[1][2]) >= 100:
                condition = True
            else:
                condition = False
        #
        data = []
        #
        mineral = "Colt"
        #
        # Molar mass
        M = x_Col*chem_columbite[1][0] + x_Tant*chem_tantalite[1][0]
        # Density
        rho = x_Col*chem_columbite[2] + x_Tant*chem_tantalite[2]
        # Bulk Modulus
        K = (x_Col*chem_columbite[3][0] + x_Tant*chem_tantalite[3][0])*10**9
        # Shear Modulus
        G = (x_Col*chem_columbite[3][1] + x_Tant*chem_tantalite[3][1])*10**9
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
        GR = x_Col*chem_columbite[5][0] + x_Tant*chem_tantalite[5][0]
        # Photoelectricity
        PE = x_Col*chem_columbite[5][1] + x_Tant*chem_tantalite[5][1]
        U = x_Col*chem_columbite[5][2] + x_Tant*chem_tantalite[5][2]
        # Electrical resistivity
        p = x_Col*chem_columbite[5][3] + x_Tant*chem_tantalite[5][3]
        #
        data.append(mineral)
        data.append([round(M, 2), x_Col, x_Tant])
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(GR, 2), round(PE, 2), round(U, 2), p])
        #
        return data