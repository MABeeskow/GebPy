#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		carbonates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		08.09.2022

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals
from modules import fluids
from modules.geophysics import Elasticity as elast
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry
from modules.minerals import Organics
from modules.oxides import Oxides
from modules.silicates import Tectosilicates, Phyllosilicates
from modules.sulfides import Sulfides
from modules.pyllosilicates import Pyllosilicates
from modules.sulfates import Sulfates
from modules.fluids import Water

class limestone:
    #
    def __init__(self, fluid=None, actualThickness=None):
        self.fluid = fluid
        self.actualThickness = actualThickness
    #
    def createLimestone(self):
        # [symbol, atomic number, atomic mass, oxidation states, melting point, boiling point, density, electronegativity]
        chemH = ["H", 1, 1.0078, 1, 13.99, 20.271, 0.084, 2.2]
        chemC = ["C", 6, 12.009, 4, 0.0, 3915, 3510, 2.55]
        chemN = ["N", 7, 14.006, -3, 63.15, 77.355, 1.170, 3.04]
        chemO = ["O", 8, 15.999, -2, 54.3, 90.188, 1.33, 3.44]
        chemS = ["S", 16, 32.059, 6, 368.4, 717.8, 2060, 2.58]
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemDolomite = minerals.carbonates.dolomite("")
        chemCalcite = minerals.carbonates.calcite("")
        chemQuartz = minerals.oxides.quartz("")
        chemAlkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chemPlagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        chemSiderite = minerals.carbonates.siderite("")
        chemKaolinite = minerals.phyllosilicates.kaolinite("")
        chemChlorite = minerals.phyllosilicates.chamosite("")
        chemIllite = minerals.phyllosilicates.illite("")
        #
        # [molar mass, density, bulk modulus, vP]
        chemWater = [2 * chemH[2] + 1 * chemO[2], 997, 2.08, 1444]
        chemOil = [0.83 * chemC[2] + 0.11 * chemH[2] + 0.05 * chemS[2] + 0.005 * (chemO[2] + chemN[2]), 0.9, 1.35, 1225]
        chemGas = [0.83 * chemC[2] + 0.11 * chemH[2] + 0.05 * chemS[2] + 0.005 * (chemO[2] + chemN[2]), 0.8, 0.081, 475]
        #
        limestone = []
        #
        cond = False
        composition = []
        while cond == False:
            magicnumber = rd.randint(0, 2)
            if magicnumber == 0:    # limestone
                xCalcite = rd.randint(60, 100)/100          # Ca-bearing carbonates
                xNonCaCO3 = rd.randint(0, 20)/100           # Non Ca-bearing carbonates
                xDolomite2 = rd.randint(0, 100)/100
                xSiderite2 = 1 - xDolomite2
                xDolomite = xNonCaCO3*xDolomite2
                xSiderite = xNonCaCO3*xSiderite2
                xImpurities = rd.randint(0, 20)/100         # Impurities
                magicnumber2 = rd.randint(0, 2)
                if magicnumber2 == 0:                       # Qz + Fsp
                    xQuartz2 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz2
                    xAlkalifeldspar3 = rd.randint(0, 100)/100
                    xPlagioclase3 = 1 - xAlkalifeldspar3
                    xAlkalifeldspar2 = xFeldspars*xAlkalifeldspar3
                    xPlagioclase2 = xFeldspars*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xKaolinite = 0.0
                    xChlorite = 0.0
                    xIllite = 0.0
                elif magicnumber2 == 1:                     # Clays
                    xKaolinite2 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite2,2)*100)
                    xChlorite2 = round(rd.randint(0, limChl)/100,2)
                    xIllite2 = round(1 - xKaolinite2 - xChlorite2,2)
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
                    xQuartz = 0.0
                    xAlkalifeldspar = 0.0
                    xPlagioclase = 0.0
                elif magicnumber2 == 2:                     # Qz + Fsp + Clays
                    xQzFsp = rd.randint(0, 100)/100
                    xQuartz3 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz3
                    xAlkalifeldspar4 = rd.randint(0, 100)/100
                    xPlagioclase4 = 1 - xAlkalifeldspar4
                    xAlkalifeldspar3 = xFeldspars*xAlkalifeldspar4
                    xPlagioclase3 = xFeldspars*xPlagioclase4
                    xQuartz2 = xQzFsp*xQuartz3
                    xAlkalifeldspar2 = xQzFsp*xAlkalifeldspar3
                    xPlagioclase2 = xQzFsp*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xClays = 1 - xQzFsp
                    xKaolinite3 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite3,2)*100)
                    xChlorite3 = round(rd.randint(0, limChl)/100,2)
                    xIllite3 = round(1 - xKaolinite3 - xChlorite3,2)
                    xKaolinite2 = xClays*xKaolinite3
                    xChlorite2 = xClays*xChlorite3
                    xIllite2 = xClays*xIllite3
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
            elif magicnumber == 1:    # clean limestone
                xCalcite = rd.randint(80, 100)/100          # Ca-bearing carbonates
                xNonCaCO3 = rd.randint(0, 15)/100           # Non Ca-bearing carbonates
                xDolomite2 = rd.randint(0, 100)/100
                xSiderite2 = 1 - xDolomite2
                xDolomite = xNonCaCO3*xDolomite2
                xSiderite = xNonCaCO3*xSiderite2
                xImpurities = rd.randint(0, 5)/100         # Impurities
                magicnumber2 = rd.randint(0, 2)
                if magicnumber2 == 0:                       # Qz + Fsp
                    xQuartz2 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz2
                    xAlkalifeldspar3 = rd.randint(0, 100)/100
                    xPlagioclase3 = 1 - xAlkalifeldspar3
                    xAlkalifeldspar2 = xFeldspars*xAlkalifeldspar3
                    xPlagioclase2 = xFeldspars*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xKaolinite = 0.0
                    xChlorite = 0.0
                    xIllite = 0.0
                elif magicnumber2 == 1:                     # Clays
                    xKaolinite2 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite2,2)*100)
                    xChlorite2 = round(rd.randint(0, limChl)/100,2)
                    xIllite2 = round(1 - xKaolinite2 - xChlorite2,2)
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
                    xQuartz = 0.0
                    xAlkalifeldspar = 0.0
                    xPlagioclase = 0.0
                elif magicnumber2 == 2:                     # Qz + Fsp + Clays
                    xQzFsp = rd.randint(0, 100)/100
                    xQuartz3 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz3
                    xAlkalifeldspar4 = rd.randint(0, 100)/100
                    xPlagioclase4 = 1 - xAlkalifeldspar4
                    xAlkalifeldspar3 = xFeldspars*xAlkalifeldspar4
                    xPlagioclase3 = xFeldspars*xPlagioclase4
                    xQuartz2 = xQzFsp*xQuartz3
                    xAlkalifeldspar2 = xQzFsp*xAlkalifeldspar3
                    xPlagioclase2 = xQzFsp*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xClays = 1 - xQzFsp
                    xKaolinite3 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite3,2)*100)
                    xChlorite3 = round(rd.randint(0, limChl)/100,2)
                    xIllite3 = round(1 - xKaolinite3 - xChlorite3,2)
                    xKaolinite2 = xClays*xKaolinite3
                    xChlorite2 = xClays*xChlorite3
                    xIllite2 = xClays*xIllite3
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
            elif magicnumber == 2:    # dirty limestone
                xCalcite = rd.randint(50, 80)/100          # Ca-bearing carbonates
                xNonCaCO3 = rd.randint(0, 10)/100           # Non Ca-bearing carbonates
                xDolomite2 = rd.randint(0, 100)/100
                xSiderite2 = 1 - xDolomite2
                xDolomite = xNonCaCO3*xDolomite2
                xSiderite = xNonCaCO3*xSiderite2
                xImpurities = rd.randint(20, 40)/100        # Impurities
                magicnumber2 = rd.randint(0, 2)
                if magicnumber2 == 0:                       # Qz + Fsp
                    xQuartz2 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz2
                    xAlkalifeldspar3 = rd.randint(0, 100)/100
                    xPlagioclase3 = 1 - xAlkalifeldspar3
                    xAlkalifeldspar2 = xFeldspars*xAlkalifeldspar3
                    xPlagioclase2 = xFeldspars*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xKaolinite = 0.0
                    xChlorite = 0.0
                    xIllite = 0.0
                elif magicnumber2 == 1:                     # Clays
                    xKaolinite2 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite2,2)*100)
                    xChlorite2 = round(rd.randint(0, limChl)/100,2)
                    xIllite2 = round(1 - xKaolinite2 - xChlorite2,2)
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
                    xQuartz = 0.0
                    xAlkalifeldspar = 0.0
                    xPlagioclase = 0.0
                elif magicnumber2 == 2:                     # Qz + Fsp + Clays
                    xQzFsp = rd.randint(0, 100)/100
                    xQuartz3 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz3
                    xAlkalifeldspar4 = rd.randint(0, 100)/100
                    xPlagioclase4 = 1 - xAlkalifeldspar4
                    xAlkalifeldspar3 = xFeldspars*xAlkalifeldspar4
                    xPlagioclase3 = xFeldspars*xPlagioclase4
                    xQuartz2 = xQzFsp*xQuartz3
                    xAlkalifeldspar2 = xQzFsp*xAlkalifeldspar3
                    xPlagioclase2 = xQzFsp*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xClays = 1 - xQzFsp
                    xKaolinite3 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite3,2)*100)
                    xChlorite3 = round(rd.randint(0, limChl)/100,2)
                    xIllite3 = round(1 - xKaolinite3 - xChlorite3,2)
                    xKaolinite2 = xClays*xKaolinite3
                    xChlorite2 = xClays*xChlorite3
                    xIllite2 = xClays*xIllite3
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
            sumMin = round(xCalcite,2) + round(xDolomite,2) + round(xSiderite,2) + round(xQuartz,2) + round(xAlkalifeldspar,2) + round(xPlagioclase,2) + round(xKaolinite,2) + round(xChlorite,2) + round(xIllite,2)
            if sumMin == 1:
                cond = True
                composition.extend([["Cal", round(xCalcite,2), round(chemCalcite[1],2)], ["Dol", round(xDolomite,2), round(chemDolomite[1],2)], ["Sd", round(xSiderite,2), round(chemSiderite[1],2)], ["Qz", round(xQuartz,2), round(chemQuartz[1],2)], ["Kfs", round(xAlkalifeldspar,2), round(chemAlkalifeldspar[1][0],2), round(chemAlkalifeldspar[1][1],2)], ["Pl", round(xPlagioclase,2), round(chemPlagioclase[1][0],2), round(chemPlagioclase[1][1],2)], ["Kln", round(xKaolinite,2), round(chemKaolinite[1],2)], ["Chl", round(xChlorite,2), round(chemChlorite[1],2)], ["Ilt", round(xIllite,2), round(chemIllite[1],2)]])
            else:
                cond = False
        xCalcite = composition[0][1]
        xDolomite = composition[1][1]
        xSiderite = composition[2][1]
        xQuartz = composition[3][1]
        xAlkalifeldspar = composition[4][1]
        xPlagioclase = composition[5][1]
        xKaolinite = composition[6][1]
        xChlorite = composition[7][1]
        xIllite = composition[8][1]
        limestone.append(composition)
        mineralogy = [chemCalcite, chemDolomite, chemSiderite, chemQuartz, chemAlkalifeldspar, chemPlagioclase, chemKaolinite, chemChlorite, chemIllite]
        #
        rhoSolid = 1.1*(xCalcite*chemCalcite[2] + xDolomite*chemDolomite[2] + xSiderite*chemSiderite[2] + xQuartz*chemQuartz[2] + xAlkalifeldspar*chemAlkalifeldspar[2] + xPlagioclase*chemPlagioclase[2] + xKaolinite*chemKaolinite[2] + xChlorite*chemChlorite[2] + xIllite*chemIllite[2]) / 1000
        X = [xCalcite, xDolomite, xSiderite, xQuartz, xAlkalifeldspar, xPlagioclase, xKaolinite, xChlorite, xIllite]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_arith = elast.calc_arithmetic_mean(self, X, K_list)
        G_arith = elast.calc_arithmetic_mean(self, X, G_list)
        K_solid = K_arith
        G_solid = G_arith
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if self.actualThickness <= 1000:
            phi = randint(30, 35)/100
        elif self.actualThickness > 1000 and self.actualThickness <= 2000:
            phi = randint(25, 30)/100
        elif self.actualThickness > 2000 and self.actualThickness <= 3000:
            phi = randint(15, 25)/100
        elif self.actualThickness > 3000 and self.actualThickness <= 4000:
            phi = randint(10, 15)/100
        elif self.actualThickness > 4000:
            phi = randint(5, 10)/100
        #
        if self.fluid == "water":
            rho = ((1 - phi) * 10/11*rhoSolid + phi * chemWater[1] / 1000)
            vP = (1-phi)*vP_solid + phi*chemWater[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = ((10/11*rhoSolid - rho) / (10/11*rhoSolid - chemWater[1] / 1000))
            phiN = ((2 * phi ** 2 - phiD ** 2) ** (0.5))
            GR = xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xSiderite*chemSiderite[5][0] + xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xKaolinite*chemKaolinite[5][0] + xChlorite*chemChlorite[5][0] + xIllite*chemIllite[5][0]
            PE = xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xSiderite*chemSiderite[5][1] + xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xKaolinite*chemKaolinite[5][1] + xChlorite*chemChlorite[5][1] + xIllite*chemIllite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xCalcite*chemCalcite[3][3] + xDolomite*chemDolomite[3][3] + xSiderite*chemSiderite[3][3] + xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xKaolinite*chemKaolinite[3][3] + xChlorite*chemChlorite[3][3] + xIllite*chemIllite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            limestone.append([round(rho, 3), round(10/11*rhoSolid, 3), round(chemWater[1] / 1000, 6)])
            limestone.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            limestone.append([round(vP, 2), round(vS, 2), round(10/12.5*vP_solid, 2), round(chemWater[3], 2)])
            limestone.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            limestone.append("water")
            limestone.append([GR, PE])
        elif self.fluid == "oil":
            rho = (1 - phi) * rhoSolid + phi * chemOil[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemOil[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemOil[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xSiderite*chemSiderite[5][0] + xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xKaolinite*chemKaolinite[5][0] + xChlorite*chemChlorite[5][0] + xIllite*chemIllite[5][0]
            PE = xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xSiderite*chemSiderite[5][1] + xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xKaolinite*chemKaolinite[5][1] + xChlorite*chemChlorite[5][1] + xIllite*chemIllite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xCalcite*chemCalcite[3][3] + xDolomite*chemDolomite[3][3] + xSiderite*chemSiderite[3][3] + xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xKaolinite*chemKaolinite[3][3] + xChlorite*chemChlorite[3][3] + xIllite*chemIllite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            limestone.append([round(rho, 3), round(rhoSolid, 3), round(chemOil[1] / 1000, 6)])
            limestone.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            limestone.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemOil[3], 2)])
            limestone.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            limestone.append("oil")
            limestone.append([GR, PE])
        elif self.fluid == "gas":
            rho = (1 - phi) * rhoSolid + phi * chemGas[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemGas[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemGas[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xSiderite*chemSiderite[5][0] + xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xKaolinite*chemKaolinite[5][0] + xChlorite*chemChlorite[5][0] + xIllite*chemIllite[5][0]
            PE = xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xSiderite*chemSiderite[5][1] + xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xKaolinite*chemKaolinite[5][1] + xChlorite*chemChlorite[5][1] + xIllite*chemIllite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xCalcite*chemCalcite[3][3] + xDolomite*chemDolomite[3][3] + xSiderite*chemSiderite[3][3] + xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xKaolinite*chemKaolinite[3][3] + xChlorite*chemChlorite[3][3] + xIllite*chemIllite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            limestone.append([round(rho, 3), round(rhoSolid, 3), round(chemGas[1] / 1000, 6)])
            limestone.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            limestone.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemGas[3], 2)])
            limestone.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            limestone.append("gas")
            limestone.append([GR, PE])
        #
        #  limestone = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return limestone
    #
    def create_simple_limestone(self, w_Na=None, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None, porosity=None, dict=False):
        #
        results = {}
        results["rock"] = "limestone"
        #
        self.w_Na = w_Na
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # mineralogy
        chem_cal = minerals.carbonates.calcite("")
        chem_arg = minerals.carbonates.aragonite("")
        chem_dol = minerals.carbonates.dolomite("")
        chem_sd = minerals.carbonates.siderite("")
        chem_qz = minerals.oxides.quartz("")
        chem_kfs = minerals.feldspars.alkalifeldspar(self, "Kfs")
        chem_pl = minerals.feldspars.plagioclase(self, "Pl")
        chem_mnt = minerals.phyllosilicates.montmorillonite("")
        chem_kln = minerals.phyllosilicates.kaolinite("")
        chem_chl = minerals.phyllosilicates.chamosite("")
        chem_py = minerals.sulfides.pyrite("")
        #
        mineralogy = [chem_cal, chem_arg, chem_dol, chem_sd, chem_qz, chem_kfs, chem_pl, chem_mnt, chem_kln, chem_chl, chem_py]
        #for i in range(len(mineralogy)):
        #    print(mineralogy[i][0], mineralogy[i][6])
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        oil = fluids.Hydrocarbons.oil("")
        gas = fluids.Hydrocarbons.natural_gas("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Na == None and self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                magicnumber = rd.randint(0, 6)
                if magicnumber == 0:    # Cal-rich
                    w_carb = round(rd.uniform(0.9, 1.0), 4)
                    w_cal2 = rd.uniform(0.85, 0.95)
                    w_arg2 = rd.randint(0, int((1-w_cal2)*100))/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 1:    # Arg-rich
                    w_carb = round(rd.uniform(0.90, 1.0), 4)
                    w_arg2 = rd.uniform(0.85, 0.95)
                    w_cal2 = rd.randint(0, int((1-w_arg2)*100))/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 2:    # Dol-rich
                    w_carb = round(rd.uniform(0.90, 1.0), 4)
                    w_cal2 = rd.randint(50, 75)/100
                    w_dol2 = rd.randint(25, int((1-w_cal2)*100))/100
                    w_arg2 = rd.randint(0, int((1-w_cal2-w_dol2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 3:    # Mixed
                    w_carb = round(rd.uniform(0.60, 0.75), 4)
                    w_cal2 = rd.randint(40, 60)/100
                    w_arg2 = rd.randint(0, 10)/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 4:    # Qz-rich
                    w_carb = round(rd.uniform(0.60, 0.75), 4)
                    w_cal2 = rd.randint(40, 60)/100
                    w_arg2 = rd.randint(0, 10)/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_qz = round(rd.uniform(0.1, 0.2), 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 5:    # Fsp-rich
                    w_carb = round(rd.uniform(0.60, 0.75), 4)
                    w_cal2 = rd.randint(40, 60)/100
                    w_arg2 = rd.randint(0, 10)/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_fsp = round(rd.uniform(0.1, 0.2), 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb-w_fsp)*100))/100, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 6:    # Clay-rich
                    w_carb = round(rd.uniform(0.60, 0.75), 4)
                    w_cal2 = rd.randint(40, 60)/100
                    w_arg2 = rd.randint(0, 10)/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_clay = round(rd.uniform(0.1, 0.2), 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb-w_clay)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_clay-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Na != None:
                condition = False
                while condition == False:
                    condition_1 = False
                    while condition_1 == False:
                        chem_kfs = minerals.feldspars.alkalifeldspar(self, "Na")
                        chem_pl = minerals.feldspars.plagioclase(self, "Na")
                        w_fsp = round(rd.randint(5, 15)/100, 4)
                        w_kfs2 = rd.randint(0, 100)/100
                        w_pl2 = 1-w_kfs2
                        w_kfs = round(w_fsp*w_kfs2, 4)
                        w_pl = round(w_fsp*w_pl2, 4)
                        w_mnt = round((self.w_Na - chem_kfs[6][1]*w_kfs - chem_pl[6][1]*w_pl)/(chem_mnt[6][2]), 4)
                        if w_mnt >= 0.0:
                            condition_1 = True
                        else:
                            condition_1 = False
                    w_clay = round(rd.randint(int(w_mnt*100+1), int(w_mnt*100+10))/100, 4)
                    w_mnt2 = w_mnt/w_clay
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    if w_fsp + w_clay <= 1.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_mnt >= 0.0 and w_kln >= 0.0 and w_chl >= 0.0:
                        condition = True
                    else:
                        condition = False
                w_carb = round(rd.randint(int(0.75*(1-w_fsp-w_clay)*100), int((1-w_fsp-w_clay)*100))/100, 4)
                w_cal2 = rd.randint(40, 60)/100
                w_arg2 = rd.randint(0, 10)/100
                w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                w_sd2 = 1-w_cal2-w_arg2-w_dol2
                w_cal = round(w_carb*w_cal2, 4)
                w_arg = round(w_carb*w_arg2, 4)
                w_dol = round(w_carb*w_dol2, 4)
                w_sd = round(w_carb*w_sd2, 4)
                w_qz = round(rd.randint(0, int((1-w_carb-w_fsp-w_clay)*100))/100, 4)
                w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Mg != None:
                condition = False
                while condition == False:
                    w_clay = round(rd.randint(0, 5)/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_dol = round((self.w_Mg - chem_mnt[6][3]*w_mnt - chem_chl[6][4]*w_chl)/(chem_dol[6][2]), 4)
                    w_carb = round(rd.randint(int(0.75*(1-w_dol-w_clay)*100), int(0.9*(1-w_dol-w_clay)*100))/100, 4)
                    w_dol2 = w_dol/w_carb
                    w_cal2 = rd.randint(int(0.75*(1-w_dol2)*100), int((1-w_dol2)*100))/100
                    w_arg2 = rd.randint(0, int((1-w_dol2-w_cal2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    if w_clay + w_carb <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_qz = round(rd.randint(0, int((1-w_carb-w_clay)*100))/100, 4)
                w_fsp = round(rd.randint(0, int((1-w_carb-w_qz-w_clay)*100))/100, 4)
                w_kfs2 = rd.randint(0, 100)/100
                w_pl2 = 1-w_kfs2
                w_kfs = round(w_fsp*w_kfs2, 4)
                w_pl = round(w_fsp*w_pl2, 4)
                w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_K != None:
                condition = False
                while condition == False:
                    w_kfs = round(self.w_K/chem_kfs[6][4], 4)
                    w_fsp = round(rd.randint(0, 20)/100, 4)
                    w_kfs2 = w_kfs/w_fsp
                    w_pl2 = 1-w_kfs2
                    w_pl = round(w_fsp*w_pl2, 4)
                    if w_fsp <= 1.0 and w_kfs2 <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_carb = round(rd.randint(60, 75)/100, 4)
                w_cal2 = rd.randint(40, 60)/100
                w_arg2 = rd.randint(0, 10)/100
                w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                w_sd2 = 1-w_cal2-w_arg2-w_dol2
                w_cal = round(w_carb*w_cal2, 4)
                w_arg = round(w_carb*w_arg2, 4)
                w_dol = round(w_carb*w_dol2, 4)
                w_sd = round(w_carb*w_sd2, 4)
                w_qz = round(rd.randint(0, int((1-w_carb-w_fsp)*100))/100, 4)
                w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                w_mnt2 = rd.randint(0, 100)/100
                w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                w_chl2 = 1-w_mnt2-w_kln2
                w_mnt = round(w_clay*w_mnt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_chl = round(w_clay*w_chl2, 4)
                w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Ca != None:
                condition = False
                while condition == False:
                    w_fsp = round(rd.randint(0, 5)/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, 5)/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_carb = round(rd.randint(60, 90)/100, 4)
                    w_arg2 = rd.randint(0, 10)/100
                    w_dol2 = rd.randint(0, 20)/100
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_cal = (self.w_Ca - chem_dol[6][3]*w_dol - chem_arg[6][2]*w_arg - chem_pl[6][4]*w_pl - chem_mnt[6][6]*w_mnt)/(chem_cal[6][2])
                    w_cal2 = w_cal/w_carb
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_sd = round(w_carb*w_sd2, 4)
                    if w_fsp + w_clay + w_carb <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_qz = round(rd.randint(0, int((1-w_carb-w_fsp-w_clay)*100))/100, 4)
                w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Fe != None:
                condition = False
                while condition == False:
                    w_carb = round(rd.randint(60, 75)/100, 4)
                    w_cal2 = rd.randint(40, 60)/100
                    w_arg2 = rd.randint(0, 10)/100
                    w_dol2 = rd.randint(0, int((1-w_cal2-w_arg2)*100))/100
                    w_sd2 = 1-w_cal2-w_arg2-w_dol2
                    w_cal = round(w_carb*w_cal2, 4)
                    w_arg = round(w_carb*w_arg2, 4)
                    w_dol = round(w_carb*w_dol2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_clay = round(rd.randint(0, 10)/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_chl2 = 1-w_mnt2-w_kln2
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_chl = round(w_clay*w_chl2, 4)
                    w_py = (self.w_Fe - chem_sd[6][2]*w_sd - chem_chl[6][5]*w_chl)/(chem_py[6][1])
                    if w_carb + w_clay + w_py <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_qz = round(rd.randint(0, int((1-w_carb-w_clay-w_py)*100))/100, 4)
                w_fsp = round(1-w_carb-w_clay-w_py-w_qz, 4)
                w_kfs2 = rd.randint(0, 100)/100
                w_pl2 = 1-w_kfs2
                w_kfs = round(w_fsp*w_kfs2, 4)
                w_pl = round(w_fsp*w_pl2, 4)
            elif type(self.amounts) is list:
                w_cal = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_arg = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_dol = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_sd = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_qz = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_mnt = round(abs(np.random.normal(self.amounts[7], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[8], 0.025)), 4)
                w_chl = round(abs(np.random.normal(self.amounts[9], 0.025)), 4)
                w_py = round(1-w_cal-w_arg-w_dol-w_sd-w_qz-w_kfs-w_pl-w_mnt-w_kln-w_chl, 4)
            #
            if w_cal >= 0.0 and w_arg >= 0.0 and w_dol >= 0.0 and w_sd >= 0.0 and w_qz >= 0.0 and w_kfs  >= 0.0 and w_pl >= 0.0 and w_mnt >= 0.0 and w_kln >= 0.0 and w_chl >= 0.0 and w_py >= 0.0:
                sumMin = round(w_cal + w_arg + w_dol + w_sd + w_qz + w_kfs + w_pl + w_mnt + w_kln + w_chl + w_py, 4)
            else:
                sumMin = 0
            #
            w_H = round(chem_mnt[6][0]*w_mnt + chem_kln[6][0]*w_kln + chem_chl[6][0]*w_chl, 4)
            w_C = round(chem_cal[6][0]*w_cal + chem_arg[6][0]*w_arg + chem_dol[6][0]*w_dol + chem_sd[6][0]*w_sd, 4)
            #w_O = round(chem_cal[6][1]*w_cal + chem_arg[6][1]*w_arg + chem_dol[6][1]*w_dol + chem_sd[6][1]*w_sd + chem_qz[6][1]*w_qz + chem_kfs[6][0]*w_kfs + chem_pl[6][0]*w_pl + chem_mnt[6][1]*w_mnt + chem_kln[6][1]*w_kln + chem_chl[6][1]*w_chl, 4)
            w_Na = round(chem_kfs[6][1]*w_kfs + chem_pl[6][1]*w_pl + chem_mnt[6][2]*w_mnt, 4)
            w_Mg = round(chem_dol[6][2]*w_dol + chem_mnt[6][3]*w_mnt + chem_chl[6][4]*w_chl, 4)
            w_Al = round(chem_kfs[6][2]*w_kfs + chem_pl[6][2]*w_pl + chem_mnt[6][4]*w_mnt + chem_kln[6][2]*w_kln + chem_chl[6][2]*w_chl, 4)
            w_Si = round(chem_qz[6][1]*w_qz + chem_kfs[6][3]*w_kfs + chem_pl[6][3]*w_pl + chem_mnt[6][5]*w_mnt + chem_kln[6][3]*w_kln + chem_chl[6][3]*w_chl, 4)
            w_S = round(chem_py[6][0]*w_py, 4)
            w_K = round(chem_kfs[6][4]*w_kfs, 4)
            w_Ca = round(chem_cal[6][2]*w_cal + chem_arg[6][2]*w_arg + chem_dol[6][3]*w_dol + chem_pl[6][4]*w_pl + chem_mnt[6][6]*w_mnt, 4)
            w_Fe = round(chem_sd[6][2]*w_sd + chem_py[6][1]*w_py + chem_chl[6][5]*w_chl, 4)
            w_O = round(1 - w_H - w_C - w_Na - w_Mg - w_Al - w_Si - w_S - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_C + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                #composition.extend((["Cal", w_cal, round(chem_cal[1], 2)], ["Arg", w_arg, round(chem_arg[1], 2)], ["Dol", w_dol, round(chem_dol[1], 2)], ["Sd", w_sd, round(chem_sd[1], 2)], ["Qz", w_qz, round(chem_qz[1], 2)], ["Kfs", w_kfs, round(chem_kfs[1][0], 2), round(chem_kfs[1][1], 2)], ["Pl", w_pl, round(chem_pl[1][0], 2), round(chem_pl[1][1], 2)], ["Mnt", w_mnt, round(chem_mnt[1][0], 2), round(chem_mnt[1][1], 2)], ["Kln", w_kln, round(chem_kln[1], 2)], ["Chl", w_chl, round(chem_chl[1], 2)], ["Py", w_py, round(chem_py[1], 2)]))
                composition.extend((["Cal", "Arg", "Dol", "Sd", "Qz", "Kfs", "Pl", "Mnt", "Kln", "Chl", "Py"]))
                concentrations = [w_H, w_C, w_O, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe]
                amounts = [w_cal, w_arg, w_dol, w_sd, w_qz, w_kfs, w_pl, w_mnt, w_kln, w_chl, w_py]
            else:
                cond = False
        #
        element_list = ["H", "C", "O", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Fe"]
        mineral_list = ["Cal", "Arg", "Dol", "Sd", "Qz", "Kfs", "Pl", "Mnt", "Kln", "Chl", "Py"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_cal*chem_cal[2] + w_arg*chem_arg[2] + w_dol*chem_dol[2] + w_sd*chem_sd[2] + w_qz*chem_qz[2] + w_kfs*chem_kfs[2] + w_pl*chem_pl[2] + w_mnt*chem_mnt[2] + w_kln*chem_kln[2] + w_chl*chem_chl[2] + w_py*chem_py[2]) / 1000
        X = [w_cal, w_arg, w_dol, w_sd, w_qz, w_kfs, w_pl, w_mnt, w_kln, w_chl, w_py]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.randint(35, 40)/100
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.randint(30, 35)/100
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.randint(20, 30)/100
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.randint(10, 20)/100
            elif self.actualThickness > 4000:
                phi = rd.randint(5, 10)/100
        else:
            phi = porosity
        #
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        if self.fluid == "water":
            rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
            vP = (1-phi)*vP_solid + phi*water[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_cal*chem_cal[5][0] + w_arg*chem_arg[5][0] + w_dol*chem_dol[5][0] + w_sd*chem_sd[5][0] + w_qz*chem_qz[5][0] + w_kfs*chem_kfs[5][0] + w_pl*chem_pl[5][0] + w_mnt*chem_mnt[5][0] + w_kln*chem_kln[5][0] + w_chl*chem_chl[5][0] + w_py*chem_py[5][0]
            PE = w_cal*chem_cal[5][1] + w_arg*chem_arg[5][1] + w_dol*chem_dol[5][1] + w_sd*chem_sd[5][1] + w_qz*chem_qz[5][1] + w_kfs*chem_kfs[5][1] + w_pl*chem_pl[5][1] + w_mnt*chem_mnt[5][1] + w_kln*chem_kln[5][1] + w_chl*chem_chl[5][1] + w_py*chem_py[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_cal*chem_cal[3][3] + w_arg*chem_arg[3][3] + w_dol*chem_dol[3][3] + w_sd*chem_sd[3][3] + w_qz*chem_qz[3][3] + w_kfs*chem_kfs[3][3] + w_pl*chem_pl[3][3] + w_mnt*chem_mnt[3][3] + w_kln*chem_kln[3][3] + w_chl*chem_chl[3][3] + w_py*chem_py[3][3]
            #
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("water")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        elif self.fluid == "oil":
            rho = (1 - phi) * rhoSolid + phi * oil[2] / 1000
            vP = (1-phi)*vP_solid + phi*oil[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - oil[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_cal*chem_cal[5][0] + w_arg*chem_arg[5][0] + w_dol*chem_dol[5][0] + w_sd*chem_sd[5][0] + w_qz*chem_qz[5][0] + w_kfs*chem_kfs[5][0] + w_pl*chem_pl[5][0] + w_mnt*chem_mnt[5][0] + w_kln*chem_kln[5][0] + w_chl*chem_chl[5][0] + w_py*chem_py[5][0]
            PE = w_cal*chem_cal[5][1] + w_arg*chem_arg[5][1] + w_dol*chem_dol[5][1] + w_sd*chem_sd[5][1] + w_qz*chem_qz[5][1] + w_kfs*chem_kfs[5][1] + w_pl*chem_pl[5][1] + w_mnt*chem_mnt[5][1] + w_kln*chem_kln[5][1] + w_chl*chem_chl[5][1] + w_py*chem_py[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_cal*chem_cal[3][3] + w_arg*chem_arg[3][3] + w_dol*chem_dol[3][3] + w_sd*chem_sd[3][3] + w_qz*chem_qz[3][3] + w_kfs*chem_kfs[3][3] + w_pl*chem_pl[3][3] + w_mnt*chem_mnt[3][3] + w_kln*chem_kln[3][3] + w_chl*chem_chl[3][3] + w_py*chem_py[3][3]
            #
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(oil[2] / 1000, 6)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(oil[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("oil")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        elif self.fluid == "gas":
            rho = (1 - phi) * rhoSolid + phi * gas[2] / 1000
            vP = (1-phi)*vP_solid + phi*gas[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - gas[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_cal*chem_cal[5][0] + w_arg*chem_arg[5][0] + w_dol*chem_dol[5][0] + w_sd*chem_sd[5][0] + w_qz*chem_qz[5][0] + w_kfs*chem_kfs[5][0] + w_pl*chem_pl[5][0] + w_mnt*chem_mnt[5][0] + w_kln*chem_kln[5][0] + w_chl*chem_chl[5][0] + w_py*chem_py[5][0]
            PE = w_cal*chem_cal[5][1] + w_arg*chem_arg[5][1] + w_dol*chem_dol[5][1] + w_sd*chem_sd[5][1] + w_qz*chem_qz[5][1] + w_kfs*chem_kfs[5][1] + w_pl*chem_pl[5][1] + w_mnt*chem_mnt[5][1] + w_kln*chem_kln[5][1] + w_chl*chem_chl[5][1] + w_py*chem_py[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_cal*chem_cal[3][3] + w_arg*chem_arg[3][3] + w_dol*chem_dol[3][3] + w_sd*chem_sd[3][3] + w_qz*chem_qz[3][3] + w_kfs*chem_kfs[3][3] + w_pl*chem_pl[3][3] + w_mnt*chem_mnt[3][3] + w_kln*chem_kln[3][3] + w_chl*chem_chl[3][3] + w_py*chem_py[3][3]
            #
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(gas[2] / 1000, 6)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(gas[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("gas")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        #
        if dict == False:
            return data
        else:
            return results
    #
class CarbonateRocks:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        #
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_aragonite = Carbonates(impurity="pure", data_type=True).create_aragonite()
        self.data_dolomite = Carbonates(impurity="pure", data_type=True).create_dolomite()
        self.data_siderite = Carbonates(impurity="pure", data_type=True).create_siderite()
        self.data_magnesite = Carbonates(impurity="pure", data_type=True).create_magnesite()
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        self.data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        self.data_water = Water.water("")
    #
    def create_limestone(self, number, porosity=None, dominance="Cal", dominant_group="Carbonates"):
        #
        data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
        data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Na")
        #
        assemblage = [self.data_calcite, self.data_aragonite, self.data_dolomite, self.data_siderite,
                      self.data_magnesite,  self.data_quartz, data_illite, self.data_kaolinite, data_montmorillonite,
                      data_alkalifeldspar, data_plagioclase]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage:
            amounts_mineralogy[mineral["mineral"]] = []
            mineral_list.append(mineral["mineral"])
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        mineral_list_carb = ["Arg", "Cal", "Dol", "Mgs", "Sd"]
        mineral_list_silic = ["Kfs", "Pl", "Qz"]
        mineral_list_clays = ["Ilt", "Kln", "Mnt"]
        #
        n = 0
        amounts_helper = []
        while n < number:
            if dominant_group == "Carbonates":
                w_carb = round(rd.uniform(0.8, 1.0), 4)
                w_silic = round(rd.uniform(0.0, (1 - w_carb)), 4)
                w_clay = round(1 - w_carb - w_silic, 4)
            elif dominant_group == "Siliciclastics":
                w_silic = round(rd.uniform(0.2, 0.3), 4)
                w_clay = round(rd.uniform(0.0, 0.05), 4)
                w_carb = round(1 - w_silic - w_clay, 4)
            elif dominant_group == "Clays":
                w_clay = round(rd.uniform(0.2, 0.3), 4)
                w_silic = round(rd.uniform(0.0, 0.05), 4)
                w_carb = round(1 - w_clay - w_silic, 4)
            mineral_fractions = {"Carbonates": w_carb, "Siliciclastics": w_silic, "Clays": w_clay}
            w_total = 0
            n_minerals = 0
            #
            for group, fraction in mineral_fractions.items():
                if group == "Carbonates":
                    condition_carb = False
                    while condition_carb == False:
                        if dominance == "Cal":
                            phi_arg = round(fraction * rd.uniform(0.0, 0.05), 4)
                            phi_cal = round(fraction*rd.uniform(0.8, (1.0 - phi_arg)), 4)
                        elif dominance == "Arg":
                            phi_cal = round(fraction * rd.uniform(0.0, 0.05), 4)
                            phi_arg = round(fraction * rd.uniform(0.8, (1.0 - phi_cal)), 4)
                        #
                        phi_dol = round(fraction * rd.uniform(0.0, (1.0 - phi_arg - phi_cal)), 4)
                        #
                        magicnumber = rd.randint(0, 1)
                        #
                        if magicnumber == 0:
                            phi_mgs = round(fraction * rd.uniform(0.0, (1.0 - phi_arg - phi_cal - phi_dol)), 4)
                            phi_sd = round(fraction - phi_arg - phi_cal - phi_dol - phi_mgs, 4)
                        elif magicnumber == 1:
                            phi_sd = round(fraction * rd.uniform(0.0, (1.0 - phi_arg - phi_cal - phi_dol)), 4)
                            phi_mgs = round(fraction - phi_arg - phi_cal - phi_dol - phi_sd, 4)
                        #
                        phi_carb = phi_arg + phi_cal + phi_dol + phi_mgs + phi_sd
                        #
                        if np.isclose(phi_carb, fraction) == True:
                            phi_list = [phi_arg, phi_cal, phi_dol, phi_mgs, phi_sd]
                            if all(x >= 0 for x in phi_list):
                                amounts_helper.extend(phi_list)
                                n_minerals += len(mineral_list_carb)
                                w_total += fraction
                                condition_carb = True
                #
                elif group == "Siliciclastics":
                    condition_silic = False
                    while condition_silic == False:
                        phi_qz = round(fraction * rd.uniform(0.0, 1.0), 4)
                        #
                        magicnumber = rd.randint(0, 1)
                        #
                        if magicnumber == 0:
                            phi_kfs = round(fraction * rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_pl = round(fraction - phi_qz - phi_kfs, 4)
                        elif magicnumber == 1:
                            phi_pl = round(fraction * rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_kfs = round(fraction - phi_qz - phi_pl, 4)
                        #
                        phi_silic = phi_qz + phi_kfs + phi_pl
                        #
                        if np.isclose(phi_silic, fraction) == True:
                            phi_list = [phi_kfs, phi_pl, phi_qz]
                            if all(x >= 0 for x in phi_list):
                                amounts_helper.extend(phi_list)
                                n_minerals += len(mineral_list_silic)
                                w_total += fraction
                                condition_silic = True
                #
                elif group == "Clays":
                    condition_clays = False
                    while condition_clays == False:
                        magicnumber = rd.randint(0, 2)
                        if magicnumber == 0:
                            phi_ilt = round(fraction * rd.uniform(0.0, 1.0), 4)
                            phi_kln = round(fraction * rd.uniform(0.0, (1.0 - phi_ilt)), 4)
                            phi_mnt = round(fraction - phi_ilt - phi_kln, 4)
                        elif magicnumber == 1:
                            phi_kln = round(fraction * rd.uniform(0.0, 1.0), 4)
                            phi_ilt = round(fraction * rd.uniform(0.0, (1.0 - phi_kln)), 4)
                            phi_mnt = round(fraction - phi_ilt - phi_kln, 4)
                        elif magicnumber == 2:
                            phi_mnt = round(fraction * rd.uniform(0.0, 1.0), 4)
                            phi_ilt = round(fraction * rd.uniform(0.0, (1.0 - phi_mnt)), 4)
                            phi_kln = round(fraction - phi_ilt - phi_mnt, 4)
                        #
                        phi_clays = phi_ilt + phi_kln + phi_mnt
                        #
                        if np.isclose(phi_clays, fraction) == True:
                            phi_list = [phi_ilt, phi_kln, phi_mnt]
                            if all(x >= 0 for x in phi_list):
                                amounts_helper.extend(phi_list)
                                n_minerals += len(mineral_list_clays)
                                w_total += fraction
                                condition_clays = True
            #
            mineral_list.clear()
            mineral_list.extend(mineral_list_carb)
            mineral_list.extend(mineral_list_silic)
            mineral_list.extend(mineral_list_clays)
            if np.sum(amounts_helper) == 1.0 and n_minerals == len(mineral_list):
                for index, mineral in enumerate(mineral_list):
                    amounts_mineralogy[mineral].append(amounts_helper[index])
                n += 1
                amounts_helper.clear()
            else:
                n += 0
                amounts_helper.clear()
        #
        n = 0
        amounts_helper = {}
        while n < number:
            w_total = 0
            n_elements = 0
            rho_s_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.4), 4)
            else:
                phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
            #
            data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
            data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase(enrichment="Na")
            #
            old_index = elements.index("O")
            elements += [elements.pop(old_index)]
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_aragonite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Arg"][n] * self.data_aragonite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_dolomite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Dol"][n] * self.data_dolomite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_siderite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Sd"][n] * self.data_siderite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_magnesite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Mgs"][n] * self.data_magnesite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Qz"][n] * self.data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_illite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Ilt"][n] * data_illite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_alkalifeldspar["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Kfs"][n] * data_alkalifeldspar["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_plagioclase["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Pl"][n] * data_plagioclase["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_montmorillonite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Mnt"][n] * data_montmorillonite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_kaolinite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Kln"][n] * self.data_kaolinite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                #
                n_elements += 1
            #
            shear_factor = 1.0
            #
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Arg":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_aragonite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_aragonite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_aragonite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_aragonite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_aragonite["PE"], 3)
                    elif mineral == "Dol":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_dolomite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_siderite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_siderite["PE"], 3)
                    elif mineral == "Sd":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_siderite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_siderite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_siderite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_siderite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_siderite["PE"], 3)
                    elif mineral == "Mgs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_magnesite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["PE"], 3)
                    elif mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                    elif mineral == "Ilt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_illite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_illite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_illite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_illite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_illite["PE"], 3)
                    elif mineral == "Mnt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_montmorillonite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["PE"], 3)
                    elif mineral == "Kln":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_kaolinite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_kaolinite["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_alkalifeldspar["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["PE"], 3)
                #
                rho_helper = round((1 - phi_helper) * rho_s_helper + phi_helper * self.data_water[2] / 1000, 3)
                youngsmod_helper = round(
                    (9 * bulkmod_helper * shearmod_helper) / (3 * bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round(
                    (3 * bulkmod_helper - 2 * shearmod_helper) / (6 * bulkmod_helper + 2 * shearmod_helper), 3)
                vP_helper = round(
                    ((bulkmod_helper * 10 ** 9 + 4 / 3 * shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vS_helper = round(((shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vPvS_helper_helper = round(vP_helper / vS_helper, 3)
                #
                bulk_properties["rho_s"].append(round(rho_s_helper, 3))
                bulk_properties["rho"].append(rho_helper)
                bulk_properties["K"].append(round(bulkmod_helper, 3))
                bulk_properties["G"].append(round(shearmod_helper, 3))
                bulk_properties["E"].append(youngsmod_helper)
                bulk_properties["nu"].append(poisson_helper)
                bulk_properties["vP"].append(vP_helper)
                bulk_properties["vS"].append(vS_helper)
                bulk_properties["vPvS"].append(vPvS_helper_helper)
                bulk_properties["GR"].append(round(gr_helper, 3))
                bulk_properties["PE"].append(round(pe_helper, 3))
                bulk_properties["phi"].append(round(phi_helper, 3))
                n += 1
        #
        results = {}
        results["rock"] = "Limestone"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = bulk_properties["phi"][0]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"][0]
            results["rho"] = bulk_properties["rho"][0]
            results["vP"] = bulk_properties["vP"][0]
            results["vS"] = bulk_properties["vS"][0]
            results["vP/vS"] = bulk_properties["vPvS"][0]
            results["K"] = bulk_properties["K"][0]
            results["G"] = bulk_properties["G"][0]
            results["E"] = bulk_properties["E"][0]
            results["nu"] = bulk_properties["nu"][0]
            results["GR"] = bulk_properties["GR"][0]
            results["PE"] = bulk_properties["PE"][0]
        #
        return results
    #
    def create_limestone_alternative(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
        data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
        #
        mineralogy = {"Cal": self.data_calcite, "Dol": self.data_dolomite, "Sd": self.data_siderite,
                      "Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase,
                      "Kln": self.data_kaolinite, "Mnt": data_montmorillonite, "Ilt": data_illite}
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral, dataset in mineralogy.items():
            amounts_mineralogy[dataset["mineral"]] = []
            mineral_list.append(dataset["mineral"])
            elements_mineral = list(dataset["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        while n < number:
            condition = False
            while condition == False:
                elements_list = []
                phi_minerals = {}
                w_minerals = {}
                w_elements = {}
                #
                if composition != None:
                    phi_cal = composition["Cal"]
                    phi_dol = composition["Dol"]
                    phi_sd = composition["Sd"]
                    phi_qz = composition["Qz"]
                    phi_kfs = composition["Kfs"]
                    phi_pl = composition["Pl"]
                    phi_kln = composition["Kln"]
                    phi_mnt = composition["Mnt"]
                    phi_ilt = composition["Ilt"]
                    #
                    phi_minerals["Cal"] = phi_cal
                    phi_minerals["Dol"] = phi_dol
                    phi_minerals["Sd"] = phi_sd
                    phi_minerals["Qz"] = phi_qz
                    phi_minerals["Kfs"] = phi_kfs
                    phi_minerals["Pl"] = phi_pl
                    phi_minerals["Kln"] = phi_kln
                    phi_minerals["Mnt"] = phi_mnt
                    phi_minerals["Ilt"] = phi_ilt
                else:
                    condition_2 = False
                    while condition_2 == False:
                        magicnumber = rd.randint(0, 12)
                        if 0 <= magicnumber <= 8:   # Carbonate-dominated
                            w_carb = round(rd.uniform(0.85, 1.0), 4)
                            w_clast = round(rd.uniform(0.0, (1.0 - w_carb)), 4)
                            w_clay = round(1 - w_carb - w_clast, 4)
                            #
                            phi_cal = round(w_carb*rd.uniform(0.9, 1.0), 4)
                            phi_dol = round(w_carb*rd.uniform(0.0, (1.0 - phi_cal)), 4)
                            phi_sd = round(w_carb - phi_cal - phi_dol, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.0, 1.0), 4)
                            phi_kfs = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_pl = round(w_clast - phi_qz - phi_kfs, 4)
                            #
                            phi_kln = round(w_clay*rd.uniform(0.0, 1.0), 4)
                            phi_mnt = round(w_clay*rd.uniform(0.0, (1.0 - phi_kln)), 4)
                            phi_ilt = round(1 - phi_cal - phi_dol - phi_sd - phi_qz - phi_kfs - phi_pl - phi_kln - phi_mnt, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_sd + phi_qz + phi_kfs + phi_pl + phi_kln + phi_mnt + phi_ilt
                        elif magicnumber in [9, 10]:   # Clastic-dominated
                            w_clast = round(rd.uniform(0.1, 0.25), 4)
                            w_carb = round(rd.uniform(0.7, (1.0 - w_clast)), 4)
                            w_clay = round(1 - w_carb - w_clast, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.0, 1.0), 4)
                            phi_kfs = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_pl = round(w_clast - phi_qz - phi_kfs, 4)
                            #
                            phi_cal = round(w_carb*rd.uniform(0.8, 1.0), 4)
                            phi_dol = round(w_carb*rd.uniform(0.0, (1.0 - phi_cal)), 4)
                            phi_sd = round(w_carb - phi_cal - phi_dol, 4)
                            #
                            phi_kln = round(w_clay*rd.uniform(0.0, 1.0), 4)
                            phi_mnt = round(w_clay*rd.uniform(0.0, (1.0 - phi_kln)), 4)
                            phi_ilt = round(1 - phi_cal - phi_dol - phi_sd - phi_qz - phi_kfs - phi_pl - phi_kln - phi_mnt, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_sd + phi_qz + phi_kfs + phi_pl + phi_kln + phi_mnt + phi_ilt
                        elif magicnumber in [11, 12]:   # Clay-dominated
                            w_clay = round(rd.uniform(0.1, 0.25), 4)
                            w_carb = round(rd.uniform(0.7, (1.0 - w_clay)), 4)
                            w_clast = round(1 - w_carb - w_clay, 4)
                            #
                            phi_cal = round(w_carb*rd.uniform(0.8, 1.0), 4)
                            phi_dol = round(w_carb*rd.uniform(0.0, (1.0 - phi_cal)), 4)
                            phi_sd = round(w_carb - phi_cal - phi_dol, 4)
                            #
                            phi_qz = round(w_clast*rd.uniform(0.0, 1.0), 4)
                            phi_kfs = round(w_clast*rd.uniform(0.0, (1.0 - phi_qz)), 4)
                            phi_pl = round(w_clast - phi_qz - phi_kfs, 4)
                            #
                            phi_kln = round(w_clay*rd.uniform(0.0, 1.0), 4)
                            phi_mnt = round(w_clay*rd.uniform(0.0, (1.0 - phi_kln)), 4)
                            phi_ilt = round(1 - phi_cal - phi_dol - phi_sd - phi_qz - phi_kfs - phi_pl - phi_kln - phi_mnt, 4)
                            #
                            phi_total = phi_cal + phi_dol + phi_sd + phi_qz + phi_kfs + phi_pl + phi_kln + phi_mnt + phi_ilt
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.8 <= phi_cal <= 1.0 and 0.0 <= phi_dol <= 0.2 and 0.0 <= phi_sd <= 0.2 \
                                    and 0.0 <= phi_qz <= 0.2 and 0.0 <= phi_kfs <= 0.2 and 0.0 <= phi_pl <= 0.2 \
                                    and 0.0 <= phi_kln <= 0.2 and 0.0 <= phi_mnt <= 0.2 and 0.0 <= phi_ilt <= 0.2:
                                condition_2 = True
                        #
                    phi_minerals["Cal"] = abs(phi_cal)
                    phi_minerals["Dol"] = abs(phi_dol)
                    phi_minerals["Sd"] = abs(phi_sd)
                    phi_minerals["Qz"] = abs(phi_qz)
                    phi_minerals["Kfs"] = abs(phi_kfs)
                    phi_minerals["Pl"] = abs(phi_pl)
                    phi_minerals["Kln"] = abs(phi_kln)
                    phi_minerals["Mnt"] = abs(phi_mnt)
                    phi_minerals["Ilt"] = abs(phi_ilt)
                #
                rho_s = 0
                for key, value in phi_minerals.items():
                    rho_s += value*mineralogy[key]["rho"]
                    for element, value in mineralogy[key]["chemistry"].items():
                        if element not in elements_list:
                            elements_list.append(element)
                            w_elements[element] = 0.0
                #
                rho_s = round(rho_s, 3)
                #
                for key, value in phi_minerals.items():
                    w_result = round((phi_minerals[key] * mineralogy[key]["rho"]) / rho_s, 4)
                    w_minerals[key] = w_result
                #
                if porosity == None:
                    phi_helper = round(rd.uniform(0.0, 0.4), 4)
                else:
                    phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
                #
                rho = round((1 - phi_helper)*rho_s + phi_helper*self.data_water[2], 3)
                #
                old_index = elements_list.index("O")
                elements_list += [elements_list.pop(old_index)]
                #
                w_elements_total = 0.0
                for element in elements_list:
                    if element != "O":
                        for mineral, w_mineral in w_minerals.items():
                            if element in mineralogy[mineral]["chemistry"]:
                                value = round(w_mineral * mineralogy[mineral]["chemistry"][element], 4)
                                w_elements[element] += value
                                w_elements_total += value
                                #
                                w_elements[element] = round(w_elements[element], 4)
                    elif element == "O":
                        w_elements[element] += round(1 - w_elements_total, 4)
                        #
                        w_elements[element] = round(w_elements[element], 4)
                #
                if sum(w_minerals.values()) == 1.0 and sum(w_elements.values()) == 1.0:
                    #
                    condition = True
                #
                bulk_mod = 0.0
                shear_mod = 0.0
                gamma_ray = 0.0
                photoelectricity = 0.0
                for key, value in phi_minerals.items():
                    bulk_mod += phi_minerals[key] * mineralogy[key]["K"]
                    shear_mod += phi_minerals[key] * mineralogy[key]["G"]
                    gamma_ray += phi_minerals[key] * mineralogy[key]["GR"]
                    photoelectricity += phi_minerals[key] * mineralogy[key]["PE"]
                    #
                    bulk_mod = round(bulk_mod, 3)
                    shear_mod = round(shear_mod, 3)
                    gamma_ray = round(gamma_ray, 3)
                    photoelectricity = round(photoelectricity, 3)
                #
                w_list = []
                K_list = []
                G_list = []
                for key, mineral in mineralogy.items():
                    w_list.append(w_minerals[key])
                    K_list.append(mineral["K"])
                    G_list.append(mineral["G"])
                K_geo = elast.calc_geometric_mean(self, w_list, K_list)
                G_geo = elast.calc_geometric_mean(self, w_list, G_list)
                bulk_mod = K_geo
                shear_mod = G_geo
                #
                vP_s = round(((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(rho_s))**0.5, 3)
                vS_s = round(((shear_mod * 10 ** 9) / (rho_s)) ** 0.5, 3)
                #
                vP = (1 - phi_helper)*vP_s + phi_helper*self.data_water[4][0]
                vS = (1 - phi_helper)*vS_s
                vPvS = round(vP / vS, 3)
                #
                shear_mod = (vS**2 * rho)*10**(-9)
                bulk_mod = (vP**2 * rho - 4/3*shear_mod)*10**(-9)
                youngs_mod = round((9 * bulk_mod * shear_mod) / (3 * bulk_mod + shear_mod), 3)
                poisson_rat = round((3 * bulk_mod - 2 * shear_mod) / (6 * bulk_mod + 2 * shear_mod), 4)
            #
            for mineral, value in w_minerals.items():
                amounts_mineralogy[mineral].append(value)
            for element, value in w_elements.items():
                amounts_chemistry[element].append(value)
            #
            bulk_properties["rho_s"].append(rho_s)
            bulk_properties["rho"].append(rho)
            bulk_properties["phi"].append(phi_helper)
            bulk_properties["K"].append(bulk_mod)
            bulk_properties["G"].append(shear_mod)
            bulk_properties["E"].append(youngs_mod)
            bulk_properties["nu"].append(poisson_rat)
            bulk_properties["vP"].append(vP)
            bulk_properties["vS"].append(vS)
            bulk_properties["vPvS"].append(vPvS)
            bulk_properties["GR"].append(gamma_ray)
            bulk_properties["PE"].append(photoelectricity)
            #
            n += 1
        #
        ## EXPORT DATA
        #
        results = {}
        results["rock"] = "Limestone"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in w_minerals.items():
                single_amounts_mineralogy[mineral] = value
            for element, value in w_elements.items():
                single_amounts_chemistry[element] = value
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = phi_helper
            results["fluid"] = "water"
            results["rho_s"] = rho_s
            results["rho"] = rho
            results["vP"] = vP
            results["vS"] = vS
            results["vP/vS"] = vPvS
            results["K"] = bulk_mod
            results["G"] = shear_mod
            results["E"] = youngs_mod
            results["nu"] = poisson_rat
            results["GR"] = gamma_ray
            results["PE"] = photoelectricity
        #
        return results
    #
    def create_dolomite(self, number, porosity=None):
        #
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()  # variable
        #
        assemblage = [self.data_calcite, self.data_dolomite, self.data_magnesite, self.data_quartz, data_chlorite]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage:
            amounts_mineralogy[mineral["mineral"]] = []
            mineral_list.append(mineral["mineral"])
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        amounts_helper = []
        while n < number:
            w_total = 0
            n_minerals = 0
            for mineral in mineral_list:
                if mineral == "Cal":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Dol":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.8, 1.0), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.8 <= value <= 1.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Mgs":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.1), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.1:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Qz":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
            #
            if np.sum(amounts_helper) == 1.0 and n_minerals == len(mineral_list):
                for index, mineral in enumerate(mineral_list):
                    amounts_mineralogy[mineral].append(amounts_helper[index])
                n += 1
                amounts_helper.clear()
            else:
                n += 0
                amounts_helper.clear()
        #
        n = 0
        amounts_helper = {}
        while n < number:
            w_total = 0
            n_elements = 0
            rho_s_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.1, 0.4), 4)
            else:
                phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
            #
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            old_index = elements.index("O")
            elements += [elements.pop(old_index)]
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_dolomite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Dol"][n] * self.data_dolomite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_magnesite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Mgs"][n] * self.data_magnesite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Qz"][n] * self.data_quartz["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in data_chlorite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Chl"][n] * data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                #
                n_elements += 1
            #
            shear_factor = 1.0
            #
            if sum(amounts_helper.values()) == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Dol":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_dolomite["G"],
                            3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["PE"], 3)
                    elif mineral == "Mgs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_magnesite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["PE"], 3)
                    elif mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_chlorite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["PE"], 3)
                #
                rho_helper = round((1 - phi_helper) * rho_s_helper + phi_helper * self.data_water[2] / 1000, 3)
                youngsmod_helper = round(
                    (9 * bulkmod_helper * shearmod_helper) / (3 * bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round(
                    (3 * bulkmod_helper - 2 * shearmod_helper) / (6 * bulkmod_helper + 2 * shearmod_helper), 3)
                vP_helper = round(
                    ((bulkmod_helper * 10 ** 9 + 4 / 3 * shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vS_helper = round(((shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vPvS_helper_helper = round(vP_helper / vS_helper, 3)
                #
                bulk_properties["rho_s"].append(round(rho_s_helper, 3))
                bulk_properties["rho"].append(rho_helper)
                bulk_properties["K"].append(round(bulkmod_helper, 3))
                bulk_properties["G"].append(round(shearmod_helper, 3))
                bulk_properties["E"].append(youngsmod_helper)
                bulk_properties["nu"].append(poisson_helper)
                bulk_properties["vP"].append(vP_helper)
                bulk_properties["vS"].append(vS_helper)
                bulk_properties["vPvS"].append(vPvS_helper_helper)
                bulk_properties["GR"].append(round(gr_helper, 3))
                bulk_properties["PE"].append(round(pe_helper, 3))
                bulk_properties["phi"].append(round(phi_helper, 3))
                n += 1
        #
        results = {}
        results["rock"] = "Dolomite"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = bulk_properties["phi"][0]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"][0]
            results["rho"] = bulk_properties["rho"][0]
            results["vP"] = bulk_properties["vP"][0]
            results["vS"] = bulk_properties["vS"][0]
            results["vP/vS"] = bulk_properties["vPvS"][0]
            results["K"] = bulk_properties["K"][0]
            results["G"] = bulk_properties["G"][0]
            results["E"] = bulk_properties["E"][0]
            results["nu"] = bulk_properties["nu"][0]
            results["GR"] = bulk_properties["GR"][0]
            results["PE"] = bulk_properties["PE"][0]
        #
        return results
    #
class dolomite:
    #
    def __init__(self, fluid, actualThickness):
        self.fluid = fluid
        self.actualThickness = actualThickness
    #
    def createDolomite(self):
        # [symbol, atomic number, atomic mass, oxidation states, melting point, boiling point, density, electronegativity]
        chemH = ["H", 1, 1.0078, 1, 13.99, 20.271, 0.084, 2.2]
        chemC = ["C", 6, 12.009, 4, 0.0, 3915, 3510, 2.55]
        chemN = ["N", 7, 14.006, -3, 63.15, 77.355, 1.170, 3.04]
        chemO = ["O", 8, 15.999, -2, 54.3, 90.188, 1.33, 3.44]
        chemMg = ["Mg", 12, 24.304, 2, 923, 1363, 1740, 1.31]
        chemS = ["S", 16, 32.059, 6, 368.4, 717.8, 2060, 2.58]
        chemAr = ["Ar", 18, 39.948, 0.0, 83.81, 87.302, 1.66, 0.0]
        chemCa = ["Ca", 20, 40.078, 2, 1115, 1757, 1540, 1.0]
        chemFe = ["Fe", 26, 55.845, 3, 1811, 3134, 7870, 1.83]
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemDolomite = minerals.carbonates.dolomite("")
        chemCalcite = minerals.carbonates.calcite("")
        chemQuartz = minerals.oxides.quartz("")
        chemAlkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chemPlagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        chemMagnesite = minerals.carbonates.magnesite("")
        chemKaolinite = minerals.phyllosilicates.kaolinite("")
        chemChlorite = minerals.phyllosilicates.chamosite("")
        chemIllite = minerals.phyllosilicates.illite("")
        #
        # [molar mass, density, bulk modulus, vP]
        chemWater = [2 * chemH[2] + 1 * chemO[2], 997, 2.08, 1444]
        chemOil = [0.83 * chemC[2] + 0.11 * chemH[2] + 0.05 * chemS[2] + 0.005 * (chemO[2] + chemN[2]), 0.9, 1.35, 1225]
        chemGas = [0.83 * chemC[2] + 0.11 * chemH[2] + 0.05 * chemS[2] + 0.005 * (chemO[2] + chemN[2]), 0.8, 0.081, 475]
        #
        dolomite = []
        #
        cond = False
        composition = []
        while cond == False:
            magicnumber = rd.randint(0, 2)
            if magicnumber == 0:    # dolomite
                xMgCO3 = rd.randint(50, 100)/100            # Mg-bearing carbonates
                xDolomite2 = rd.randint(75, 100)/100
                xMagnesite2 = 1 - xDolomite2
                xDolomite = xMgCO3*xDolomite2
                xMagnesite = xMgCO3*xMagnesite2
                xCalcite = rd.randint(0, 25)/100            # Ca-bearing carbonates
                xImpurities = rd.randint(0, 25)/100         # Impurities
                magicnumber2 = rd.randint(0, 2)
                if magicnumber2 == 0:                       # Qz + Fsp
                    xQuartz2 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz2
                    xAlkalifeldspar3 = rd.randint(0, 100)/100
                    xPlagioclase3 = 1 - xAlkalifeldspar3
                    xAlkalifeldspar2 = xFeldspars*xAlkalifeldspar3
                    xPlagioclase2 = xFeldspars*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xKaolinite = 0.0
                    xChlorite = 0.0
                    xIllite = 0.0
                elif magicnumber2 == 1:                     # Clays
                    xKaolinite2 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite2,2)*100)
                    xChlorite2 = round(rd.randint(0, limChl)/100,2)
                    xIllite2 = round(1 - xKaolinite2 - xChlorite2,2)
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
                    xQuartz = 0.0
                    xAlkalifeldspar = 0.0
                    xPlagioclase = 0.0
                elif magicnumber2 == 2:                     # Qz + Fsp + Clays
                    xQzFsp = rd.randint(0, 100)/100
                    xQuartz3 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz3
                    xAlkalifeldspar4 = rd.randint(0, 100)/100
                    xPlagioclase4 = 1 - xAlkalifeldspar4
                    xAlkalifeldspar3 = xFeldspars*xAlkalifeldspar4
                    xPlagioclase3 = xFeldspars*xPlagioclase4
                    xQuartz2 = xQzFsp*xQuartz3
                    xAlkalifeldspar2 = xQzFsp*xAlkalifeldspar3
                    xPlagioclase2 = xQzFsp*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xClays = 1 - xQzFsp
                    xKaolinite3 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite3,2)*100)
                    xChlorite3 = round(rd.randint(0, limChl)/100,2)
                    xIllite3 = round(1 - xKaolinite3 - xChlorite3,2)
                    xKaolinite2 = xClays*xKaolinite3
                    xChlorite2 = xClays*xChlorite3
                    xIllite2 = xClays*xIllite3
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
            elif magicnumber == 1:    # clean dolomite
                xMgCO3 = rd.randint(80, 100)/100            # Mg-bearing carbonates
                xDolomite2 = rd.randint(90, 100)/100
                xMagnesite2 = 1 - xDolomite2
                xDolomite = xMgCO3*xDolomite2
                xMagnesite = xMgCO3*xMagnesite2
                xCalcite = rd.randint(0, 10)/100            # Ca-bearing carbonates
                xImpurities = rd.randint(0, 10)/100         # Impurities
                magicnumber2 = rd.randint(0, 2)
                if magicnumber2 == 0:                       # Qz + Fsp
                    xQuartz2 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz2
                    xAlkalifeldspar3 = rd.randint(0, 100)/100
                    xPlagioclase3 = 1 - xAlkalifeldspar3
                    xAlkalifeldspar2 = xFeldspars*xAlkalifeldspar3
                    xPlagioclase2 = xFeldspars*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xKaolinite = 0.0
                    xChlorite = 0.0
                    xIllite = 0.0
                elif magicnumber2 == 1:                     # Clays
                    xKaolinite2 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite2,2)*100)
                    xChlorite2 = round(rd.randint(0, limChl)/100,2)
                    xIllite2 = round(1 - xKaolinite2 - xChlorite2,2)
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
                    xQuartz = 0.0
                    xAlkalifeldspar = 0.0
                    xPlagioclase = 0.0
                elif magicnumber2 == 2:                     # Qz + Fsp + Clays
                    xQzFsp = rd.randint(0, 100)/100
                    xQuartz3 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz3
                    xAlkalifeldspar4 = rd.randint(0, 100)/100
                    xPlagioclase4 = 1 - xAlkalifeldspar4
                    xAlkalifeldspar3 = xFeldspars*xAlkalifeldspar4
                    xPlagioclase3 = xFeldspars*xPlagioclase4
                    xQuartz2 = xQzFsp*xQuartz3
                    xAlkalifeldspar2 = xQzFsp*xAlkalifeldspar3
                    xPlagioclase2 = xQzFsp*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xClays = 1 - xQzFsp
                    xKaolinite3 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite3,2)*100)
                    xChlorite3 = round(rd.randint(0, limChl)/100,2)
                    xIllite3 = round(1 - xKaolinite3 - xChlorite3,2)
                    xKaolinite2 = xClays*xKaolinite3
                    xChlorite2 = xClays*xChlorite3
                    xIllite2 = xClays*xIllite3
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
            elif magicnumber == 2:    # dirty dolomite
                xMgCO3 = rd.randint(50, 80)/100            # Mg-bearing carbonates
                xDolomite2 = rd.randint(75, 100)/100
                xMagnesite2 = 1 - xDolomite2
                xDolomite = xMgCO3*xDolomite2
                xMagnesite = xMgCO3*xMagnesite2
                xCalcite = rd.randint(0, 10)/100            # Ca-bearing carbonates
                xImpurities = rd.randint(20, 40)/100         # Impurities
                magicnumber2 = rd.randint(0, 2)
                if magicnumber2 == 0:                       # Qz + Fsp
                    xQuartz2 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz2
                    xAlkalifeldspar3 = rd.randint(0, 100)/100
                    xPlagioclase3 = 1 - xAlkalifeldspar3
                    xAlkalifeldspar2 = xFeldspars*xAlkalifeldspar3
                    xPlagioclase2 = xFeldspars*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xKaolinite = 0.0
                    xChlorite = 0.0
                    xIllite = 0.0
                elif magicnumber2 == 1:                     # Clays
                    xKaolinite2 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite2,2)*100)
                    xChlorite2 = round(rd.randint(0, limChl)/100,2)
                    xIllite2 = round(1 - xKaolinite2 - xChlorite2,2)
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
                    xQuartz = 0.0
                    xAlkalifeldspar = 0.0
                    xPlagioclase = 0.0
                elif magicnumber2 == 2:                     # Qz + Fsp + Clays
                    xQzFsp = rd.randint(0, 100)/100
                    xQuartz3 = rd.randint(0, 100)/100
                    xFeldspars = 1 - xQuartz3
                    xAlkalifeldspar4 = rd.randint(0, 100)/100
                    xPlagioclase4 = 1 - xAlkalifeldspar4
                    xAlkalifeldspar3 = xFeldspars*xAlkalifeldspar4
                    xPlagioclase3 = xFeldspars*xPlagioclase4
                    xQuartz2 = xQzFsp*xQuartz3
                    xAlkalifeldspar2 = xQzFsp*xAlkalifeldspar3
                    xPlagioclase2 = xQzFsp*xPlagioclase3
                    xQuartz = xImpurities*xQuartz2
                    xAlkalifeldspar = xImpurities*xAlkalifeldspar2
                    xPlagioclase = xImpurities*xPlagioclase2
                    xClays = 1 - xQzFsp
                    xKaolinite3 = round(rd.randint(0, 100)/100,2)
                    limChl = int(round(1 - xKaolinite3,2)*100)
                    xChlorite3 = round(rd.randint(0, limChl)/100,2)
                    xIllite3 = round(1 - xKaolinite3 - xChlorite3,2)
                    xKaolinite2 = xClays*xKaolinite3
                    xChlorite2 = xClays*xChlorite3
                    xIllite2 = xClays*xIllite3
                    xKaolinite = xImpurities*xKaolinite2
                    xChlorite = xImpurities*xChlorite2
                    xIllite = xImpurities*xIllite2
            sumMin = round(xCalcite,2) + round(xDolomite,2) + round(xMagnesite,2) + round(xQuartz,2) + round(xAlkalifeldspar,2) + round(xPlagioclase,2) + round(xKaolinite,2) + round(xChlorite,2) + round(xIllite,2)
            if sumMin == 1:
                cond = True
                composition.extend([["Cal", round(xCalcite,2), round(chemCalcite[1],2)], ["Dol", round(xDolomite,2), round(chemDolomite[1],2)], ["Mgs", round(xMagnesite,2), round(chemMagnesite[1],2)], ["Qz", round(xQuartz,2), round(chemQuartz[1],2)], ["Kfs", round(xAlkalifeldspar, 2), round(chemAlkalifeldspar[1][0],2), round(chemAlkalifeldspar[1][1],2)], ["Pl", round(xPlagioclase, 2), round(chemPlagioclase[1][0],2), round(chemPlagioclase[1][1],2)], ["Kln", round(xKaolinite,2), round(chemKaolinite[1],2)], ["Chl", round(xChlorite,2), round(chemChlorite[1],2)], ["Ilt", round(xIllite,2), round(chemIllite[1],2)]])
            else:
                cond = False
        xCalcite = composition[0][1]
        xDolomite = composition[1][1]
        xMagnesite = composition[2][1]
        xQuartz = composition[3][1]
        xAlkalifeldspar = composition[4][1]
        xPlagioclase = composition[5][1]
        xKaolinite = composition[6][1]
        xChlorite = composition[7][1]
        xIllite = composition[8][1]
        dolomite.append(composition)
        mineralogy = [chemCalcite, chemDolomite, chemMagnesite, chemQuartz, chemAlkalifeldspar, chemPlagioclase, chemKaolinite, chemChlorite, chemIllite]
        #
        rhoSolid = 1.25*(xCalcite*chemCalcite[2] + xDolomite*chemDolomite[2] + xMagnesite*chemMagnesite[2] + xQuartz*chemQuartz[2] + xAlkalifeldspar*chemAlkalifeldspar[2] + xPlagioclase*chemPlagioclase[2] + xKaolinite*chemKaolinite[2] + xChlorite*chemChlorite[2] + xIllite*chemIllite[2]) / 1000
        X = [xCalcite, xDolomite, xMagnesite, xQuartz, xAlkalifeldspar, xPlagioclase, xKaolinite, xChlorite, xIllite]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_arith = elast.calc_arithmetic_mean(self, X, K_list)
        G_arith = elast.calc_arithmetic_mean(self, X, G_list)
        K_solid = K_arith
        G_solid = G_arith
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if self.actualThickness <= 1000:
            phi = randint(25, 30)/100
        elif self.actualThickness > 1000 and self.actualThickness <= 2000:
            phi = randint(20, 25)/100
        elif self.actualThickness > 2000 and self.actualThickness <= 3000:
            phi = randint(15, 20)/100
        elif self.actualThickness > 3000 and self.actualThickness <= 4000:
            phi = randint(10, 15)/100
        elif self.actualThickness > 4000:
            phi = randint(5, 10)/100
        #
        if self.fluid == "water":
            rho = (1 - phi) * rhoSolid + phi * chemWater[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemWater[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemWater[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xMagnesite*chemMagnesite[5][0] + xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xKaolinite*chemKaolinite[5][0] + xChlorite*chemChlorite[5][0] + xIllite*chemIllite[5][0]
            PE = xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xMagnesite*chemMagnesite[5][1] + xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xKaolinite*chemKaolinite[5][1] + xChlorite*chemChlorite[5][1] + xIllite*chemIllite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xCalcite*chemCalcite[3][3] + xDolomite*chemDolomite[3][3] + xMagnesite*chemMagnesite[3][3] + xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xKaolinite*chemKaolinite[3][3] + xChlorite*chemChlorite[3][3] + xIllite*chemIllite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            dolomite.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[1] / 1000, 6)])
            dolomite.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            dolomite.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[3], 2)])
            dolomite.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            dolomite.append("water")
            dolomite.append([GR, PE])
        elif self.fluid == "oil":
            rho = (1 - phi) * rhoSolid + phi * chemOil[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemOil[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemOil[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xMagnesite*chemMagnesite[5][0] + xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xKaolinite*chemKaolinite[5][0] + xChlorite*chemChlorite[5][0] + xIllite*chemIllite[5][0]
            PE = xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xMagnesite*chemMagnesite[5][1] + xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xKaolinite*chemKaolinite[5][1] + xChlorite*chemChlorite[5][1] + xIllite*chemIllite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xCalcite*chemCalcite[3][3] + xDolomite*chemDolomite[3][3] + xMagnesite*chemMagnesite[3][3] + xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xKaolinite*chemKaolinite[3][3] + xChlorite*chemChlorite[3][3] + xIllite*chemIllite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            dolomite.append([round(rho, 3), round(rhoSolid, 3), round(chemOil[1] / 1000, 6)])
            dolomite.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            dolomite.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemOil[3], 2)])
            dolomite.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            dolomite.append("oil")
            dolomite.append([GR, PE])
        elif self.fluid == "gas":
            rho = (1 - phi) * rhoSolid + phi * chemGas[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemGas[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemGas[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xMagnesite*chemMagnesite[5][0] + xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xKaolinite*chemKaolinite[5][0] + xChlorite*chemChlorite[5][0] + xIllite*chemIllite[5][0]
            PE = xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xMagnesite*chemMagnesite[5][1] + xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xKaolinite*chemKaolinite[5][1] + xChlorite*chemChlorite[5][1] + xIllite*chemIllite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xCalcite*chemCalcite[3][3] + xDolomite*chemDolomite[3][3] + xMagnesite*chemMagnesite[3][3] + xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xKaolinite*chemKaolinite[3][3] + xChlorite*chemChlorite[3][3] + xIllite*chemIllite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            dolomite.append([round(rho, 3), round(rhoSolid, 3), round(chemGas[1] / 1000, 6)])
            dolomite.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            dolomite.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemGas[3], 2)])
            dolomite.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            dolomite.append("gas")
            dolomite.append([GR, PE])
        #
        #  dolomite = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return dolomite
    #
    def create_simple_dolomite(self, w_Mg=None, w_Ca=None, w_Fe=None, amounts=None, porosity=None, dict=None):
        #
        results = {}
        results["rock"] = "dolomite rock"
        #
        self.w_Mg = w_Mg
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # mineralogy
        chem_dol = minerals.carbonates.dolomite("")
        chem_ank = minerals.carbonates.ankerite(self, "Pure")
        chem_sd = minerals.carbonates.siderite("")
        chem_cal = minerals.carbonates.calcite("")
        chem_qz = minerals.oxides.quartz("")
        chem_kfs = minerals.feldspars.alkalifeldspar(self, "Kfs")
        chem_pl = minerals.feldspars.plagioclase(self, "Pl")
        chem_kln = minerals.phyllosilicates.kaolinite("")
        chem_py = minerals.sulfides.pyrite("")
        #
        mineralogy = [chem_dol, chem_ank, chem_sd, chem_cal, chem_qz, chem_kfs, chem_pl, chem_kln, chem_py]
        #for i in range(len(mineralogy)):
        #    print(mineralogy[i][0], mineralogy[i][6])
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Mg == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                magicnumber = rd.randint(0, 6)
                if magicnumber == 0:    # Dol-rich
                    w_carb = round(rd.randint(90, 100)/100, 4)
                    w_dol2 = rd.randint(90, 100)/100
                    w_ank2 = rd.randint(0, int((1-w_dol2)*100))/100
                    w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_kln2 = w_clay
                    w_kln = round(w_clay*w_kln2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 1:    # Ank-rich
                    w_carb = round(rd.randint(90, 100)/100, 4)
                    w_dol2 = rd.randint(60, 80)/100
                    w_ank2 = rd.randint(10, int((1-w_dol2)*100))/100
                    w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_kln2 = w_clay
                    w_kln = round(w_clay*w_kln2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 2:    # Sd-rich
                    w_carb = round(rd.randint(90, 100)/100, 4)
                    w_dol2 = rd.uniform(0.6, 0.8)
                    w_sd2 = rd.uniform(0.1, (1-w_dol2))
                    w_ank2 = rd.randint(0, int((1-w_dol2-w_sd2)*100))/100
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_kln2 = w_clay
                    w_kln = round(w_clay*w_kln2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 3:    # Mixed
                    w_carb = round(rd.randint(75, 85)/100, 4)
                    w_dol2 = rd.randint(80, 90)/100
                    w_ank2 = rd.randint(0, int((1-w_dol2)*100))/100
                    w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_kln2 = w_clay
                    w_kln = round(w_clay*w_kln2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 4:    # Qz-rich
                    w_carb = round(rd.randint(75, 85)/100, 4)
                    w_dol2 = rd.randint(80, 90)/100
                    w_ank2 = rd.randint(0, int((1-w_dol2)*100))/100
                    w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_qz = round(rd.randint(10, int((1-w_carb)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_kln2 = w_clay
                    w_kln = round(w_clay*w_kln2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 5:    # Fsp-rich
                    w_carb = round(rd.randint(75, 85)/100, 4)
                    w_dol2 = rd.randint(80, 90)/100
                    w_ank2 = rd.randint(0, int((1-w_dol2)*100))/100
                    w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_fsp = round(rd.randint(10, int((1-w_carb)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_qz = round(rd.randint(0, int((1-w_carb-w_fsp)*100))/100, 4)
                    w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                    w_kln2 = w_clay
                    w_kln = round(w_clay*w_kln2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
                elif magicnumber == 6:    # Clay-rich
                    w_carb = round(rd.uniform(0.75, 0.85), 4)
                    w_dol2 = rd.uniform(0.8, 0.9)
                    w_ank2 = rd.uniform(0, (1-w_dol2))
                    w_sd2 = rd.uniform(0, (1-w_dol2-w_ank2))
                    w_cal2 = 1-w_dol2-w_ank2-w_sd2
                    w_dol = round(w_carb*w_dol2, 4)
                    w_ank = round(w_carb*w_ank2, 4)
                    w_sd = round(w_carb*w_sd2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_clay = round(rd.uniform(0.1, (1-w_carb)), 4)
                    w_kln = w_clay
                    w_qz = round(rd.randint(0, int((1-w_carb-w_clay)*100))/100, 4)
                    w_fsp = round(rd.randint(0, int((1-w_carb-w_qz-w_clay)*100))/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Mg != None:
                condition = False
                while condition == False:
                    w_dol = round(self.w_Mg/chem_dol[6][2], 4)
                    if w_dol <= 1.0:
                        w_carb = round(rd.randint(75, 100)/100, 4)
                    else:
                        print("Amount of Dol is larger than 1.")
                        break
                    w_dol2 = w_dol/w_carb
                    if w_dol2 <= 1.0:
                        condition = True
                    else:
                        continue
                w_ank2 = rd.randint(0, int((1-w_dol2)*100))/100
                w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                w_cal2 = 1-w_dol2-w_ank2-w_sd2
                w_ank = round(w_carb*w_ank2, 4)
                w_sd = round(w_carb*w_sd2, 4)
                w_cal = round(w_carb*w_cal2, 4)
                w_qz = round(rd.randint(0, int((1-w_carb)*100))/100, 4)
                w_fsp = round(rd.randint(0, int((1-w_carb-w_qz)*100))/100, 4)
                w_kfs2 = rd.randint(0, 100)/100
                w_pl2 = 1-w_kfs2
                w_kfs = round(w_fsp*w_kfs2, 4)
                w_pl = round(w_fsp*w_pl2, 4)
                w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                w_kln2 = w_clay
                w_kln = round(w_clay*w_kln2, 4)
                w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Ca != None:
                condition = False
                while condition == False:
                    w_fsp = round(rd.randint(0, 10)/100, 4)
                    w_kfs2 = rd.randint(0, 100)/100
                    w_pl2 = 1-w_kfs2
                    w_kfs = round(w_fsp*w_kfs2, 4)
                    w_pl = round(w_fsp*w_pl2, 4)
                    w_carb = round(rd.randint(75, 90)/100, 4)
                    w_ank2 = rd.randint(0, 15)/100
                    w_cal2 = rd.randint(0, 15)/100
                    w_ank = round(w_carb*w_ank2, 4)
                    w_cal = round(w_carb*w_cal2, 4)
                    w_dol = round((self.w_Ca - w_ank*chem_ank[6][3] - w_cal*chem_cal[6][2] - w_pl*chem_pl[6][4])/chem_dol[6][3], 4)
                    w_dol2 = w_dol/w_carb
                    w_sd2 = 1-w_dol2-w_ank2-w_cal2
                    w_sd = round(w_carb*w_sd2, 4)
                    if w_dol <= 1.0 and w_fsp+w_carb <= 1.0:
                        condition = True
                    else:
                        print("Amount of Fsp+Carb is larger than 1.")
                        continue
                w_qz = round(rd.randint(0, int((1-w_carb-w_fsp)*100))/100, 4)
                w_clay = round(rd.randint(0, int((1-w_carb-w_qz-w_fsp)*100))/100, 4)
                w_kln2 = w_clay
                w_kln = round(w_clay*w_kln2, 4)
                w_py = round(1-w_carb-w_qz-w_fsp-w_clay, 4)
            elif self.w_Fe != None:
                w_carb = round(rd.randint(75, 85)/100, 4)
                w_dol2 = rd.randint(80, 90)/100
                w_ank2 = rd.randint(0, int((1-w_dol2)*100))/100
                w_sd2 = rd.randint(0, int((1-w_dol2-w_ank2)*100))/100
                w_cal2 = 1-w_dol2-w_ank2-w_sd2
                w_dol = round(w_carb*w_dol2, 4)
                w_ank = round(w_carb*w_ank2, 4)
                w_sd = round(w_carb*w_sd2, 4)
                w_cal = round(w_carb*w_cal2, 4)
                w_py = round((self.w_Fe - w_ank*chem_ank[6][5] - w_sd*chem_sd[6][2])/chem_py[6][1], 4)
                w_qz = round(rd.randint(0, int((1-w_carb-w_py)*100))/100, 4)
                w_fsp = round(rd.randint(0, int((1-w_carb-w_qz-w_py)*100))/100, 4)
                w_kfs2 = rd.randint(0, 100)/100
                w_pl2 = 1-w_kfs2
                w_kfs = round(w_fsp*w_kfs2, 4)
                w_pl = round(w_fsp*w_pl2, 4)
                w_clay = round(1-w_py-w_carb-w_qz-w_fsp, 4)
                w_kln2 = w_clay
                w_kln = round(w_clay*w_kln2, 4)
            elif type(self.amounts) is list:
                w_dol = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_ank = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_sd = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_qz = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[7], 0.025)), 4)
                w_py = round(1-w_dol-w_ank-w_sd-w_cal-w_qz-w_kfs-w_pl-w_kln, 4)
            #
            if w_dol >= 0.0 and w_ank >= 0.0 and w_sd >= 0.0 and w_cal >= 0.0 and w_qz >= 0.0 and w_kfs  >= 0.0 and w_pl >= 0.0 and w_kln >= 0.0 and w_py >= 0.0:
                sumMin = round(w_dol + w_ank + w_sd + w_cal + w_qz + w_kfs + w_pl + w_kln + w_py, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_kln*chem_kln[6][0], 4)
            w_C = round(w_dol*chem_dol[6][0] + w_ank*chem_ank[6][0] + w_sd*chem_sd[6][0] + w_cal*chem_cal[6][0], 4)
            #w_O = round(w_dol*chem_dol[6][1] + w_ank*chem_ank[6][1] + w_sd*chem_sd[6][1] + w_cal*chem_cal[6][1] + w_qz*chem_qz[6][0] + w_kfs*chem_kfs[6][0] + w_pl*chem_pl[6][0] + w_kln*chem_kln[6][1], 4)
            w_Na = round(w_kfs*chem_kfs[6][1] + w_pl*chem_pl[6][1], 4)
            w_Mg = round(w_dol*chem_dol[6][2], 4)
            w_Al = round(w_kfs*chem_kfs[6][2] + w_pl*chem_pl[6][2] + w_kln*chem_kln[6][2], 4)
            w_Si = round(w_qz*chem_qz[6][1] + w_kfs*chem_kfs[6][3] + w_pl*chem_pl[6][3] + w_kln*chem_kln[6][3], 4)
            w_S = round(w_py*chem_py[6][0], 4)
            w_K = round(w_kfs*chem_kfs[6][4], 4)
            w_Ca = round(w_dol*chem_dol[6][3] + w_ank*chem_ank[6][3] + w_cal*chem_cal[6][2] + w_pl*chem_pl[6][4], 4)
            w_Fe = round(w_ank*chem_ank[6][5] + w_sd*chem_sd[6][2] + w_py*chem_py[6][1], 4)
            w_O = round(1 - w_H - w_C - w_Na - w_Mg - w_Al - w_Si - w_S - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_C + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                #composition.extend((["Dol", w_dol, round(chem_dol[1], 2)], ["Ank", w_ank, round(chem_ank[1][0], 2)], ["Sd", w_sd, round(chem_sd[1], 2)], ["Cal", w_cal, round(chem_cal[1], 2)], ["Qz", w_qz, round(chem_qz[1], 2)], ["Kfs", w_kfs, round(chem_kfs[1][0], 2), round(chem_kfs[1][1], 2)], ["Pl", w_pl, round(chem_pl[1][0], 2), round(chem_pl[1][1], 2)], ["Kln", w_kln, round(chem_kln[1], 2)], ["Py", w_py, round(chem_py[1], 2)]))
                composition.extend((["Dol", "Ank", "Sd", "Cal", "Qz", "Kfs", "Pl", "Kln", "Py"]))
                concentrations = [w_H, w_C, w_O, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe]
                amounts = [w_dol, w_ank, w_sd, w_cal, w_qz, w_kfs, w_pl, w_kln, w_py]
            else:
                cond = False
        #
        element_list = ["H", "C", "O", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Fe"]
        mineral_list = ["Dol", "Ank", "Sd", "Cal", "Qz", "Kfs", "Pl", "Kln", "Py"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_dol*chem_dol[2] + w_ank*chem_ank[2] + w_sd*chem_sd[2] + w_cal*chem_cal[2] + w_qz*chem_qz[2] + w_kfs*chem_kfs[2] + w_pl*chem_pl[2] + w_kln*chem_kln[2] + w_py*chem_py[2]) / 1000
        X = [w_dol, w_ank, w_sd, w_cal, w_qz, w_kfs, w_pl, w_kln, w_py]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.randint(25, 30)/100
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.randint(20, 25)/100
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.randint(15, 20)/100
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.randint(10, 15)/100
            elif self.actualThickness > 4000:
                phi = rd.randint(5, 10)/100
        else:
            phi = porosity
        #
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_dol*chem_dol[5][0] + w_ank*chem_ank[5][0] + w_sd*chem_sd[5][0] + w_cal*chem_cal[5][0] + w_qz*chem_qz[5][0] + w_kfs*chem_kfs[5][0] + w_pl*chem_pl[5][0] + w_kln*chem_kln[5][0] + w_py*chem_py[5][0]
        PE = w_dol*chem_dol[5][1] + w_ank*chem_ank[5][1] + w_sd*chem_sd[5][1] + w_cal*chem_cal[5][1] + w_qz*chem_qz[5][1] + w_kfs*chem_kfs[5][1] + w_pl*chem_pl[5][1] + w_kln*chem_kln[5][1] + w_py*chem_py[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_dol*chem_dol[3][3] + w_ank*chem_ank[3][3] + w_sd*chem_sd[3][3] + w_cal*chem_cal[3][3] + w_qz*chem_qz[3][3] + w_kfs*chem_kfs[3][3] + w_pl*chem_pl[3][3] + w_kln*chem_kln[3][3] + w_py*chem_py[3][3]
        #
        results["rho"] = round(rho*1000, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vP/vS, 4)
        results["G"] = round(G_bulk*10**(-6), 4)
        results["K"] = round(K_bulk*10**(-6), 4)
        results["E"] = round(E_bulk*10**(-6), 4)
        results["nu"] = round(poisson_mineralogical, 4)
        results["GR"] = round(GR, 4)
        results["PE"] = round(PE, 4)
        #
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append("water")
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        #
        if dict == False:
            return data
        else:
            return results
#
# Carbonates
class Carbonates():
    """ Class that generates geophysical and geochemical data of carbonate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1):
        if self.mineral in ["Cal", "Calcite"]:
            if number > 1:
                data = [self.create_calcite() for n in range(number)]
            else:
                data = self.create_calcite()
        elif self.mineral in ["Dol", "Dolomite"]:
            if number > 1:
                data = [self.create_dolomite() for n in range(number)]
            else:
                data = self.create_dolomite()
        elif self.mineral in ["Mgs", "Magnesite"]:
            if number > 1:
                data = [self.create_magnesite() for n in range(number)]
            else:
                data = self.create_magnesite()
        elif self.mineral in ["Sd", "Siderite"]:
            if number > 1:
                data = [self.create_siderite() for n in range(number)]
            else:
                data = self.create_siderite()
        elif self.mineral in ["Org", "Organic Matter"]:
            if number > 1:
                data = [self.create_organic_matter() for n in range(number)]
            else:
                data = self.create_organic_matter()
        elif self.mineral in ["Rdc", "Rhodochrosite"]:
            if number > 1:
                data = [self.create_rhodochrosite() for n in range(number)]
            else:
                data = self.create_rhodochrosite()
        elif self.mineral in ["Arg", "Aragonite"]:
            if number > 1:
                data = [self.create_aragonite() for n in range(number)]
            else:
                data = self.create_aragonite()
        elif self.mineral in ["Cer", "Cerussite"]:
            if number > 1:
                data = [self.create_cerussite() for n in range(number)]
            else:
                data = self.create_cerussite()
        elif self.mineral in ["Ank", "Ankerite"]:
            if number > 1:
                data = [self.create_ankerite() for n in range(number)]
            else:
                data = self.create_ankerite()
        elif self.mineral in ["Az", "Azurite"]:
            if number > 1:
                data = [self.create_azurite() for n in range(number)]
            else:
                data = self.create_azurite()
        elif self.mineral in ["Mal", "Malachite"]:
            if number > 1:
                data = [self.create_malachite() for n in range(number)]
            else:
                data = self.create_malachite()
        #
        return data
    #
    def create_calcite(self):   # CaCO3
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Ca"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
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
                minors = ["Mn", "Fe", "Zn", "Co", "Ba", "Sr", "Pb", "Mg", "Cu", "Al", "Ni", "V", "Cr", "Mo"]
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
        mineral = "Cal"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + carbon[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.99, 17.06], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
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
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 2*10**12
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
    def create_dolomite(self):   # CaMg(CO3)2
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Mg", "Ca"]
        majors_data = np.array([["C", carbon[1], 2, carbon[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Mg", magnesium[1], 1, magnesium[2]], ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
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
                minors = ["Fe", "Mn", "Co", "Pb", "Zn"]
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
        mineral = "Dol"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + magnesium[2] + 2*(carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.81, 16.01], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 3
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
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
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 5*10**3
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
    def create_magnesite(self):   # Mg(CO3)
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        majors_name = ["C", "O", "Mg"]
        majors_data = np.array([["C", carbon[1], 2, carbon[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Mg", magnesium[1], 1, magnesium[2]]], dtype=object)
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
                minors = ["Fe", "Mn", "Ca", "Co", "Ni"]
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
        mineral = "Mgs"
        #
        # Molar mass
        molar_mass_pure = magnesium[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.63, 15.03], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
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
    def create_siderite(self):   # Fe(CO3)
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["C", "O", "Fe"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Fe", iron[1], 1, iron[2]]], dtype=object)
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
                minors = ["Mn", "Mg", "Ca", "Zn", "Co"]
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
        mineral = "Sd"
        #
        # Molar mass
        molar_mass_pure = iron[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.69, 15.38], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
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
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 70
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
    def create_organic_matter(self):
        # CHEMISTRY
        carbohydrates = Organics.carbohydrates("")
        lignin = Organics.lignin("")
        lipid = Organics.lipid("")
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
        molar_mass_pure = w_ch*(0.06*hydrogen[2] + 0.44*carbon[2] + 0.50*oxygen[2]) \
                          + w_lg*(0.06*hydrogen[2] + 0.63*carbon[2] + 0.003*nitrogen[2] + 0.31*oxygen[2] + 0.001*sulfur[2]) \
                          + w_lp*(0.10*hydrogen[2] + 0.80*carbon[2] + 0.10*oxygen[2])
        #
        majors_data = np.array([["H", hydrogen[1], w_ch*0.06*hydrogen[2] + w_lg*0.06*hydrogen[2] + w_lp*0.10*hydrogen[2], hydrogen[2]],
                                ["C", carbon[1], w_ch*0.44*carbon[2] + w_lg*0.63*carbon[2] + w_lp*0.80*carbon[2], carbon[2]],
                                ["N", nitrogen[1], w_lp*0.003*nitrogen[2], nitrogen[2]],
                                ["O", oxygen[1], w_ch*0.50*oxygen[2] + w_lg*0.31*oxygen[2] + w_lp*0.10*oxygen[2], oxygen[2]],
                                ["S", sulfur[1], w_lp*0.001*sulfur[2], sulfur[2]]], dtype=object)
        #
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
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
            data.append([round(GR, 2), round(PE, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_rhodochrosite(self):   # Mn(CO3)
        #
        name = "Rdc"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["C", "O", "Mn"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
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
                minors = ["Fe", "Ca", "Mg", "Zn", "Co", "Cd"]
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
        molar_mass_pure = manganese[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.777, 15.67], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 99*10**9
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
    def create_aragonite(self):   # Ca(CO3)
        #
        name = "Arg"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Ca"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
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
                minors = ["Sr", "Pb", "Zn"]
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
        molar_mass_pure = calcium[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.959, 7.968, 5.741], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
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
    def create_cerussite(self):   # Pb(CO3)
        #
        name = "Cer"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["C", "O", "Pb"]
        majors_data = np.array([["C", carbon[1], 1, carbon[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Pb", lead[1], 1, lead[2]]], dtype=object)
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
        molar_mass_pure = lead[2] + (carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.195, 8.436, 6.152], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 56*10**9
        # Shear modulus
        G = 21*10**9
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
    def create_ankerite(self):   # Ca(Fe,Mg)(CO3)2
        #
        name = "Ank"
        #
        # Major elements
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["C", "O", "Mg", "Fe"]
        #
        x = round(rd.uniform(0, 1.0), 4)
        #
        majors_data = np.array([["C", carbon[1], 2, carbon[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Mg", magnesium[1], (1-x), magnesium[2]], ["Ca", calcium[1], 1, calcium[2]],
                                ["Fe", iron[1], x, iron[2]]], dtype=object)
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
                minors = ["Mn"]
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
        molar_mass_pure = calcium[2] + (x*iron[2] + (1-x)*magnesium[2]) + 2*(carbon[2] + 3*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[4.83, 16.167], [], "trigonal"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 3
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        dataV_Mg = CrystalPhysics([[4.83, 16.167], [], "trigonal"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 3
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Fe.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        V = x*V_Fe + (1-x)*V_Mg
        V_m = x*V_m_Fe + (1-x)*V_m_Mg
        rho = x*rho_Fe + (1-x)*rho_Mg
        rho_e = x*rho_e_Fe + (1-x)*rho_e_Mg
        # Bulk modulus
        K_Fe = 73*10**9
        K_Mg = 89*10**9
        K = x*K_Fe + (1-x)*K_Mg
        # Shear modulus
        G_Fe = 32*10**9
        G_Mg = 44*10**9
        G = x*G_Fe + (1-x)*G_Mg
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
    def create_azurite(self):   # Cu3(CO3)2(OH)2
        #
        name = "Az"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["H", "C", "O", "Cu"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["C", carbon[1], 2, carbon[2]],
                                ["O", oxygen[1], 8, oxygen[2]], ["Cu", copper[1], 3, copper[2]]], dtype=object)
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
        molar_mass_pure = 3*copper[2] + 2*(carbon[2] + 3*oxygen[2]) + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.008, 5.844, 10.336], [92.333], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 112.35*10**9
        # Shear modulus
        G = 49.33*10**9
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
    def create_malachite(self):   # Cu2(CO3)(OH)2
        #
        name = "Mal"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["H", "C", "O", "Cu"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["C", carbon[1], 2, carbon[2]],
                                ["O", oxygen[1], 8, oxygen[2]], ["Cu", copper[1], 3, copper[2]]], dtype=object)
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
                minors = ["Zn", "Co", "Ni"]
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
        molar_mass_pure = 2*copper[2] + (carbon[2] + 3*oxygen[2]) + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.502, 11.974, 3.24], [98.75], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
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
    def create_ikaite(self):   # CaCO3 * 6*H2O
        #
        name = "Ika"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "C", "O", "Ca"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["C", carbon[1], 2, carbon[2]],
                                ["O", oxygen[1], 8, oxygen[2]], ["Ca", calcium[1], 3, calcium[2]]], dtype=object)
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
        molar_mass_pure = calcium[2] + (carbon[2] + 3*oxygen[2]) + 6*(2*hydrogen[2] + oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.87, 8.23, 11.02], [110.2], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = np.mean([30.9524, 40.3962])*10**9
        # Shear modulus
        G = np.mean([16.4415, 18.714])*10**9
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
class CustomCarbonates:
    #
    def __init__(self, fluid, actualThickness, output_type=False, porosity=None):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.output_type = output_type
        self.porosity = porosity
    #
    def create_custom_rock_01(self, amounts=None):
        #
        self.amounts = amounts
        #
        # Mineralogy + Fluids
        org = Carbonates(impurity="pure", dict=True).create_organic_matter()
        quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        calcite = Carbonates(impurity="pure", dict=True).create_calcite()
        dolomite = Carbonates(impurity="pure", dict=True).create_dolomite()
        siderite = Carbonates(impurity="pure", dict=True).create_siderite()
        anhydrite = Sulfates(data_type=True).create_anhydrite()
        gypsum = Sulfates(data_type=True).create_gypsum()
        pyrite = Sulfides(impurity="pure", dict=True).create_pyrite()
        illite = Pyllosilicates(impurity="pure", dict=True).create_illite()
        #
        mineralogy = [org, quartz, alkalifeldspar, plagioclase, calcite, dolomite, siderite, anhydrite, gypsum, pyrite,
                      illite]
        #
        water = fluids.Water.water("")
        #
        data = []
        results = {}
        #
        cond = False
        composition = []
        while cond == False:
            if self.amounts == None:
                w_org = round(0.1/100, 4)
                w_qz = round(1.3/100, 4)
                w_kfs = round(0/100, 4)
                w_pl = round(0/100, 4)
                w_cal = round(8.1/100, 4)
                w_dol = round(87/100, 4)
                w_sd = round(0/100, 4)
                w_anh = round(0/100, 4)
                w_gyp = round(0/100, 4)
                w_py = round(0.1/100, 4)
                w_ilt = round(1 - w_org - w_qz - w_kfs - w_pl - w_cal - w_dol - w_sd - w_anh - w_gyp - w_py, 4)
            elif type(self.amounts) is list:
                w_org = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_qz = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_kfs = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_pl = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_dol = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_sd = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_anh = round(abs(np.random.normal(self.amounts[7], 0.025)), 4)
                w_gyp = round(abs(np.random.normal(self.amounts[8], 0.025)), 4)
                w_py = round(abs(np.random.normal(self.amounts[9], 0.025)), 4)
                w_ilt = round(1 - w_org - w_qz - w_kfs - w_pl - w_cal - w_dol - w_sd - w_anh - w_gyp - w_py, 4)
            #
            if w_org >= 0.0 and w_qz >= 0.0 and w_kfs >= 0.0 and w_pl >= 0.0 and w_cal >= 0.0 and w_dol >= 0.0 \
                    and w_sd >= 0.0 and w_anh >= 0.0 and w_gyp >= 0.0 and w_py >= 0.0 and w_ilt >= 0.0:
                sumMin = round(w_org + w_qz + w_kfs + w_pl + w_cal + w_dol + w_sd + w_anh + w_gyp + w_py + w_ilt, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_org*org["chemistry"]["H"] + w_gyp*gypsum["chemistry"]["H"] + w_ilt*illite["chemistry"]["H"], 4)
            w_C = round(w_org*org["chemistry"]["C"] + w_cal*calcite["chemistry"]["C"] + w_dol*dolomite["chemistry"]["C"] + w_sd*siderite["chemistry"]["C"], 4)
            w_N = round(w_org*org["chemistry"]["N"], 4)
            w_Na = round(w_kfs*alkalifeldspar["chemistry"]["Na"] + w_pl*plagioclase["chemistry"]["Na"], 4)
            w_Mg = round(w_dol*dolomite["chemistry"]["Mg"] + w_ilt*illite["chemistry"]["Mg"], 4)
            w_Al = round(w_kfs*alkalifeldspar["chemistry"]["Al"] + w_pl*plagioclase["chemistry"]["Al"] + w_ilt*illite["chemistry"]["Al"], 4)
            w_Si = round(w_qz*quartz["chemistry"]["Si"] + w_kfs*alkalifeldspar["chemistry"]["Si"] + w_pl*plagioclase["chemistry"]["Si"] + w_ilt*illite["chemistry"]["Si"], 4)
            w_S = round(w_org*org["chemistry"]["S"] + w_py*pyrite["chemistry"]["S"], 4)
            w_K = round(w_kfs*alkalifeldspar["chemistry"]["K"] + w_ilt*illite["chemistry"]["K"], 4)
            w_Ca = round(w_pl*plagioclase["chemistry"]["Ca"] + w_cal*calcite["chemistry"]["Ca"] + w_dol*dolomite["chemistry"]["Ca"] + w_anh*anhydrite["chemistry"]["Ca"] + w_gyp*gypsum["chemistry"]["Ca"], 4)
            w_Fe = round(w_sd*siderite["chemistry"]["Fe"] + w_py*pyrite["chemistry"]["Fe"] + w_ilt*illite["chemistry"]["Fe"], 4)
            w_O = round(1 - w_H - w_C - w_N - w_Na - w_Mg - w_Al - w_Si - w_S - w_K - w_Ca - w_Fe, 4)
            sumConc = round(w_H + w_C + w_N + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Org", "Qz", "Kfs", "Pl", "Cal", "Dol", "Sd", "Anh", "Gyp", "Py", "Ilt"]))
                concentrations = [w_H, w_C, w_N, w_O, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe]
                amounts = [w_org, w_qz, w_kfs, w_pl, w_cal, w_dol, w_sd, w_anh, w_gyp, w_py, w_ilt]
            else:
                cond = False
        #
        element_list = ["H", "C", "N", "O", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Fe"]
        mineral_list = ["Org", "Qz", "Kfs", "Pl", "Cal", "Dol", "Sd", "Anh", "Gyp", "Py", "Ilt"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_org*org["rho"] + w_qz*quartz["rho"] + w_kfs*alkalifeldspar["rho"] + w_pl*plagioclase["rho"]
                    + w_cal*calcite["rho"] + w_dol*dolomite["rho"] + w_sd*siderite["rho"] + w_anh*anhydrite["rho"]
                    + w_gyp*gypsum["rho"] + w_py*pyrite["rho"] + w_ilt*illite["rho"]) / 1000
        rhoSolid = 0.975*rhoSolid
        X = [w_org, w_qz, w_kfs, w_pl, w_cal, w_dol, w_sd, w_anh, w_gyp, w_py, w_ilt]
        K_list = [mineralogy[i]["K"] for i in range(len(mineralogy))]
        G_list = [mineralogy[i]["G"] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = 0.35*K_geo
        G_solid = 0.275*G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if self.porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 4000:
                phi = rd.uniform(0.0, 0.025)
        else:
            phi = self.porosity
        #
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_org*org["GR"] + w_qz*quartz["GR"] + w_kfs*alkalifeldspar["GR"] + w_pl*plagioclase["GR"] \
             + w_cal*calcite["GR"] + w_dol*dolomite["GR"] + w_sd*siderite["GR"] + w_anh*anhydrite["GR"] \
             + w_gyp*gypsum["GR"] + w_py*pyrite["GR"] + w_ilt*illite["GR"]
        PE = w_org*org["PE"] + w_qz*quartz["PE"] + w_kfs*alkalifeldspar["PE"] + w_pl*plagioclase["PE"] \
             + w_cal*calcite["PE"] + w_dol*dolomite["PE"] + w_sd*siderite["PE"] + w_anh*anhydrite["PE"] \
             + w_gyp*gypsum["PE"] + w_py*pyrite["PE"] + w_ilt*illite["PE"]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_org*org["nu"] + w_qz*quartz["nu"] + w_kfs*alkalifeldspar["nu"] + w_pl*plagioclase["nu"] \
                                + w_cal*calcite["nu"] + w_dol*dolomite["nu"] + w_sd*siderite["nu"] + w_anh*anhydrite["nu"] \
                                + w_gyp*gypsum["nu"] + w_py*pyrite["nu"] + w_ilt*illite["nu"]
        #
        results["rho"] = round(rho*1000, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vP/vS, 4)
        results["G"] = round(G_bulk*10**(-6), 4)
        results["K"] = round(K_bulk*10**(-6), 4)
        results["E"] = round(E_bulk*10**(-6), 4)
        results["nu"] = round(poisson_mineralogical, 4)
        results["GR"] = round(GR, 4)
        results["PE"] = round(PE, 4)
        #
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append("water")
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        #
        if dict == False:
            return data
        else:
            return results