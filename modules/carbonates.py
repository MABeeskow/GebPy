#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		carbonates.py
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
from modules import fluids
from modules.geophysics import Elasticity as elast

class limestone:
    #
    def __init__(self, fluid, actualThickness):
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