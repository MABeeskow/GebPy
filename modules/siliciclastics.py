#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		siliciclastics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		21.02.2020

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals, oxides, fluids
from modules.geophysics import Elasticity as elast

class Soil:
    #
    def __init__(self):
        pass
    #
    def create_simple_soil(self, w_C=None, amounts=None, grainsize_list=False):
        self.w_C = w_C
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        illite = minerals.phyllosilicates.illite("")
        kaolinite = minerals.phyllosilicates.kaolinite("")
        organic = minerals.natives.organic_matter("")
        #
        mineralogy = [quartz, illite, kaolinite, organic]
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        air = fluids.Gas.air("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_C == None and self.amounts == None:
                w_org = 0.05
                w_qz = round(abs(rd.uniform(0.5, 0.95)), 4)
                w_ilt = round(abs(rd.uniform(0.0, (1-w_org-w_qz))), 4)
                w_kln = round(abs(1-w_qz-w_ilt-w_org), 4)
            elif self.w_C != None:
                w_org = round(self.w_C, 4)
                w_mineral = round(1-w_org, 4)
                w_qz = round(abs(w_mineral*rd.uniform(0.25, 1)), 4)
                w_ilt = round(abs(w_mineral*rd.uniform(0, (1-w_qz))), 4)
                w_kln = round(abs(w_mineral*(1-w_qz-w_ilt)), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_org = round(1-w_qz-w_ilt-w_kln, 4)
            #
            if 0.0 <= w_qz <= 1.0 and 0.0 <= w_ilt <= 1.0 and 0.0 <= w_kln <= 1.0 and 0.0 <= w_org <= 1.0:
                sumMin = round(w_qz + w_ilt + w_kln + w_org, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ilt*illite[6][0] + w_kln*kaolinite[6][0], 4)
            w_C = round(w_org, 4)
            w_O = round(w_qz*quartz[6][0] + w_ilt*illite[6][1] + w_kln*kaolinite[6][1], 4)
            w_Mg = round(w_ilt*illite[6][2], 4)
            w_Al = round(w_ilt*illite[6][3] + w_kln*kaolinite[6][2], 4)
            w_Si = round(w_qz*quartz[6][1] + w_ilt*illite[6][4] + w_kln*kaolinite[6][3], 4)
            w_K = round(w_ilt*illite[6][5], 4)
            w_Fe = round(w_ilt*illite[6][6], 4)
            sumConc = w_H + w_C + w_O + w_Mg + w_Al + w_Si + w_K + w_Fe
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Ilt", w_ilt, round(illite[1], 2)], ["Kln", w_kln, round(kaolinite[1], 2)], ["Org", w_org, round(organic[1], 2)]))
                composition.extend((["Qz", "Ilt", "Kln", "Org"]))
                concentrations = [w_H, w_C, w_O, w_Mg, w_Al, w_Si, w_K, w_Fe]
                amounts = [w_qz, w_ilt, w_kln, w_org]
            else:
                cond = False
        #
        grainsize = []
        n_Grains = 100
        if w_qz > 0:
            grainsize.append([rd.randint(2, 2000) for i in range(int(round(w_qz*n_Grains, 0)))])
        if w_ilt > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_ilt*n_Grains, 0)))], 2)))
        if w_kln > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_kln*n_Grains, 0)))], 2)))
        if w_org > 0:
            grainsize.append([rd.randint(2, 63) for i in range(int(round(w_org*n_Grains, 0)))])
        #
        rhoSolid = (w_qz*quartz[2] + w_ilt*illite[2] + w_kln*kaolinite[2] + w_org*organic[2])/1000
        X = [w_qz, w_ilt, w_kln, w_org]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        #vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        #vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        vP_solid = w_qz*quartz[4][0] + w_ilt*illite[4][0] + w_kln*kaolinite[4][0] + w_org*organic[4][0]
        vS_solid = w_qz*quartz[4][1] + w_ilt*illite[4][1] + w_kln*kaolinite[4][1] + w_org*organic[4][1]
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = rd.uniform(0.5, 0.65)
        #
        rho = (1 - phi) * rhoSolid + phi * (0.5*water[2]+0.5*air[2]/1000) / 1000
        vP = ((1-phi)*vP_solid + phi*(0.5*water[4][0] + 0.5*air[4][0]))/3
        vS = ((1 - phi) * vS_solid)/3
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - (0.5*water[2]+0.5*air[2]/1000) / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_ilt*illite[5][0] + w_kln*kaolinite[5][0] + w_org*organic[5][0]
        PE = w_qz*quartz[5][1] + w_ilt*illite[5][1] + w_kln*kaolinite[5][1] + w_org*organic[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_ilt*illite[3][3] + w_kln*kaolinite[3][3] + w_org*organic[3][3]
        #
        data.append(composition)
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6), round(air[2]/1000, 3)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2), round(air[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append(["water", "air"])
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        if grainsize_list == True:
            data.append(grainsize)
        #
        return data
    #
    def create_simple_sand(self, w_C=None, amounts=None, grainsize_list=False):
        self.w_C = w_C
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        quartz = minerals.oxides.quartz("")
        illite = minerals.phyllosilicates.illite("")
        kaolinite = minerals.phyllosilicates.kaolinite("")
        organic = minerals.natives.organic_matter("")
        #
        mineralogy = [quartz, illite, kaolinite, organic]
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        air = fluids.Gas.air("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_C == None and self.amounts == None:
                w_org = 0.025
                w_qz = round(abs(rd.uniform(0.85, 0.975)), 4)
                w_ilt = round(abs(rd.uniform(0.0, (1-w_org-w_qz))), 4)
                w_kln = round(abs(1-w_qz-w_ilt-w_org), 4)
            elif self.w_C != None:
                w_org = round(self.w_C, 4)
                w_mineral = round(1-w_org, 4)
                w_qz = round(abs(w_mineral*rd.uniform(0.25, 1)), 4)
                w_ilt = round(abs(w_mineral*rd.uniform(0, (1-w_qz))), 4)
                w_kln = round(abs(w_mineral*(1-w_qz-w_ilt)), 4)
            elif type(self.amounts) is list:
                w_qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_org = round(1-w_qz-w_ilt-w_kln, 4)
            #
            if 0.0 <= w_qz <= 1.0 and 0.0 <= w_ilt <= 1.0 and 0.0 <= w_kln <= 1.0 and 0.0 <= w_org <= 1.0:
                sumMin = round(w_qz + w_ilt + w_kln + w_org, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ilt*illite[6][0] + w_kln*kaolinite[6][0], 4)
            w_C = round(w_org, 4)
            w_O = round(w_qz*quartz[6][0] + w_ilt*illite[6][1] + w_kln*kaolinite[6][1], 4)
            w_Mg = round(w_ilt*illite[6][2], 4)
            w_Al = round(w_ilt*illite[6][3] + w_kln*kaolinite[6][2], 4)
            w_Si = round(w_qz*quartz[6][1] + w_ilt*illite[6][4] + w_kln*kaolinite[6][3], 4)
            w_K = round(w_ilt*illite[6][5], 4)
            w_Fe = round(w_ilt*illite[6][6], 4)
            sumConc = w_H + w_C + w_O + w_Mg + w_Al + w_Si + w_K + w_Fe
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                #composition.extend((["Qz", w_qz, round(quartz[1], 2)], ["Ilt", w_ilt, round(illite[1], 2)], ["Kln", w_kln, round(kaolinite[1], 2)], ["Org", w_org, round(organic[1], 2)]))
                composition.extend((["Qz", "Ilt", "Kln", "Org"]))
                concentrations = [w_H, w_C, w_O, w_Mg, w_Al, w_Si, w_K, w_Fe]
                amounts = [w_qz, w_ilt, w_kln, w_org]
            else:
                cond = False
        #
        grainsize = []
        n_Grains = 100
        if w_qz > 0:
            grainsize.append([rd.randint(2, 2000) for i in range(int(round(w_qz*n_Grains, 0)))])
        if w_ilt > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_ilt*n_Grains, 0)))], 2)))
        if w_kln > 0:
            grainsize.append(list(np.around([rd.uniform(1, 2) for i in range(int(round(w_kln*n_Grains, 0)))], 2)))
        if w_org > 0:
            grainsize.append([rd.randint(2, 63) for i in range(int(round(w_org*n_Grains, 0)))])
        #
        rhoSolid = (w_qz*quartz[2] + w_ilt*illite[2] + w_kln*kaolinite[2] + w_org*organic[2])/1000
        X = [w_qz, w_ilt, w_kln, w_org]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        #vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        #vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        vP_solid = w_qz*quartz[4][0] + w_ilt*illite[4][0] + w_kln*kaolinite[4][0] + w_org*organic[4][0]
        vS_solid = w_qz*quartz[4][1] + w_ilt*illite[4][1] + w_kln*kaolinite[4][1] + w_org*organic[4][1]
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = rd.uniform(0.45, 0.55)
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = ((1-phi)*vP_solid + phi*water[4][0])/3
        vS = ((1 - phi) * vS_solid)/3
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2]/1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_qz*quartz[5][0] + w_ilt*illite[5][0] + w_kln*kaolinite[5][0] + w_org*organic[5][0]
        PE = w_qz*quartz[5][1] + w_ilt*illite[5][1] + w_kln*kaolinite[5][1] + w_org*organic[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_qz*quartz[3][3] + w_ilt*illite[3][3] + w_kln*kaolinite[3][3] + w_org*organic[3][3]
        #
        data.append(composition)
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append("water")
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        if grainsize_list == True:
            data.append(grainsize)
        #
        return data
#
class sandstone:
    #
    def __init__(self, fluid, actualThickness):
        self.fluid = fluid
        self.actualThickness = actualThickness
    #
    def createSandstone(self):
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
        chem_quartz = minerals.oxides.quartz("")
        chem_dolomite = minerals.carbonates.dolomite("")
        chem_calcite = minerals.carbonates.calcite("")
        chemBiotite = minerals.Biotites.biotite_group(self, "Biotite")
        chemGlauconite = minerals.phyllosilicates.glauconite("")
        chemAlkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chemPlagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        #
        # [molar mass, density, bulk modulus, vP]
        chemWater = [18.0146, 997, 2.08, 1444]
        water = fluids.Water.water("")
        chemOil = [0.83 * chemC[2] + 0.11 * chemH[2] + 0.05 * chemS[2] + 0.005 * (chemO[2] + chemN[2]), 0.9, 1.35, 1225]
        chemGas = [0.83 * chemC[2] + 0.11 * chemH[2] + 0.05 * chemS[2] + 0.005 * (chemO[2] + chemN[2]), 0.8, 0.081, 475]
        #
        sandstone = []
        #
        cond = False
        composition = []
        while cond == False:
            magicnumber = rd.randint(0, 13)
            if magicnumber < 3:    # Quartz arenite
                xQuartz = rd.randint(90,100)/100                # Quartz
                xFeldspars = rd.randint(0,5)/100                # Feldspars
                xAlkalifeldspar2 = rd.randint(0,100)/100
                xPlagioclase2 = 1 - xAlkalifeldspar2
                xAlkalifeldspar = xFeldspars*xAlkalifeldspar2
                xPlagioclase = xFeldspars*xPlagioclase2
                xRockFragments = rd.randint(0,5)/100            # Rock Fragments
                xCalcite2 = rd.randint(0,100)/100
                xDolomite2 = 1-xCalcite2
                xCalcite = xRockFragments*xCalcite2
                xDolomite = xRockFragments*xDolomite2
                xAccessory = rd.randint(0,5)/100                # Accessory Minerals
                xBiotite2 = rd.randint(0,100)/100
                xGlauconite2 = 1-xBiotite2
                xBiotite = xAccessory*xBiotite2
                xGlauconite = xAccessory*xGlauconite2
            elif magicnumber > 2 and magicnumber < 6:    # Arkose
                xQuartz = rd.randint(0,75)/100                  # Quartz
                xFeldspars = rd.randint(25,100)/100             # Feldspars
                xAlkalifeldspar2 = rd.randint(0,100)/100
                xPlagioclase2 = 1 - xAlkalifeldspar2
                xAlkalifeldspar = xFeldspars*xAlkalifeldspar2
                xPlagioclase = xFeldspars*xPlagioclase2
                xRockFragments = rd.randint(0,10)/100           # Rock Fragments
                xCalcite2 = rd.randint(0,100)/100
                xDolomite2 = 1-xCalcite2
                xCalcite = xRockFragments*xCalcite2
                xDolomite = xRockFragments*xDolomite2
                xAccessory = rd.randint(0,5)/100                # Accessory Minerals
                xBiotite2 = rd.randint(0, 100) / 100
                xGlauconite2 = 1 - xBiotite2
                xBiotite = xAccessory*xBiotite2
                xGlauconite = xAccessory*xGlauconite2
            elif magicnumber == 6:    # Glauconite
                xQuartz = rd.randint(40,75)/100                 # Quartz
                xFeldspars = rd.randint(10,25)/100              # Feldspars
                xAlkalifeldspar2 = rd.randint(0,100)/100
                xPlagioclase2 = 1 - xAlkalifeldspar2
                xAlkalifeldspar = xFeldspars*xAlkalifeldspar2
                xPlagioclase = xFeldspars*xPlagioclase2
                xRockFragments = rd.randint(5,10)/100           # Rock Fragments
                xCalcite2 = rd.randint(0,100)/100
                xDolomite2 = 1-xCalcite2
                xCalcite = xRockFragments*xCalcite2
                xDolomite = xRockFragments*xDolomite2
                xAccessory = rd.randint(10,25)/100              # Accessory Minerals
                xGlauconite2 = rd.randint(95,100)/100
                xBiotite2 = 1-xGlauconite2
                xBiotite = xAccessory*xBiotite2
                xGlauconite = xAccessory*xGlauconite2
            elif magicnumber > 6 and magicnumber < 10:    # Calcic sandstone
                xQuartz = rd.randint(40,70)/100                  # Quartz
                xFeldspars = rd.randint(15,25)/100               # Feldspars
                xAlkalifeldspar2 = rd.randint(0,100)/100
                xPlagioclase2 = 1 - xAlkalifeldspar2
                xAlkalifeldspar = xFeldspars*xAlkalifeldspar2
                xPlagioclase = xFeldspars*xPlagioclase2
                xRockFragments = rd.randint(15,30)/100           # Rock Fragments
                xCalcite2 = rd.randint(90,100)/100
                xDolomite2 = 1-xCalcite2
                xCalcite = xRockFragments*xCalcite2
                xDolomite = xRockFragments*xDolomite2
                xAccessory = rd.randint(0,5)/100                 # Accessory Minerals
                xBiotite2 = rd.randint(0, 100) / 100
                xGlauconite2 = 1 - xBiotite2
                xBiotite = xAccessory*xBiotite2
                xGlauconite = xAccessory*xGlauconite2
            elif 9 < magicnumber < 12:    # Dolomitic sandstone
                xQuartz = rd.randint(40,70)/100                  # Quartz
                xFeldspars = rd.randint(15,25)/100               # Feldspars
                xAlkalifeldspar2 = rd.randint(0,100)/100
                xPlagioclase2 = 1 - xAlkalifeldspar2
                xAlkalifeldspar = xFeldspars*xAlkalifeldspar2
                xPlagioclase = xFeldspars*xPlagioclase2
                xRockFragments = rd.randint(15,30)/100           # Rock Fragments
                xDolomite2 = rd.randint(90,100)/100
                xCalcite2 = 1-xDolomite2
                xCalcite = xRockFragments*xCalcite2
                xDolomite = xRockFragments*xDolomite2
                xAccessory = rd.randint(0,5)/100                 # Accessory Minerals
                xBiotite2 = rd.randint(0, 100) / 100
                xGlauconite2 = 1 - xBiotite2
                xBiotite = xAccessory*xBiotite2
                xGlauconite = xAccessory*xGlauconite2
            elif 11 < magicnumber < 14:    # Subarkose
                xQuartz = rd.randint(50,100)/100                 # Quartz
                xFeldspars = rd.randint(0,25)/100                # Feldspars
                xAlkalifeldspar2 = rd.randint(0,100)/100
                xPlagioclase2 = 1 - xAlkalifeldspar2
                xAlkalifeldspar = xFeldspars*xAlkalifeldspar2
                xPlagioclase = xFeldspars*xPlagioclase2
                xRockFragments = rd.randint(0,25)/100            # Rock Fragments
                xCalcite2 = rd.randint(0,100)/100
                xDolomite2 = 1-xCalcite2
                xCalcite = xRockFragments*xCalcite2
                xDolomite = xRockFragments*xDolomite2
                xAccessory = rd.randint(0,5)/100                 # Accessory Minerals
                xBiotite2 = rd.randint(0, 100) / 100
                xGlauconite2 = 1 - xBiotite2
                xBiotite = xAccessory*xBiotite2
                xGlauconite = xAccessory*xGlauconite2
            sumMin = round(xQuartz,2) + round(xCalcite,2) + round(xDolomite,2) + round(xAlkalifeldspar,2) + round(xPlagioclase,2) + round(xBiotite,2) + round(xGlauconite,2)
            if sumMin == 1:
                cond = True
                composition.extend([["Qz", round(xQuartz,2), round(chem_quartz[1],2)], ["Cal", round(xCalcite,2), round(chem_calcite[1],2)], ["Dol", round(xDolomite,2), round(chem_dolomite[1],2)], ["Kfs", round(xAlkalifeldspar,2), round(chemAlkalifeldspar[1][0],2), chemAlkalifeldspar[1][1]], ["Pl", round(xPlagioclase,2), round(chemPlagioclase[1][0],2), chemPlagioclase[1][1]], ["Bt", round(xBiotite,2), round(chemBiotite[1][0],2), chemBiotite[1][1], chemBiotite[1][2]], ["Glt", round(xGlauconite,2), round(chemGlauconite[1],2)]])
            else:
                cond = False
        xQuartz = composition[0][1]
        xCalcite = composition[1][1]
        xDolomite = composition[2][1]
        xAlkalifeldspar = composition[3][1]
        xPlagioclase = composition[4][1]
        xBiotite = composition[5][1]
        xGlauconite = composition[6][1]
        sandstone.append(composition)
        mineralogy = [chem_quartz, chem_calcite, chem_dolomite, chemAlkalifeldspar, chemPlagioclase, chemBiotite, chemGlauconite]
        #
        rhoSolid = 1.1*(xQuartz*chem_quartz[2] + xCalcite*chem_calcite[2] + xDolomite*chem_dolomite[2] + xAlkalifeldspar*chemAlkalifeldspar[2] + xPlagioclase*chemPlagioclase[2] + xBiotite *chemBiotite[2] + xGlauconite *chemGlauconite[2]) / 1000
        X = [xQuartz, xCalcite, xDolomite, xAlkalifeldspar, xPlagioclase, xBiotite, xGlauconite]
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
        if self.actualThickness <= 1000:
            phi = randint(25, 35)/100
        elif self.actualThickness > 1000 and self.actualThickness <= 2000:
            phi = randint(20, 25)/100
        elif self.actualThickness > 2000 and self.actualThickness <= 3000:
            phi = randint(15, 20)/100
        elif self.actualThickness > 3000 and self.actualThickness <= 4000:
            phi = randint(10, 15)/100
        elif self.actualThickness > 4000:
            phi = randint(0, 10)/100
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
            GR = xQuartz*chem_quartz[5][0] + xCalcite*chem_calcite[5][0] + xDolomite*chem_dolomite[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xBiotite*chemBiotite[5][0] + xGlauconite*chemGlauconite[5][0]
            PE = xQuartz*chem_quartz[5][1] + xCalcite*chem_calcite[5][1] + xDolomite*chem_dolomite[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xBiotite*chemBiotite[5][1] + xGlauconite*chemGlauconite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xQuartz*chem_quartz[3][3] + xCalcite*chem_calcite[3][3] + xDolomite*chem_dolomite[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xBiotite*chemBiotite[3][3] + xGlauconite*chemGlauconite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            sandstone.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
            sandstone.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            sandstone.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[3], 2)])
            sandstone.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            sandstone.append("water")
            sandstone.append([GR, PE])
            sandstone.append(composition)
        elif self.fluid == "oil":
            rho = (1 - phi) * rhoSolid + phi * chemOil[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemOil[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemOil[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xQuartz*chem_quartz[5][0] + xCalcite*chem_calcite[5][0] + xDolomite*chem_dolomite[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xBiotite * chemBiotite[5][0] + xGlauconite * chemGlauconite[5][0]
            PE = xQuartz*chem_quartz[5][1] + xCalcite*chem_calcite[5][1] + xDolomite*chem_dolomite[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xBiotite*chemBiotite[5][1] + xGlauconite*chemGlauconite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xQuartz*chem_quartz[3][3] + xCalcite*chem_calcite[3][3] + xDolomite*chem_dolomite[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xBiotite*chemBiotite[3][3] + xGlauconite*chemGlauconite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            sandstone.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
            sandstone.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            sandstone.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemOil[3], 2)])
            sandstone.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            sandstone.append("oil")
            sandstone.append([GR, PE])
            sandstone.append(composition)
        elif self.fluid == "gas":
            rho = (1 - phi) * rhoSolid + phi * chemGas[1] / 1000
            vP = (1-phi)*vP_solid + phi*chemGas[3]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - chemGas[1] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = xQuartz*chem_quartz[5][0] + xCalcite*chem_calcite[5][0] + xDolomite*chem_dolomite[5][0] + xAlkalifeldspar*chemAlkalifeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xBiotite * chemBiotite[5][0] + xGlauconite * chemGlauconite[5][0]
            PE = xQuartz*chem_quartz[5][1] + xCalcite*chem_calcite[5][1] + xDolomite*chem_dolomite[5][1] + xAlkalifeldspar*chemAlkalifeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xBiotite*chemBiotite[5][1] + xGlauconite*chemGlauconite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = xQuartz*chem_quartz[3][3] + xCalcite*chem_calcite[3][3] + xDolomite*chem_dolomite[3][3] + xAlkalifeldspar*chemAlkalifeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xBiotite*chemBiotite[3][3] + xGlauconite*chemGlauconite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            sandstone.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[1] / 1000, 6)])
            sandstone.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            sandstone.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemGas[3], 2)])
            sandstone.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            sandstone.append("gas")
            sandstone.append([GR, PE])
            sandstone.append(composition)
        #
        #  sandstone = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return sandstone
    
    def create_simple_sandstone(self, w_Fe=None, amounts=None, porosity=None, pure=False):
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
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chem_quartz = minerals.oxides.quartz("")
        chem_alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chem_plagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        chem_calcite = minerals.carbonates.calcite("")
        chem_chlorite = minerals.phyllosilicates.chamosite("")
        chem_muscovite = minerals.phyllosilicates.muscovite("")
        chem_hematite = minerals.oxides.hematite("")
        #
        w_OQz = chem_quartz[6][0]
        w_SiQz = chem_quartz[6][1]
        w_OAfs = chem_alkalifeldspar[6][0]
        w_NaAfs = chem_alkalifeldspar[6][1]
        #w_AlAfs
        w_FeChl = chem_chlorite[6][5]
        w_FeHem = chem_hematite[6][1]
        #
        # [molar mass, density, bulk modulus, vP]
        air = fluids.Gas.air("")
        water = fluids.Water.water("")
        oil = fluids.Hydrocarbons.oil("")
        gas = fluids.Hydrocarbons.natural_gas("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Fe == None and self.amounts == None:
                magicnumber = rd.randint(0, 6)
                if magicnumber == 0 or pure == True:    # Qz-rich
                    w_ore = round(0.0, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
                    w_Qz = round(rd.randint(90, int((1-w_ore)*100))/100, 4)
                    w_Fsp = round(rd.randint(0, int((1-w_Qz)*100))/100, 4)
                    w_Afs2 = rd.randint(0, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_rf = round(1 - w_Qz - w_Fsp, 4)
                    magicnumber2 = rd.randint(0, 2)
                    if magicnumber2 == 0:
                        w_Cal2 = rd.randint(75, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 1:
                        magicnumber3 = rd.randint(0, 1)
                        if magicnumber3 == 0:
                            w_Chl2 = rd.randint(75, 100)/100
                            w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                        else:
                            w_Ms2 = rd.randint(75, 100)/100
                            w_Chl2 = rd.randint(0, int((1-w_Ms2)*100))/100
                        w_Cal2 = 1 - w_Chl2 - w_Ms2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 2:
                        w_Cal2 = rd.randint(0, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                elif magicnumber == 1:    # Afs-rich
                    w_ore = round(0.0, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
                    w_Fsp = round(rd.randint(25, int((1-w_ore)*75))/100, 4)
                    w_Afs2 = rd.randint(75, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Qz = round(rd.randint(25, int((1-w_ore-w_Fsp)*100))/100, 4)
                    w_rf = round(1 - w_Qz - w_Fsp - w_ore, 4)
                    magicnumber2 = rd.randint(0, 2)
                    if magicnumber2 == 0:
                        w_Cal2 = rd.randint(75, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 1:
                        magicnumber3 = rd.randint(0, 1)
                        if magicnumber3 == 0:
                            w_Chl2 = rd.randint(75, 100)/100
                            w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                        else:
                            w_Ms2 = rd.randint(75, 100)/100
                            w_Chl2 = rd.randint(0, int((1-w_Ms2)*100))/100
                        w_Cal2 = 1 - w_Chl2 - w_Ms2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 2:
                        w_Cal2 = rd.randint(0, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                elif magicnumber == 2:    # Pl-rich
                    w_ore = round(0.0, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
                    w_Fsp = round(rd.randint(25, int((1-w_ore)*75))/100, 4)
                    w_Pl2 = rd.randint(75, 100)/100
                    w_Afs2 = 1 - w_Pl2
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Qz = round(rd.randint(25, int((1-w_ore-w_Fsp)*100))/100, 4)
                    w_rf = round(1 - w_Qz - w_Fsp, 4)
                    magicnumber2 = rd.randint(0, 2)
                    if magicnumber2 == 0:
                        w_Cal2 = rd.randint(75, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 1:
                        magicnumber3 = rd.randint(0, 1)
                        if magicnumber3 == 0:
                            w_Chl2 = rd.randint(75, 100)/100
                            w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                        else:
                            w_Ms2 = rd.randint(75, 100)/100
                            w_Chl2 = rd.randint(0, int((1-w_Ms2)*100))/100
                        w_Cal2 = 1 - w_Chl2 - w_Ms2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 2:
                        w_Cal2 = rd.randint(0, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                elif magicnumber == 3:    # Cal-rich
                    w_ore = round(0.0, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
                    w_rf = round(rd.randint(25, 50)/100, 4)
                    w_Qz = round(rd.randint(25, int((1-w_rf)*100))/100, 4)
                    w_Fsp = round(1 - w_rf - w_Qz, 4)
                    w_Afs2 = rd.randint(0, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Cal2 = rd.randint(75, 100)/100
                    w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                    w_Ms2 = 1 - w_Cal2 - w_Chl2
                    w_Cal = round(w_rf*w_Cal2, 4)
                    w_Chl = round(w_rf*w_Chl2, 4)
                    w_Ms = round(w_rf*w_Ms2, 4)
                elif magicnumber == 4:    # Chl+Ms-rich
                    w_ore = round(0.0, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
                    w_rf = round(rd.randint(25, 50)/100, 4)
                    w_Qz = round(rd.randint(25, int((1-w_rf)*100))/100, 4)
                    w_Fsp = round(1 - w_rf - w_Qz, 4)
                    w_Afs2 = rd.randint(0, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Chl2 = rd.randint(50, 100)/100
                    w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                    w_Cal2 = 1 - w_Chl2 - w_Ms2
                    w_Cal = round(w_rf*w_Cal2, 4)
                    w_Chl = round(w_rf*w_Chl2, 4)
                    w_Ms = round(w_rf*w_Ms2, 4)
                elif magicnumber == 5:    # Hem-rich
                    w_ore = round(rd.randint(10, 20)/100, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
                    w_Qz = round(rd.randint(50, 80)/100, 4)
                    w_Fsp = round(rd.randint(0, int((1-w_ore-w_Qz)*100))/100, 4)
                    w_Afs2 = rd.randint(0, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_rf = round(1 - w_Qz - w_Fsp - w_ore, 4)
                    magicnumber2 = rd.randint(0, 2)
                    if magicnumber2 == 0:
                        w_Cal2 = rd.randint(75, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 1:
                        magicnumber3 = rd.randint(0, 1)
                        if magicnumber3 == 0:
                            w_Chl2 = rd.randint(75, 100)/100
                            w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                        else:
                            w_Ms2 = rd.randint(75, 100)/100
                            w_Chl2 = rd.randint(0, int((1-w_Ms2)*100))/100
                        w_Cal2 = 1 - w_Chl2 - w_Ms2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 2:
                        w_Cal2 = rd.randint(0, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                elif magicnumber == 6:    # intermediate
                    w_Qz = round(rd.randint(50, 75)/100, 4)
                    w_Fsp = round(rd.randint(0, 25)/100, 4)
                    w_Afs2 = rd.randint(0, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_rf = round(1 - w_Qz - w_Fsp, 4)
                    magicnumber2 = rd.randint(0, 2)
                    if magicnumber2 == 0:
                        w_Cal2 = rd.randint(75, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 1:
                        magicnumber3 = rd.randint(0, 1)
                        if magicnumber3 == 0:
                            w_Chl2 = rd.randint(75, 100)/100
                            w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                        else:
                            w_Ms2 = rd.randint(75, 100)/100
                            w_Chl2 = rd.randint(0, int((1-w_Ms2)*100))/100
                        w_Cal2 = 1 - w_Chl2 - w_Ms2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 2:
                        w_Cal2 = rd.randint(0, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    w_ore = round(0.0, 4)
                    w_Hem2 = rd.randint(0, 100)/100
                    w_Hem = round(w_ore*w_Hem2, 4)
            elif self.w_Fe != None:
                condition = False
                while condition == False:
                    w_rf = round(rd.randint(0, 10)/100, 4)
                    magicnumber2 = rd.randint(0, 2)
                    if magicnumber2 == 0:
                        w_Cal2 = rd.randint(75, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 1:
                        magicnumber3 = rd.randint(0, 1)
                        if magicnumber3 == 0:
                            w_Chl2 = rd.randint(75, 100)/100
                            w_Ms2 = rd.randint(0, int((1-w_Chl2)*100))/100
                        else:
                            w_Ms2 = rd.randint(75, 100)/100
                            w_Chl2 = rd.randint(0, int((1-w_Ms2)*100))/100
                        w_Cal2 = 1 - w_Chl2 - w_Ms2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    elif magicnumber2 == 2:
                        w_Cal2 = rd.randint(0, 100)/100
                        w_Chl2 = rd.randint(0, int((1-w_Cal2)*100))/100
                        w_Ms2 = 1 - w_Cal2 - w_Chl2
                        w_Cal = round(w_rf*w_Cal2, 4)
                        w_Chl = round(w_rf*w_Chl2, 4)
                        w_Ms = round(w_rf*w_Ms2, 4)
                    w_Hem = round((self.w_Fe - w_FeChl*w_Chl)/(w_FeHem), 4)
                    w_ore = round(w_Hem, 4)
                    if w_rf + w_ore <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_Qz = round(rd.randint(0, int((1-w_ore-w_rf)*100))/100, 4)
                w_Fsp = round(1 - w_ore - w_Qz - w_rf, 4)
                w_Afs2 = rd.randint(0, 100)/100
                w_Pl2 = 1 - w_Afs2
                w_Afs = round(w_Fsp*w_Afs2, 4)
                w_Pl = round(w_Fsp*w_Pl2, 4)
            elif type(self.amounts) is list:
                w_Qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_Afs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_Pl = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_Cal = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_Chl = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_Ms = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_Hem = round(1-w_Qz-w_Afs-w_Pl-w_Cal-w_Chl-w_Ms, 4)
            #
            if w_Qz >= 0.0 and w_Afs >= 0.0 and w_Pl >= 0.0 and w_Cal >= 0.0 and w_Chl >= 0.0 and w_Ms  >= 0.0 and w_Hem >= 0.0:
                sumMin = round(w_Qz + w_Afs + w_Pl + w_Cal + w_Chl + w_Ms + w_Hem, 4)
            else:
                sumMin = 0
            #
            w_H = round(chem_muscovite[6][0]*w_Ms + chem_chlorite[6][0]*w_Chl, 4)
            w_C = round(chem_calcite[6][0]*w_Cal, 4)
            w_O = round(chem_quartz[6][0]*w_Qz + chem_alkalifeldspar[6][0]*w_Afs + chem_plagioclase[6][0]*w_Pl + chem_calcite[6][1]*w_Cal + chem_chlorite[6][1]*w_Chl + chem_muscovite[6][1]*w_Ms + chem_hematite[6][0]*w_Hem, 4)
            w_F = round(chem_muscovite[6][2]*w_Ms, 4)
            w_Na = round(chem_alkalifeldspar[6][1]*w_Afs + chem_plagioclase[6][1]*w_Pl, 4)
            w_Mg = round(chem_chlorite[6][2]*w_Chl, 4)
            w_Al = round(chem_alkalifeldspar[6][2]*w_Afs + chem_plagioclase[6][2]*w_Pl + chem_chlorite[6][3]*w_Chl + chem_muscovite[6][3]*w_Ms, 4)
            w_Si = round(chem_quartz[6][1]*w_Qz + chem_alkalifeldspar[6][3]*w_Afs + chem_plagioclase[6][3]*w_Pl + chem_chlorite[6][4]*w_Chl + chem_muscovite[6][4]*w_Ms, 4)
            w_K = round(chem_alkalifeldspar[6][4]*w_Afs + chem_muscovite[6][5]*w_Ms, 4)
            w_Ca = round(chem_plagioclase[6][4]*w_Pl + chem_calcite[6][2]*w_Cal, 4)
            w_Fe_calc = round(chem_chlorite[6][5]*w_Chl + chem_hematite[6][1]*w_Hem, 4)
            sumConc = w_H + w_C + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_K + w_Ca + w_Fe_calc
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                #composition.extend((["Qz", w_Qz, round(chem_quartz[1], 2)], ["Kfs", w_Afs, round(chem_alkalifeldspar[1][0], 2), round(chem_alkalifeldspar[1][1], 2)], ["Pl", w_Pl, round(chem_plagioclase[1][0], 2), round(chem_plagioclase[1][1], 2)], ["Cal", w_Cal, round(chem_calcite[1], 2)], ["Chl", w_Chl, round(chem_chlorite[1], 2)], ["Ms", w_Ms, round(chem_muscovite[1], 2)], ["Hem", w_Hem, round(chem_hematite[1], 2)]))
                composition.extend((["Qz", "Kfs", "Pl", "Cal", "Chl", "Ms", "Hem"]))
                concentrations = [w_H, w_C, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_K, w_Ca, w_Fe_calc]
                amounts = [w_Qz, w_Afs, w_Pl, w_Cal, w_Chl, w_Ms, w_Hem]
            else:
                cond = False
        element_list = ["H", "C", "O", "F", "Na", "Mg", "Al", "Si", "K", "Ca", "Fe"]
        mineral_list = ["Qz", "Kfs", "Pl", "Cal", "Chl", "Ms", "Hem"]
        #data.append(composition)
        data.append([element_list, mineral_list])
        #
        mineralogy = [chem_quartz, chem_alkalifeldspar, chem_plagioclase, chem_calcite, chem_chlorite, chem_muscovite, chem_hematite]
        #
        rhoSolid = (w_Qz*chem_quartz[2] + w_Afs*chem_alkalifeldspar[2] + w_Pl*chem_plagioclase[2] + w_Cal*chem_calcite[2] + w_Chl*chem_chlorite[2] + w_Ms*chem_muscovite[2] + w_Hem*chem_hematite[2]) / 1000
        X = [w_Qz, w_Afs, w_Pl, w_Cal, w_Chl, w_Ms, w_Hem]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        #K_solid = K_geo/2
        #G_solid = G_geo/2
        K_solid = 0.5*K_geo/2
        G_solid = 0.5*G_geo/2
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if porosity == None:
            if self.actualThickness <= 1000 and self.w_Fe == None:
                phi = rd.uniform(0.25, 0.35)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000 and self.w_Fe == None:
                phi = rd.uniform(0.20, 0.25)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000 and self.w_Fe == None:
                phi = rd.uniform(0.15, 0.20)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000 and self.w_Fe == None:
                phi = rd.uniform(0.10, 0.15)
            elif self.actualThickness > 4000 or self.w_Fe != None:
                phi = rd.uniform(0.05, 0.10)
        else:
            phi = porosity
        #
        if self.fluid == "water" or self.w_Fe != None:
            rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
            vP = (1-phi)*vP_solid + phi*water[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0] + w_Cal*chem_calcite[5][0] + w_Chl*chem_chlorite[5][0] + w_Ms*chem_muscovite[5][0] + w_Hem*chem_hematite[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1] + w_Cal*chem_calcite[5][1] + w_Chl*chem_chlorite[5][1] + w_Ms*chem_muscovite[5][1] + w_Hem*chem_hematite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3] + w_Cal*chem_calcite[3][3] + w_Chl*chem_chlorite[3][3] + w_Ms*chem_muscovite[3][3] + w_Hem*chem_hematite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 3)])
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
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0] + w_Cal*chem_calcite[5][0] + w_Chl*chem_chlorite[5][0] + w_Ms*chem_muscovite[5][0] + w_Hem*chem_hematite[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1] + w_Cal*chem_calcite[5][1] + w_Chl*chem_chlorite[5][1] + w_Ms*chem_muscovite[5][1] + w_Hem*chem_hematite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3] + w_Cal*chem_calcite[3][3] + w_Chl*chem_chlorite[3][3] + w_Ms*chem_muscovite[3][3] + w_Hem*chem_hematite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(oil[2] / 1000, 3)])
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
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0] + w_Cal*chem_calcite[5][0] + w_Chl*chem_chlorite[5][0] + w_Ms*chem_muscovite[5][0] + w_Hem*chem_hematite[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1] + w_Cal*chem_calcite[5][1] + w_Chl*chem_chlorite[5][1] + w_Ms*chem_muscovite[5][1] + w_Hem*chem_hematite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3] + w_Cal*chem_calcite[3][3] + w_Chl*chem_chlorite[3][3] + w_Ms*chem_muscovite[3][3] + w_Hem*chem_hematite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(gas[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(gas[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("gas")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        elif self.fluid == "air":
            rho = (1 - phi) * rhoSolid + phi * air[2] / 1000
            vP = (1-phi)*vP_solid + phi*air[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - air[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0] + w_Cal*chem_calcite[5][0] + w_Chl*chem_chlorite[5][0] + w_Ms*chem_muscovite[5][0] + w_Hem*chem_hematite[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1] + w_Cal*chem_calcite[5][1] + w_Chl*chem_chlorite[5][1] + w_Ms*chem_muscovite[5][1] + w_Hem*chem_hematite[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3] + w_Cal*chem_calcite[3][3] + w_Chl*chem_chlorite[3][3] + w_Ms*chem_muscovite[3][3] + w_Hem*chem_hematite[3][3]
            #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(gas[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(gas[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("air")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        #
        return data
    #
    def create_feldspathic_sandstone(self, amounts=None, porosity=None, pure=False):
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chem_quartz = minerals.oxides.quartz("")
        chem_alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chem_plagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        #
        w_OQz = chem_quartz[6][0]
        w_SiQz = chem_quartz[6][1]
        w_OAfs = chem_alkalifeldspar[6][0]
        w_NaAfs = chem_alkalifeldspar[6][1]
        #
        # [molar mass, density, bulk modulus, vP]
        air = fluids.Gas.air("")
        water = fluids.Water.water("")
        oil = fluids.Hydrocarbons.oil("")
        gas = fluids.Hydrocarbons.natural_gas("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.amounts == None:
                magicnumber = rd.randint(0, 2)
                if magicnumber == 0 or pure == True:    # Qz-rich
                    w_Qz = round(rd.uniform(0.5, 1.0), 4)
                    w_Fsp = round(rd.uniform(0, 1-w_Qz), 4)
                    w_Afs2 = rd.uniform(0, 1)
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                elif magicnumber == 1:    # Afs-rich
                    w_Fsp = round(rd.uniform(0.25, 0.5), 4)
                    w_Afs2 = rd.uniform(0.5, 1)
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Qz = round(1-w_Fsp, 4)
                elif magicnumber == 2:    # Pl-rich
                    w_Fsp = round(rd.uniform(0.25, 0.5), 4)
                    w_Pl2 = rd.uniform(0.5, 1)
                    w_Afs2 = 1 - w_Pl2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Qz = round(1-w_Fsp, 4)
            elif type(self.amounts) is list:
                w_Qz = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_Afs = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_Pl = round(1-w_Qz-w_Afs, 4)
            #
            if w_Qz >= 0.0 and w_Afs >= 0.0 and w_Pl >= 0.0:
                sumMin = round(w_Qz + w_Afs + w_Pl, 4)
            else:
                sumMin = 0
            #
            w_O = round(w_Qz*chem_quartz[6][0] + w_Afs*chem_alkalifeldspar[6][0] + w_Pl*chem_plagioclase[6][0], 4)
            w_Na = round(w_Afs*chem_alkalifeldspar[6][1] + w_Pl*chem_plagioclase[6][1], 4)
            w_Al = round(w_Afs*chem_alkalifeldspar[6][2] + w_Pl*chem_plagioclase[6][2], 4)
            w_Si = round(w_Qz*chem_quartz[6][1] + w_Afs*chem_alkalifeldspar[6][3] + w_Pl*chem_plagioclase[6][3], 4)
            w_K = round(w_Afs*chem_alkalifeldspar[6][4], 4)
            w_Ca = round(w_Pl*chem_plagioclase[6][4], 4)
            sumConc = w_O + w_Na + w_Al + w_Si + w_K + w_Ca
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Qz", "Kfs", "Pl"]))
                concentrations = [w_O, w_Na, w_Al, w_Si, w_K, w_Ca]
                amounts = [w_Qz, w_Afs, w_Pl]
            else:
                cond = False
        element_list = ["O", "Na", "Al", "Si", "K", "Ca"]
        mineral_list = ["Qz", "Kfs", "Pl"]
        data.append([element_list, mineral_list])
        #
        mineralogy = [chem_quartz, chem_alkalifeldspar, chem_plagioclase]
        #
        rhoSolid = (w_Qz*chem_quartz[2] + w_Afs*chem_alkalifeldspar[2] + w_Pl*chem_plagioclase[2]) / 1000
        X = [w_Qz, w_Afs, w_Pl]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = 0.5*K_geo/2
        G_solid = 0.5*G_geo/2
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        if porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.uniform(0.25, 0.35)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.uniform(0.20, 0.25)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.uniform(0.15, 0.20)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.uniform(0.10, 0.15)
            elif self.actualThickness > 4000:
                phi = rd.uniform(0.05, 0.10)
        else:
            phi = porosity
        #
        if self.fluid == "water":
            w_fluid = round((phi)/(1-phi)*(w_Qz/chem_quartz[2] + w_Afs/chem_alkalifeldspar[2] + w_Pl/chem_plagioclase[2])*water[2], 4)
            rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
            vP = (1-phi)*vP_solid + phi*water[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3]
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(["water", w_fluid])
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        elif self.fluid == "oil":
            w_fluid = round((phi)/(1-phi)*(w_Qz/chem_quartz[2] + w_Afs/chem_alkalifeldspar[2] + w_Pl/chem_plagioclase[2])*oil[2], 4)
            rho = (1 - phi) * rhoSolid + phi * oil[2] / 1000
            vP = (1-phi)*vP_solid + phi*oil[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - oil[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3]
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(oil[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(oil[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(["oil", w_fluid])
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        elif self.fluid == "gas":
            w_fluid = round((phi)/(1-phi)*(w_Qz/chem_quartz[2] + w_Afs/chem_alkalifeldspar[2] + w_Pl/chem_plagioclase[2])*gas[2], 4)
            rho = (1 - phi) * rhoSolid + phi * gas[2] / 1000
            vP = (1-phi)*vP_solid + phi*gas[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - gas[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3]
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(gas[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(gas[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(["gas", w_fluid])
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        elif self.fluid == "air":
            w_fluid = round((phi)/(1-phi)*(w_Qz/chem_quartz[2] + w_Afs/chem_alkalifeldspar[2] + w_Pl/chem_plagioclase[2])*air[2], 4)
            rho = (1 - phi) * rhoSolid + phi * air[2] / 1000
            vP = (1-phi)*vP_solid + phi*air[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - air[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            GR = w_Qz*chem_quartz[5][0] + w_Afs*chem_alkalifeldspar[5][0] + w_Pl*chem_plagioclase[5][0]
            PE = w_Qz*chem_quartz[5][1] + w_Afs*chem_alkalifeldspar[5][1] + w_Pl*chem_plagioclase[5][1]
            poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
            poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
            poisson_mineralogical = w_Qz*chem_quartz[3][3] + w_Afs*chem_alkalifeldspar[3][3] + w_Pl*chem_plagioclase[3][3]
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(gas[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(gas[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(["air", w_fluid])
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
        #
        return data
    #
class shale:
    #
    def __init__(self):
        pass
    #
    def create_shale(self, w_C=None, w_Na=None, w_Mg=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        # Parameters
        self.w_C = w_C
        self.w_Na = w_Na
        self.w_Mg = w_Mg
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # Minerals + Fluids
        org = minerals.Organics.organic_matter("")
        qz = minerals.oxides.quartz("")
        cal = minerals.carbonates.calcite("")
        ilt = minerals.phyllosilicates.illite("")
        kln = minerals.phyllosilicates.kaolinite("")
        mnt = minerals.phyllosilicates.montmorillonite("")
        water = fluids.Water.water("")
        #
        mineralogy = [org, qz, cal, ilt, kln, mnt]
        #
        w_CCal = cal[6][0]
        w_COrg = org[6][1]
        w_NaMnt = mnt[6][2]
        w_MgIlt = ilt[6][2]
        w_MgMnt = mnt[6][3]
        w_KIlt = ilt[6][5]
        w_CaCal = cal[6][2]
        w_CaMnt = mnt[6][6]
        w_FeIlt = ilt[6][6]
        #
        data = []
        #
        condition = False
        composition = []
        while condition == False:
            if self.w_C == None and self.w_Na == None and self.w_Mg == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                #print("Default")
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ilt2 = rd.uniform(0.5, 1.0)
                w_kln2 = rd.uniform(0.0, float(1.0-w_ilt2))
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(rd.uniform(0.0, float(1.0-w_clay)), 4)
                w_cal = round(rd.uniform(0.0, float(1.0-w_clay-w_qz)), 4)
                w_org = round(1-w_clay-w_qz-w_cal, 4)
            elif self.w_C != None:
                #print("w_C != None")
                w_cal = round(rd.uniform(0.0, 0.15), 4)
                w_org = round((w_C - w_cal*w_CCal)/w_COrg, 4)
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ilt2 = rd.uniform(0.5, 1.0)
                w_kln2 = rd.uniform(0.0, float(1.0-w_ilt2))
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(1-w_cal-w_org-w_clay, 4)
            elif self.w_Na != None:
                #print("w_Na != None")
                w_mnt = round(w_Na/w_NaMnt, 4)
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_mnt2 = w_mnt/w_clay
                w_ilt2 = rd.uniform(0.0, float(1-w_mnt2))
                w_kln2 = 1-w_ilt2-w_mnt2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_qz = round(rd.uniform(0.2, float(1.0-w_clay)), 4)
                w_cal = round(rd.uniform(0.15, float(1.0-w_clay-w_qz)), 4)
                w_org = round(1-w_clay-w_qz-w_cal, 4)
            elif self.w_Mg != None:
                #print("w_Mg != None")
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_mnt2 = rd.uniform(0.0, 1.0)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_ilt = round((w_Mg - w_mnt*w_MgMnt)/w_MgIlt, 4)
                w_ilt2 = w_ilt/w_clay
                w_kln2 = 1-w_mnt2-w_ilt2
                w_kln = round(w_clay*w_kln2, 4)
                w_qz = round(rd.uniform(0.2, float(1.0-w_clay)), 4)
                w_cal = round(rd.uniform(0.15, float(1.0-w_clay-w_qz)), 4)
                w_org = round(1-w_clay-w_qz-w_cal, 4)
            elif self.w_K != None:
                #print("w_K != None")
                w_ilt = round(w_K/w_KIlt, 4)
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ilt2 = w_ilt/w_clay
                w_kln2 = rd.uniform(0.0, float(1.0-w_ilt2))
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(rd.uniform(0.2, float(1.0-w_clay)), 4)
                w_cal = round(rd.uniform(0.15, float(1.0-w_clay-w_qz)), 4)
                w_org = round(1-w_clay-w_qz-w_cal, 4)
            elif self.w_Ca != None:
                #print("w_Ca != None")
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ilt2 = rd.uniform(0.5, 1.0)
                w_kln2 = rd.uniform(0.0, float(1.0-w_ilt2))
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_cal = round((w_Ca - w_mnt*w_CaMnt)/w_CaCal, 4)
                w_qz = round(rd.uniform(0.2, float(1.0-w_clay-w_cal)), 4)
                w_org = round(1-w_clay-w_qz-w_cal, 4)
            elif self.w_Fe != None:
                #print("w_Fe != None")
                w_ilt = round(w_Fe/w_FeIlt, 4)
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ilt2 = w_ilt/w_clay
                w_kln2 = rd.uniform(0.0, float(1.0-w_ilt2))
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(rd.uniform(0.2, float(1.0-w_clay)), 4)
                w_cal = round(rd.uniform(0.15, float(1.0-w_clay-w_qz)), 4)
                w_org = round(1-w_clay-w_qz-w_cal, 4)
            elif type(self.amounts) is list:
                #print("Amount")
                w_org = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_qz = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_mnt = round(1-w_org-w_qz-w_cal-w_ilt-w_kln, 4)
            #
            #print("Minerals:", "Org", w_org, "Qz", w_qz, "Cal", w_cal, "Ilt", w_ilt, "Kln", w_kln, "Mnt", w_mnt, "Clay:", round(w_ilt+w_kln+w_mnt, 4) ,"Sum:", round(w_org + w_qz + w_cal + w_ilt + w_kln + w_mnt, 4))
            if w_org >= 0.0 and w_qz >= 0.0 and w_cal >= 0.0 and w_ilt >= 0.0 and w_kln >= 0.0 and w_mnt >= 0.0:
                sum_min = round(w_org + w_qz + w_cal + w_ilt + w_kln + w_mnt, 4)
                if sum_min == 1:
                    #
                    w_H = round(w_ilt*ilt[6][0] + w_kln*kln[6][0] + w_mnt*mnt[6][0] + w_org*org[6][0], 4)
                    w_C = round(w_org*org[6][1] + w_cal*cal[6][0], 4)
                    w_N = round(w_org*org[6][2], 4)
                    w_O = round(w_org*org[6][3] + w_qz*qz[6][0] + w_cal*cal[6][1] + w_ilt*ilt[6][1] + w_kln*kln[6][1] + w_mnt*mnt[6][1], 4)
                    w_Na = round(w_mnt*mnt[6][2], 4)
                    w_Mg = round(w_ilt*ilt[6][2] + w_mnt*mnt[6][3], 4)
                    w_Al = round(w_ilt*ilt[6][3] + w_kln*kln[6][2] + w_mnt*mnt[6][4], 4)
                    w_Si = round(w_qz*qz[6][1] + w_ilt*ilt[6][4] + w_kln*kln[6][3] + w_mnt*mnt[6][5], 4)
                    w_S = round(w_org*org[6][4], 4)
                    w_K = round(w_ilt*ilt[6][5], 4)
                    w_Ca = round(w_cal*cal[6][2] + w_mnt*mnt[6][6], 4)
                    w_Fe = round(w_ilt*ilt[6][6], 4)
                    sum_conc = round(w_H + w_C + w_N + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe, 4)
                    #
                    GR = w_org*org[5][0] + w_qz*qz[5][0] + w_cal*cal[5][0] + w_ilt*ilt[5][0] + w_kln*kln[5][0] + w_mnt*mnt[5][0]
                    #print("sum conc:", sum_conc, "GR", GR)
                    #
                    if sum_conc == 1 and GR <= 300:
                        #condition = True
                        composition.extend((["Org", "Qz", "Cal", "Ilt", "Kln", "Mnt"]))
                        concentrations = [abs(w_H), abs(w_C), abs(w_N), abs(w_O), abs(w_Na), abs(w_Mg), abs(w_Al), abs(w_Si), abs(w_S), abs(w_K), abs(w_Ca), abs(w_Fe)]
                        amounts = [abs(w_org), abs(w_qz), abs(w_cal), abs(w_ilt), abs(w_kln), abs(w_mnt)]
                        rhoSolid = (w_org*org[2] + w_qz*qz[2] + w_cal*cal[2] + w_ilt*ilt[2] + w_kln*kln[2] + w_mnt*mnt[2])/1000
                        condition_2 = False
                        i = 0
                        while condition_2 == False and i < 10:
                            i += 1
                            phi = rd.uniform(0.0, 0.1)
                            rho = (1 - phi)*rhoSolid + phi*water[2]/1000
                            if rho > 2.0:
                                condition_2 = True
                                condition = True
                            else:
                                condition_2 = False
                    else:
                        condition = False
            else:
                condition = False
        #
        X = [w_org, w_qz, w_cal, w_ilt, w_kln, w_mnt]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo/6
        G_solid = G_geo/6
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        #GR = w_org*org[5][0] + w_qz*qz[5][0] + w_cal*cal[5][0] + w_ilt*ilt[5][0] + w_kln*kln[5][0] + w_mnt*mnt[5][0]
        PE = w_org*org[5][1] + w_qz*qz[5][1] + w_cal*cal[5][1]  + w_ilt*ilt[5][1] + w_kln*kln[5][1] + w_mnt*mnt[5][1]
        #poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        #poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_org*org[3][3] + w_qz*qz[3][3] + w_cal*cal[3][3] + w_ilt*ilt[3][3] + w_kln*kln[3][3] + w_mnt*mnt[3][3]
        #
        data.append(composition)
        data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
        data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
        data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        data.append("water")
        data.append([round(GR, 3), round(PE, 3)])
        data.append(concentrations)
        data.append(amounts)
        #
        return data
    #
    def create_simple_shale(self, w_C=None, w_F=None, w_Na=None, w_Mg=None, w_S=None, w_K=None, w_Ca=None, w_Fe=None, amounts=None):
        #
        self.w_C = w_C
        self.w_F = w_F
        self.w_Na = w_Na
        self.w_Mg = w_Mg
        self.w_S = w_S
        self.w_K = w_K
        self.w_Ca = w_Ca
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        # mineralogy
        chem_org = minerals.natives.organic_matter("")
        chem_qz = minerals.oxides.quartz("")
        chem_ms = minerals.phyllosilicates.muscovite("")
        chem_bt = minerals.Biotites.biotite_group(self, "simple")
        chem_cal = minerals.carbonates.calcite("")
        chem_ilt = minerals.phyllosilicates.illite("")
        chem_kln = minerals.phyllosilicates.kaolinite("")
        chem_mnt = minerals.phyllosilicates.montmorillonite("")
        chem_py = minerals.sulfides.pyrite("")
        chem_urn = minerals.oxides.uraninite("")
        #
        mineralogy = [chem_org, chem_qz, chem_cal, chem_py, chem_ilt, chem_kln, chem_mnt, chem_bt, chem_ms, chem_urn]
        #for i in range(len(mineralogy)):
        #    print(mineralogy[i][0], mineralogy[i][6])
        #
        w_CCal = chem_cal[6][0]
        w_COrg = chem_org[6][0]
        w_FBt = chem_bt[6][2]
        w_FMs = chem_ms[6][2]
        w_NaMnt = chem_mnt[6][2]
        w_MgBt = chem_bt[6][3]
        w_MgIlt = chem_ilt[6][2]
        w_MgMnt = chem_mnt[6][3]
        w_SPy = chem_py[6][0]
        w_KMs = chem_ms[6][5]
        w_KBt = chem_bt[6][6]
        w_KIlt = chem_ilt[6][5]
        w_CaCal = chem_cal[6][2]
        w_CaMnt = chem_mnt[6][6]
        w_FePy = chem_py[6][1]
        w_FeBt = chem_bt[6][7]
        w_FeIlt = chem_ilt[6][6]
        #
        # [molar mass, density, bulk modulus, vP]
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_C == None and self.w_F == None and self.w_Na == None and self.w_Mg == None and self.w_S == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                magicnumber = rd.randint(0, 2)
                #magicnumber = 0
                if magicnumber == 0:    # Clay-rich
                    w_ore = round(rd.randint(0, 5)/100, 4)
                    w_Py = w_ore
                    w_clay = round(rd.randint(45, int((0.70-w_ore)*100))/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_mnt2 = 1-w_ilt2-w_kln2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_qz = round(rd.randint(0, int((1-w_ore-w_clay)*100))/100, 4)
                    w_mica = round(rd.randint(0, int((1-w_ore-w_clay-w_qz)*100))/100, 4)
                    w_ms2 = rd.randint(0, 100)/100
                    w_bt2 = 1 - w_ms2
                    w_ms = round(w_mica*w_ms2, 4)
                    w_bt = round(w_mica*w_bt2, 4)
                    w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                    w_urn = round(rd.randint(0, 5)/3000000, 6)
                    w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
                elif magicnumber == 1:    # Qz-rich
                    w_ore = round(rd.randint(0, 5)/100, 4)
                    w_Py = w_ore
                    w_qz = round(rd.randint(25, 50)/100, 4)
                    w_clay = round(rd.randint(40, int((1-w_ore-w_qz)*100))/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_mnt2 = 1-w_ilt2-w_kln2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_mica = round(rd.randint(0, int((1-w_ore-w_clay-w_qz)*100))/100, 4)
                    w_ms2 = rd.randint(0, 100)/100
                    w_bt2 = 1 - w_ms2
                    w_ms = round(w_mica*w_ms2, 4)
                    w_bt = round(w_mica*w_bt2, 4)
                    w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                    w_urn = round(rd.randint(0, 5)/3000000, 6)
                    w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
                elif magicnumber == 2:    # Mica-rich
                    w_ore = round(rd.randint(0, 5)/100, 4)
                    w_Py = w_ore
                    w_mica = round(rd.randint(15, 30)/100, 4)
                    w_ms2 = rd.randint(0, 100)/100
                    w_bt2 = 1 - w_ms2
                    w_ms = round(w_mica*w_ms2, 4)
                    w_bt = round(w_mica*w_bt2, 4)
                    w_clay = round(rd.randint(50, int((1-w_ore-w_mica)*100))/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_mnt2 = 1-w_ilt2-w_kln2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_mica)*100))/100, 4)
                    w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                    w_urn = round(rd.randint(0, 5)/3000000, 6)
                    w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
            elif self.w_C != None:
                condition = False
                while condition == False:
                    w_cal = round(rd.randint(0, 15)/100, 4)
                    w_org = round((self.w_C - w_CCal*w_cal)/(w_COrg), 4)
                    if w_cal + w_org <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_ore = round(rd.randint(1, int((0.25-w_cal-w_org)*100))/100, 4)
                w_Py = w_ore
                w_clay = round(rd.randint(40, int((1-w_ore-w_org-w_cal)*100))/100, 4)
                w_ilt2 = rd.randint(0, 100)/100
                w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_org-w_cal)*100))/100, 4)
                w_urn = round(rd.randint(0, 5)/3000000, 6)
                w_mica = round(1-w_ore-w_clay-w_qz-w_org-w_cal-w_urn, 4)
                w_ms2 = rd.randint(0, 100)/100
                w_bt2 = 1 - w_ms2
                w_ms = round(w_mica*w_ms2, 4)
                w_bt = round(w_mica*w_bt2, 4)
            elif self.w_F != None:
                condition = False
                while condition == False:
                    w_ms = round(rd.randint(5, 20)/100, 4)
                    w_bt = round((self.w_F - w_FMs*w_ms)/(w_FBt), 4)
                    w_mica = w_bt + w_ms
                    if w_ms + w_bt <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_ore = round(rd.randint(0, 5)/100, 4)
                w_Py = w_ore
                w_clay = round(rd.randint(50, int((1-w_ore-w_mica)*100))/100, 4)
                w_ilt2 = rd.randint(50, 100)/100
                w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_mica)*100))/100, 4)
                w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                w_urn = round(rd.randint(0, 5)/3000000, 6)
                w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
            elif self.w_Na != None:
                condition = False
                while condition == False:
                    w_mnt = round((self.w_Na)/(w_NaMnt), 4)
                    if w_mnt <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_clay = round(rd.randint(int(w_mnt*100)+1, 100)/100, 4)
                w_mnt2 = w_mnt/w_clay
                w_ilt2 = rd.randint(0, int((1-w_mnt2)*100))/100
                w_kln2 = rd.randint(0, int((1-w_ilt2-w_mnt2)*100))/100
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_qz = round(rd.randint(0, int((1-w_clay)*100))/100, 4)
                w_mica = round(rd.randint(0, int((1-w_clay-w_qz)*100))/100, 4)
                w_ms2 = rd.randint(0, 100)/100
                w_bt2 = 1 - w_ms2
                w_ms = round(w_mica*w_ms2, 4)
                w_bt = round(w_mica*w_bt2, 4)
                w_cal = round(rd.randint(0, int((1-w_clay-w_qz-w_mica)*100))/100, 4)
                w_org = round(1-w_clay-w_qz-w_mica-w_cal, 4)
                w_urn = round(rd.randint(0, 5)/3000000, 6)
                w_ore = round(rd.randint(0, int((1-w_clay-w_qz-w_mica-w_org-w_urn)*100))/100, 4)
                w_Py = w_ore
            elif self.w_Mg != None:
                condition = False
                while condition == False:
                    w_clay = round(rd.randint(50, 100)/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_mnt2 = 1-w_ilt2-w_kln2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_bt = round((self.w_Mg - w_MgIlt*w_ilt - w_MgMnt*w_mnt)/(w_MgBt), 4)
                    if w_clay + w_bt <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_ore = round(rd.randint(0, 5)/100, 4)
                w_Py = w_ore
                w_qz = round(rd.randint(0, int((1-w_ore-w_clay)*100))/100, 4)
                w_mica = round(rd.randint(0, int((1-w_ore-w_clay-w_qz)*100))/100, 4)
                w_bt2 = w_bt/w_mica
                w_ms2 = 1-w_bt2
                w_ms = round(w_mica*w_ms2, 4)
                w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                w_urn = round(rd.randint(0, 5)/3000000, 6)
                w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
            elif self.w_S != None:
                condition = False
                while condition == False:
                    w_Py = round((self.w_S)/(w_SPy), 4)
                    if w_Py <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_ore = w_Py
                w_org = round(rd.randint(2, 5)/100, 4)
                w_clay = round(rd.randint(50, int((1-w_ore-w_org)*100))/100, 4)
                w_ilt2 = rd.randint(50, 100)/100
                w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                w_mnt2 = 1-w_ilt2-w_kln2
                w_ilt = round(w_clay*w_ilt2, 4)
                w_kln = round(w_clay*w_kln2, 4)
                w_mnt = round(w_clay*w_mnt2, 4)
                w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_org)*100))/100, 4)
                w_mica = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_org)*100))/100, 4)
                w_ms2 = rd.randint(0, 100)/100
                w_bt2 = 1 - w_ms2
                w_ms = round(w_mica*w_ms2, 4)
                w_bt = round(w_mica*w_bt2, 4)
                w_urn = round(rd.randint(0, 5)/3000000, 6)
                w_cal = round(1-w_ore-w_clay-w_qz-w_mica-w_org-w_urn, 4)
            elif self.w_K != None:
                condition = False
                while condition == False:
                    w_ore = round(rd.randint(0, 5)/100, 4)
                    w_Py = w_ore
                    w_mica = round(rd.randint(15, 30)/100, 4)
                    w_bt2 = rd.randint(0, 100)/100
                    w_bt = round(w_mica*w_bt2, 4)
                    w_clay = round(rd.randint(50, int((1-w_ore-w_mica)*100))/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_mnt2 = 1-w_ilt2-w_kln2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_ms = round((self.w_K - w_KIlt*w_ilt - w_KBt*w_bt)/(w_KMs), 4)
                    if w_ore + w_bt + w_ms + w_clay <= 1.0:
                        condition = True
                    else:
                        condition = False
                w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_mica)*100))/100, 4)
                w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                w_urn = round(rd.randint(0, 5)/3000000, 6)
                w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
            elif self.w_Ca != None:
                condition_1 = False
                while condition_1 == False:
                    w_ore = round(rd.randint(0, 10)/100, 4)
                    w_Py = w_ore
                    w_clay = round(rd.randint(40, int((0.60-w_ore)*100))/100, 4)
                    w_mnt2 = rd.randint(0, 100)/100
                    w_ilt2 = rd.randint(0, int((1-w_mnt2)*100))/100
                    w_kln2 = 1-w_ilt2-w_mnt2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_cal = round((self.w_Ca - w_CaMnt*w_mnt)/(w_CaCal), 4)
                    if w_ore + w_clay + w_cal <= 1.0:
                        condition_1 = True
                    else:
                        condition_1 = False
                condition_2 = False
                while condition_2 == False:
                    w_qz = round(rd.randint(5, int((0.9-w_ore-w_clay-w_cal)*100))/100, 4)
                    w_mica = round(rd.randint(5, int((1-w_ore-w_clay-w_qz-w_cal)*100))/100, 4)
                    w_ms2 = rd.randint(0, 100)/100
                    w_bt2 = 1 - w_ms2
                    w_ms = round(w_mica*w_ms2, 4)
                    w_bt = round(w_mica*w_bt2, 4)
                    w_urn = round(rd.randint(0, 5)/3000000, 6)
                    w_org = round(1-w_ore-w_clay-w_qz-w_mica-w_cal-w_urn, 4)
                    if w_org >= 0:
                        condition_2 = True
                    else:
                        condition_2 = False
            elif self.w_Fe != None:
                condition = False
                while condition == False:
                    w_org = round(rd.randint(2, 8)/100, 4)
                    w_clay = round(rd.randint(35, 70)/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_kln2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_mnt2 = 1-w_ilt2-w_kln2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_mica = round(rd.randint(0, int((1-w_clay-w_org)*100))/100, 4)
                    w_ms2 = rd.randint(0, 100)/100
                    w_bt2 = 1 - w_ms2
                    w_ms = round(w_mica*w_ms2, 4)
                    w_bt = round(w_mica*w_bt2, 4)
                    w_Py = round((self.w_Fe - w_FeIlt*w_ilt - w_FeBt*w_bt)/(w_FePy), 4)
                    w_ore = w_Py
                    if w_org + w_clay + w_mica + w_ore <= 1.0:
                        condition = True
                    else:
                        condition = False
                condition = False
                while condition == False:
                    w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_org)*100))/100, 4)
                    w_urn = round(rd.randint(0, 5)/3000000, 6)
                    w_cal = round(1-w_ore-w_clay-w_qz-w_mica-w_org-w_urn, 4)
                    if w_qz >= 0 and w_cal >= 0:
                        condition = True
                    else:
                        condition = False
            elif type(self.amounts) is list:
                w_urn = round(abs(np.random.normal(self.amounts[9], 5e-07)), 6)
                w_org = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_qz = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_Py = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_mnt = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[7], 0.025)), 4)
                w_ms = round(1-w_org-w_qz-w_cal-w_ilt-w_kln-w_mnt-w_bt-w_urn, 4)
            #
            if w_org >= 0.0 and w_qz >= 0.0 and w_cal >= 0.0 and w_Py >= 0.0 and w_ilt >= 0.0 and w_kln >= 0.0 and w_mnt >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_urn >= 0.0:
                sumMin = round(w_org + w_qz + w_cal + w_Py + w_ilt + w_kln + w_mnt + w_bt + w_ms + w_urn, 4)
            else:
                sumMin = 0
            #
            w_H = round(chem_ilt[6][0]*w_ilt + chem_kln[6][0]*w_kln + chem_mnt[6][0]*w_mnt + chem_bt[6][0]*w_bt + chem_ms[6][0]*w_ms, 4)
            w_C = round(chem_org[6][0]*w_org + chem_cal[6][0]*w_cal, 4)
            w_O = round(chem_qz[6][0]*w_qz + chem_cal[6][1]*w_cal + chem_ilt[6][1]*w_ilt + chem_kln[6][1]*w_kln + chem_mnt[6][1]*w_mnt + chem_bt[6][1]*w_bt + chem_ms[6][1]*w_ms + chem_urn[6][0]*w_urn, 4)
            w_F = round(chem_bt[6][2]*w_bt + chem_ms[6][2]*w_ms, 4)
            w_Na = round(chem_mnt[6][2]*w_mnt, 4)
            w_Mg = round(chem_ilt[6][2]*w_ilt + chem_mnt[6][3]*w_mnt + chem_bt[6][3]*w_bt, 4)
            w_Al = round(chem_ilt[6][3]*w_ilt + chem_kln[6][2]*w_kln + chem_mnt[6][4]*w_mnt + chem_bt[6][4]*w_bt + chem_ms[6][3]*w_ms, 4)
            w_Si = round(chem_qz[6][1]*w_qz + chem_ilt[6][4]*w_ilt + chem_kln[6][3]*w_kln + chem_mnt[6][5]*w_mnt + chem_bt[6][5]*w_bt + chem_ms[6][4]*w_ms, 4)
            w_S = round(chem_py[6][0]*w_Py, 4)
            w_K = round(chem_ilt[6][5]*w_ilt + chem_bt[6][6]*w_bt + chem_ms[6][5]*w_ms, 4)
            w_Ca = round(chem_cal[6][2]*w_cal + chem_mnt[6][6]*w_mnt, 4)
            w_Fe = round(chem_py[6][1]*w_Py + chem_ilt[6][6]*w_ilt + chem_bt[6][7]*w_bt, 4)
            w_U = round(chem_urn[6][1]*w_urn, 6)
            sumConc = round(w_H + w_C + w_O + w_F + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe + w_U, 4)
            #
            GR = w_org*chem_org[5][0] + w_qz*chem_qz[5][0] + w_cal*chem_cal[5][0] + w_Py*chem_py[5][0] + w_ilt*chem_ilt[5][0] + w_kln*chem_kln[5][0] + w_mnt*chem_mnt[5][0] + w_bt*chem_bt[5][0] + w_ms*chem_ms[5][0] + w_urn*chem_urn[5][0]
            #
            if sumMin == 1 and sumConc == 1 and GR <= 300:
                cond = True
                #composition.extend((["Org", w_org, round(chem_org[1], 2)], ["Qz", w_qz, round(chem_qz[1], 2)], ["Cal", w_cal, round(chem_cal[1], 2)], ["Py", w_Py, round(chem_py[1], 2)], ["Ilt", w_ilt, round(chem_ilt[1], 2)], ["Kln", w_kln, round(chem_kln[1], 2)], ["Mnt", w_mnt, round(chem_mnt[1][0], 2), round(chem_mnt[1][1], 2)], ["Bt", w_bt, round(chem_bt[1][0], 2), round(chem_bt[1][1], 2), round(chem_bt[1][2], 2)], ["Ms", w_ms, round(chem_ms[1], 2)], ["Urn", w_urn, round(chem_urn[1], 2)]))
                composition.extend((["Org", "Qz", "Cal", "Py", "Ilt", "Kln", "Mnt", "Bt", "Ms", "Urn"]))
                concentrations = [w_H, w_C, w_O, w_F, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe, w_U]
                amounts = [w_org, w_qz, w_cal, w_Py, w_ilt, w_kln, w_mnt, w_bt, w_ms, w_urn]
            else:
                cond = False
        data.append(composition)
        #
        rhoSolid = (w_org*chem_org[2] + w_qz*chem_qz[2] + w_cal*chem_cal[2] + w_Py*chem_py[2] + w_ilt*chem_ilt[2] + w_kln*chem_kln[2] + w_mnt*chem_mnt[2] + w_bt*chem_bt[2] + w_ms*chem_ms[2] + w_urn*chem_urn[2]) / 1000
        X = [w_org, w_qz, w_cal, w_Py, w_ilt, w_kln, w_mnt, w_bt, w_ms, w_urn]
        K_list = [mineralogy[i][3][0] for i in range(len(mineralogy))]
        G_list = [mineralogy[i][3][1] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo/6
        G_solid = G_geo/6
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        E_solid = (9*K_solid*G_solid)/(3*K_solid+G_solid)
        nu_solid = (3*K_solid-2*G_solid)/(2*(3*K_solid+G_solid))
        #
        phi = randint(0, 5)/100
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_org*chem_org[5][0] + w_qz*chem_qz[5][0] + w_cal*chem_cal[5][0] + w_Py*chem_py[5][0] + w_ilt*chem_ilt[5][0] + w_kln*chem_kln[5][0] + w_mnt*chem_mnt[5][0] + w_bt*chem_bt[5][0] + w_ms*chem_ms[5][0] + w_urn*chem_urn[5][0]
        PE = w_org*chem_org[5][1] + w_qz*chem_qz[5][1] + w_cal*chem_cal[5][1] + w_Py*chem_py[5][1] + w_ilt*chem_ilt[5][1] + w_kln*chem_kln[5][1] + w_mnt*chem_mnt[5][1] + w_bt*chem_bt[5][1] + w_ms*chem_ms[5][1] + w_urn*chem_urn[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_org*chem_org[3][3] + w_qz*chem_qz[3][3] + w_cal*chem_cal[3][3] + w_Py*chem_py[3][3] + w_ilt*chem_ilt[3][3] + w_kln*chem_kln[3][3] + w_mnt*chem_mnt[3][3] + w_bt*chem_bt[3][3] + w_ms*chem_ms[3][3] + w_urn*chem_urn[3][3]
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
        return data
    #
class ore:
    #
    def __init__(self):
        pass
    #
    def create_ore_Fe(self):
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        chem_quartz = minerals.oxides.quartz("")
        chem_magnetite = minerals.oxides.magnetite("")
        chem_hematite = minerals.oxides.hematite("")
        chem_calcite = minerals.carbonates.calcite("")
        chem_alkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chem_plagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        #
        # [molar mass, density, bulk modulus, vP]
        chem_water = [18.0146, 997, 2.08, 1444]
        #
        ore_Fe = []
        #
        cond = False
        composition = []
        while cond == False:
            magicnumber = rd.randint(0, 2)
            if magicnumber == 0:    # Magnetite-rich
                x_ore = rd.randint(50,80)/100
                x_magnetite2 = rd.randint(75,100)/100
                x_hematite2 = 1-x_magnetite2
                x_magnetite = x_ore*x_magnetite2
                x_hematite = x_ore*x_hematite2
                x_gangue = 1-x_ore
                x_quartz2 = rd.randint(0,100)/100
                x_calcite2 = rd.randint(0,int((1-x_quartz2)*100))/100
                x_quartz = x_gangue*x_quartz2
                x_calcite = x_gangue*x_calcite2
                x_feldspars2 = 1-x_quartz2-x_calcite2
                x_alkalifeldspar2 = rd.randint(0,100)/100
                x_plagioclase2 = 1-x_alkalifeldspar2
                x_alkalifeldspar = x_gangue*x_feldspars2*x_alkalifeldspar2
                x_plagioclase = x_gangue*x_feldspars2*x_plagioclase2
            elif magicnumber == 1:    # Hematite-rich
                x_ore = rd.randint(50,80)/100
                x_hematite2 = rd.randint(75,100)/100
                x_magnetite2 = 1-x_hematite2
                x_hematite = x_ore*x_hematite2
                x_magnetite = x_ore*x_magnetite2
                x_gangue = 1-x_ore
                x_quartz2 = rd.randint(0,100)/100
                x_calcite2 = rd.randint(0,int((1-x_quartz2)*100))/100
                x_quartz = x_gangue*x_quartz2
                x_calcite = x_gangue*x_calcite2
                x_feldspars2 = 1-x_quartz2-x_calcite2
                x_alkalifeldspar2 = rd.randint(0,100)/100
                x_plagioclase2 = 1-x_alkalifeldspar2
                x_alkalifeldspar = x_gangue*x_feldspars2*x_alkalifeldspar2
                x_plagioclase = x_gangue*x_feldspars2*x_plagioclase2
            elif magicnumber == 2:    # No preference
                x_ore = rd.randint(50,80)/100
                x_magnetite2 = rd.randint(0,100)/100
                x_hematite2 = 1-x_magnetite2
                x_magnetite = x_ore*x_magnetite2
                x_hematite = x_ore*x_hematite2
                x_gangue = 1-x_ore
                x_quartz2 = rd.randint(0,100)/100
                x_calcite2 = rd.randint(0,int((1-x_quartz2)*100))/100
                x_quartz = x_gangue*x_quartz2
                x_calcite = x_gangue*x_calcite2
                x_feldspars2 = 1-x_quartz2-x_calcite2
                x_alkalifeldspar2 = rd.randint(0,100)/100
                x_plagioclase2 = 1-x_alkalifeldspar2
                x_alkalifeldspar = x_gangue*x_feldspars2*x_alkalifeldspar2
                x_plagioclase = x_gangue*x_feldspars2*x_plagioclase2
            sumMin = round(x_magnetite,2) + round(x_hematite,2) + round(x_quartz,2) + round(x_calcite,2) + round(x_alkalifeldspar,2) + round(x_plagioclase,2)
            if sumMin == 1:
                cond = True
                composition.extend((["Mag", round(x_magnetite,2), round(chem_magnetite[1],2)], ["Hem", round(x_hematite,2), round(chem_hematite[1],2)], ["Qz", round(x_quartz,2), round(chem_quartz[1],2)], ["Cal", round(x_calcite,2), round(chem_calcite[1],2)], ["Kfs", round(x_alkalifeldspar,2), round(chem_alkalifeldspar[1][0],2), round(chem_alkalifeldspar[1][1],2)], ["Pl", round(x_plagioclase,2), round(chem_plagioclase[1][0],2), round(chem_plagioclase[1][1],2)]))
            else:
                cond = False
        x_magnetite = composition[0][1]
        x_hematite = composition[1][1]
        x_quartz = composition[2][1]
        x_calcite = composition[3][1]
        x_alkalifeldspar = composition[4][1]
        x_plagioclase = composition[5][1]
        ore_Fe.append(composition)
        mineralogy = [chem_magnetite, chem_hematite, chem_quartz, chem_calcite, chem_alkalifeldspar, chem_plagioclase]
        #
        rhoSolid = (x_magnetite*chem_magnetite[2] + x_hematite*chem_hematite[2] + x_quartz*chem_quartz[2] + x_calcite*chem_calcite[2] + x_alkalifeldspar*chem_alkalifeldspar[2] + x_plagioclase*chem_plagioclase[2]) / 1000
        X = [x_magnetite, x_hematite, x_quartz, x_calcite, x_alkalifeldspar, x_plagioclase]
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
        phi = randint(0, 3) / 100
        rho = (1 - phi) * rhoSolid + phi * chem_water[1] / 1000
        vP = ((1 - phi) * vP_solid + phi * chem_water[3])/3
        vS = ((1 - phi) * vS_solid)/3
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chem_water[1] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = (x_magnetite*chem_magnetite[5][0] + x_hematite*chem_hematite[5][0] + x_quartz*chem_quartz[5][0] + x_calcite*chem_calcite[5][0] + x_alkalifeldspar*chem_alkalifeldspar[5][0] + x_plagioclase*chem_plagioclase[5][0])
        PE = x_magnetite*chem_magnetite[5][1] + x_hematite*chem_hematite[5][1] + x_quartz*chem_quartz[5][1] + x_calcite*chem_calcite[5][1] + x_alkalifeldspar*chem_alkalifeldspar[5][1] + x_plagioclase*chem_plagioclase[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = x_magnetite*chem_magnetite[3][3] + x_hematite*chem_hematite[3][3] + x_quartz*chem_quartz[3][3] + x_calcite*chem_calcite[3][3] + x_alkalifeldspar*chem_alkalifeldspar[3][3] + x_plagioclase*chem_plagioclase[3][3]
        #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
        #
        ore_Fe.append([round(rho, 3), round(rhoSolid, 3), round(chem_water[1] / 1000, 6)])
        ore_Fe.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        ore_Fe.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chem_water[3], 2)])
        ore_Fe.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        ore_Fe.append("water")
        ore_Fe.append([GR, PE])
        #
        #  ore_Fe = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, [GR,PE]]
        #
        return ore_Fe

class Sandstone():
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
    def create_sandstone(self):
        # Major elements
        quartz = oxides.Quartz(impurity="random").create_quartz()
