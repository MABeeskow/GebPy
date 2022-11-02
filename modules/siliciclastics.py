#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		siliciclastics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		22.08.2022

#-----------------------------------------------

## MODULES
import datetime
import sys

import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals, oxides, fluids
from modules.geophysics import Elasticity as elast
from modules import oxides, carbonates, silicates
from modules.oxides import Oxides
from modules.carbonates import Carbonates
from modules.silicates import Phyllosilicates
from modules.silicates import Tectosilicates
from modules.sulfides import Sulfides
from modules.organics import Organics
from modules.fluids import Water

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
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
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
    
    def create_simple_sandstone(self, w_Fe=None, amounts=None, porosity=None, pure=False, dict_output=False):
        #
        results = {}
        results["rock"] = "Sandstone"
        #
        self.w_Fe = w_Fe
        self.amounts = amounts
        #
        #
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
        data_muscovite = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
        #
        minerals_list = [self.data_quartz, self.data_hematite, data_alkalifeldspar, data_plagioclase, data_chlorite,
                         data_muscovite, self.data_calcite]
        #
        elements_list = []
        for mineral in minerals_list:
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements_list:
                    elements_list.append(element)
        elements_list.sort()
        #
        w_FeChl = data_chlorite["chemistry"]["Fe"]
        w_FeHem = self.data_hematite["chemistry"]["Fe"]
        #
        air = fluids.Gas.air("")
        water = fluids.Water.water("")
        oil = fluids.Hydrocarbons.oil("")
        gas = fluids.Hydrocarbons.natural_gas("")
        #
        data = []
        #
        cond = False
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
                    w_Fsp = round(rd.randint(25, int((1-w_ore)*50))/100, 4)
                    w_Afs2 = rd.randint(75, 100)/100
                    w_Pl2 = 1 - w_Afs2
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Qz = round(rd.randint(50, int((1-w_ore-w_Fsp)*100))/100, 4)
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
                    w_Fsp = round(rd.randint(25, int((1-w_ore)*50))/100, 4)
                    w_Pl2 = rd.randint(75, 100)/100
                    w_Afs2 = 1 - w_Pl2
                    w_Pl = round(w_Fsp*w_Pl2, 4)
                    w_Afs = round(w_Fsp*w_Afs2, 4)
                    w_Qz = round(rd.randint(50, int((1-w_ore-w_Fsp)*100))/100, 4)
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
                    w_Qz = round(rd.randint(50, int((1-w_rf)*100))/100, 4)
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
                    w_Qz = round(rd.randint(50, int((1-w_rf)*100))/100, 4)
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
                    w_ore = round(rd.uniform(0.1, 0.2), 4)
                    w_Hem = w_ore
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
            elif isinstance(self.amounts, dict):
                w_Qz = round(abs(np.random.normal(self.amounts["Qz"], 0.025)), 4)
                w_Afs = round(abs(np.random.normal(self.amounts["Kfs"], 0.025)), 4)
                w_Pl = round(abs(np.random.normal(self.amounts["Pl"], 0.025)), 4)
                w_Cal = round(abs(np.random.normal(self.amounts["Cal"], 0.025)), 4)
                w_Chl = round(abs(np.random.normal(self.amounts["Chl"], 0.025)), 4)
                w_Ms = round(abs(np.random.normal(self.amounts["Ms"], 0.025)), 4)
                w_Hem = round(1-w_Qz-w_Afs-w_Pl-w_Cal-w_Chl-w_Ms, 4)
            #
            if w_Qz >= 0.0 and w_Afs >= 0.0 and w_Pl >= 0.0 and w_Cal >= 0.0 and w_Chl >= 0.0 and w_Ms  >= 0.0 and w_Hem >= 0.0:
                sumMin = round(w_Qz + w_Afs + w_Pl + w_Cal + w_Chl + w_Ms + w_Hem, 4)
            else:
                sumMin = 0
            #
            mineral_amounts = {}
            mineral_amounts["Qz"] = w_Qz
            mineral_amounts["Hem"] = w_Hem
            mineral_amounts["Kfs"] = w_Afs
            mineral_amounts["Pl"] = w_Pl
            mineral_amounts["Chl"] = w_Chl
            mineral_amounts["Ms"] = w_Ms
            mineral_amounts["Cal"] = w_Cal
            #
            element_amounts = {}
            w_O = 1
            w_sum = 0
            for element in elements_list:
                if element != "O":
                    element_amounts[element] = 0
                    for mineral in minerals_list:
                        if element in mineral["chemistry"]:
                            element_amounts[element] += round(mineral["chemistry"][element]*mineral_amounts[mineral["mineral"]], 4)
                    element_amounts[element] = round(element_amounts[element], 4)
                    w_sum += element_amounts[element]
                    w_O -= element_amounts[element]
            element_amounts["O"] = round(w_O, 4)
            w_sum += element_amounts["O"]
            if sumMin == 1 and w_sum == 1:
                cond = True
            else:
                cond = False
        #
        results["mineralogy"] = mineral_amounts
        results["chemistry"] = element_amounts
        #
        rhoSolid = 0
        K_list = []
        G_list = []
        for mineral in minerals_list:
            rhoSolid += mineral["rho"]*mineral_amounts[mineral["mineral"]]/1000
            K_list.append(round(mineral["K"]*mineral_amounts[mineral["mineral"]], 3))
            G_list.append(round(mineral["G"]*mineral_amounts[mineral["mineral"]], 3))
        X = [w_Qz, w_Hem, w_Afs, w_Pl, w_Chl, w_Ms, w_Cal]

        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        #K_solid = K_geo/2
        #G_solid = G_geo/2
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
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
        phi = round(phi, 4)
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        if self.fluid == "water" or self.w_Fe != None:
            rho = (1 - phi)*rhoSolid + phi*water[2]/1000
            vP = (1 - phi)*vP_solid + phi*water[4][0]
            vS = (1 - phi)*vS_solid
            #
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            #
            GR = 0
            PE = 0
            poisson_mineralogical = 0
            for mineral in minerals_list:
                GR += mineral["GR"]*mineral_amounts[mineral["mineral"]]
                PE += mineral["PE"]*mineral_amounts[mineral["mineral"]]
                poisson_mineralogical += mineral["nu"]*mineral_amounts[mineral["mineral"]]
            #
        elif self.fluid == "oil":
            rho = (1 - phi) * rhoSolid + phi * oil[2] / 1000
            vP = (1-phi)*vP_solid + phi*oil[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - oil[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            #
            GR = 0
            PE = 0
            poisson_mineralogical = 0
            for mineral in minerals_list:
                GR += mineral["GR"]*mineral_amounts[mineral["mineral"]]
                PE += mineral["PE"]*mineral_amounts[mineral["mineral"]]
                poisson_mineralogical += mineral["nu"]*mineral_amounts[mineral["mineral"]]
            #
        elif self.fluid == "gas":
            rho = (1 - phi) * rhoSolid + phi * gas[2] / 1000
            vP = (1-phi)*vP_solid + phi*gas[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - gas[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            #
            GR = 0
            PE = 0
            poisson_mineralogical = 0
            for mineral in minerals_list:
                GR += mineral["GR"]*mineral_amounts[mineral["mineral"]]
                PE += mineral["PE"]*mineral_amounts[mineral["mineral"]]
                poisson_mineralogical += mineral["nu"]*mineral_amounts[mineral["mineral"]]
            #
        elif self.fluid == "air":
            rho = (1 - phi) * rhoSolid + phi * air[2] / 1000
            vP = (1-phi)*vP_solid + phi*air[4][0]
            vS = (1 - phi) * vS_solid
            G_bulk = vS**2 * rho
            K_bulk = vP**2 * rho - 4/3*G_bulk
            E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
            phiD = (rhoSolid - rho) / (rhoSolid - air[2] / 1000)
            phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
            #
            GR = 0
            PE = 0
            poisson_mineralogical = 0
            for mineral in minerals_list:
                GR += mineral["GR"]*mineral_amounts[mineral["mineral"]]
                PE += mineral["PE"]*mineral_amounts[mineral["mineral"]]
                poisson_mineralogical += mineral["nu"]*mineral_amounts[mineral["mineral"]]
        #
        if dict_output == False:
            data.append([round(rho, 3), round(rhoSolid, 3), round(gas[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(gas[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(self.fluid)
            data.append([round(GR, 3), round(PE, 3)])
            data.append(amounts)
            #
            return data
        else:
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            return results
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
            #w_O = round(w_Qz*chem_quartz[6][0] + w_Afs*chem_alkalifeldspar[6][0] + w_Pl*chem_plagioclase[6][0], 4)
            w_Na = round(w_Afs*chem_alkalifeldspar[6][1] + w_Pl*chem_plagioclase[6][1], 4)
            w_Al = round(w_Afs*chem_alkalifeldspar[6][2] + w_Pl*chem_plagioclase[6][2], 4)
            w_Si = round(w_Qz*chem_quartz[6][1] + w_Afs*chem_alkalifeldspar[6][3] + w_Pl*chem_plagioclase[6][3], 4)
            w_K = round(w_Afs*chem_alkalifeldspar[6][4], 4)
            w_Ca = round(w_Pl*chem_plagioclase[6][4], 4)
            w_O = round(1 - w_Na - w_Al - w_Si - w_K - w_Ca, 4)
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
    def __init__(self, fluid=None):
        self.fluid = fluid
        #
        self.data_qz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_cal = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_kln = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        self.data_py = Sulfides(impurity="pure", data_type=True).create_pyrite()
        self.data_urn = Oxides(impurity="pure", data_type=True).create_uraninite()
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
    def create_simple_shale(self, w_C=None, w_F=None, w_Na=None, w_Mg=None, w_S=None, w_K=None, w_Ca=None, w_Fe=None,
                            amounts=None, porosity=None, dict_output=False):
        #
        results = {}
        results["rock"] = "Shale"
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
        data_org = Organics(data_type=True).create_organics_matter()
        data_ms = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
        data_bt = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        data_ilt = Phyllosilicates(impurity="pure", data_type=True).create_illite()
        data_mnt = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
        #
        minerals_list = [data_org, self.data_qz, data_ms, data_bt, self.data_cal, data_ilt, self.data_kln,
                         data_mnt, self.data_py, self.data_urn]
        #
        elements_list = []
        for mineral in minerals_list:
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements_list:
                    elements_list.append(element)
        elements_list.sort()
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        while cond == False:
            if self.w_C == None and self.w_F == None and self.w_Na == None and self.w_Mg == None and self.w_S == None and self.w_K == None and self.w_Ca == None and self.w_Fe == None and self.amounts == None:
                magicnumber = rd.randint(0, 2)
                if magicnumber == 0:    # Clay-rich
                    w_ore = round(rd.randint(0, 5)/100, 4)
                    w_Py = w_ore
                    w_clay = round(rd.randint(45, int((0.70-w_ore)*100))/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_mnt2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_kln2 = 1-w_ilt2-w_mnt2
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
                    w_urn = round(rd.uniform(0.000001, 0.00001), 6)
                    w_org = round(1 - w_Py - w_ilt - w_kln - w_mnt - w_qz - w_ms - w_bt - w_cal - w_urn, 6)
                    #
                elif magicnumber == 1:    # Qz-rich
                    w_ore = round(rd.randint(0, 5)/100, 4)
                    w_Py = w_ore
                    w_qz = round(rd.randint(25, 50)/100, 4)
                    w_clay = round(rd.randint(40, int((1-w_ore-w_qz)*100))/100, 4)
                    w_ilt2 = rd.randint(50, 100)/100
                    w_mnt2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_kln2 = 1-w_ilt2-w_mnt2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_mica = round(rd.randint(0, int((1-w_ore-w_clay-w_qz)*100))/100, 4)
                    w_ms2 = rd.randint(0, 100)/100
                    w_bt2 = 1 - w_ms2
                    w_ms = round(w_mica*w_ms2, 4)
                    w_bt = round(w_mica*w_bt2, 4)
                    w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                    w_urn = round(rd.uniform(0.000001, 0.00001), 6)
                    w_org = round(1 - w_Py - w_ilt - w_kln - w_mnt - w_qz - w_ms - w_bt - w_cal - w_urn, 6)
                    #
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
                    w_mnt2 = rd.randint(0, int((1-w_ilt2)*100))/100
                    w_kln2 = 1-w_ilt2-w_mnt2
                    w_ilt = round(w_clay*w_ilt2, 4)
                    w_kln = round(w_clay*w_kln2, 4)
                    w_mnt = round(w_clay*w_mnt2, 4)
                    w_qz = round(rd.randint(0, int((1-w_ore-w_clay-w_mica)*100))/100, 4)
                    w_cal = round(rd.randint(0, int((1-w_ore-w_clay-w_qz-w_mica)*100))/100, 4)
                    w_urn = round(rd.uniform(0.000001, 0.00001), 6)
                    w_org = round(1 - w_Py - w_ilt - w_kln - w_mnt - w_qz - w_ms - w_bt - w_cal - w_urn, 6)
            #
            elif type(self.amounts) is list:
                w_urn = round(abs(np.random.normal(self.amounts[9], 5e-07)), 6)
                w_org = round(abs(np.random.normal(self.amounts[0], 0.025)), 6)
                w_qz = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_Py = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts[4], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts[5], 0.025)), 4)
                w_mnt = round(abs(np.random.normal(self.amounts[6], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts[7], 0.025)), 4)
                w_ms = round(1 - w_urn - w_org - w_qz - w_cal - w_Py - w_ilt - w_kln - w_mnt - w_bt, 4)
            elif isinstance(self.amounts, dict):
                w_urn = round(abs(np.random.normal(self.amounts["Urn"], 5e-07)), 6)
                w_org = round(abs(np.random.normal(self.amounts["Org"], 0.025)), 6)
                w_qz = round(abs(np.random.normal(self.amounts["Qz"], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts["Cal"], 0.025)), 4)
                w_Py = round(abs(np.random.normal(self.amounts["Py"], 0.025)), 4)
                w_ilt = round(abs(np.random.normal(self.amounts["Ilt"], 0.025)), 4)
                w_kln = round(abs(np.random.normal(self.amounts["Kln"], 0.025)), 4)
                w_mnt = round(abs(np.random.normal(self.amounts["Mnt"], 0.025)), 4)
                w_bt = round(abs(np.random.normal(self.amounts["Bt"], 0.025)), 4)
                w_ms = round(1 - w_urn - w_org - w_qz - w_cal - w_Py - w_ilt - w_kln - w_mnt - w_bt, 4)
            #
            if w_org >= 0.0 and w_qz >= 0.0 and w_cal >= 0.0 and w_Py >= 0.0 and w_ilt >= 0.0 and w_kln >= 0.0 \
                    and w_mnt >= 0.0 and w_bt >= 0.0 and w_ms >= 0.0 and w_urn >= 0.0:
                sumMin = round(w_org + w_qz + w_cal + w_Py + w_ilt + w_kln + w_mnt + w_bt + w_ms + w_urn, 6)
            else:
                sumMin = 0
            #
            mineral_amounts = {}
            mineral_amounts["Org"] = w_org
            mineral_amounts["Qz"] = w_qz
            mineral_amounts["Ms"] = w_ms
            mineral_amounts["Bt"] = w_bt
            mineral_amounts["Cal"] = w_cal
            mineral_amounts["Ilt"] = w_ilt
            mineral_amounts["Kln"] = w_kln
            mineral_amounts["Mnt"] = w_mnt
            mineral_amounts["Py"] = w_Py
            mineral_amounts["Urn"] = w_urn
            #
            element_amounts = {}
            w_O = 1
            w_sum = 0
            for element in elements_list:
                if element != "O":
                    element_amounts[element] = 0
                    for mineral in minerals_list:
                        if element in mineral["chemistry"]:
                            element_amounts[element] += round(mineral["chemistry"][element]*abs(mineral_amounts[mineral["mineral"]]), 6)
                    if element != "U":
                        element_amounts[element] = round(element_amounts[element], 4)
                    else:
                        element_amounts[element] = round(element_amounts[element], 6)
                    w_sum += element_amounts[element]
                    w_O -= element_amounts[element]
                else:
                    pass
            element_amounts["O"] = round(w_O, 6)
            w_sum += element_amounts["O"]
            if sumMin == 1 and w_sum == 1:
                cond = True
            else:
                cond = False
        #
        results["mineralogy"] = mineral_amounts
        results["chemistry"] = element_amounts
        #
        rhoSolid = 0
        K_list = []
        G_list = []
        for mineral in minerals_list:
            rhoSolid += mineral["rho"]*mineral_amounts[mineral["mineral"]]/1000
            K_list.append(round(mineral["K"]*mineral_amounts[mineral["mineral"]], 3))
            G_list.append(round(mineral["G"]*mineral_amounts[mineral["mineral"]], 3))
        X = [w_org, w_qz, w_ms, w_bt, w_cal, w_ilt, w_kln, w_mnt, w_Py, w_urn]
        #
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo/1.5
        G_solid = G_geo/1.5
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        #
        if porosity == None:
            phi = rd.uniform(0.01, 0.05)
        else:
            phi = porosity
        phi = round(phi, 4)
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        rho = (1 - phi)*rhoSolid + phi*water[2]/1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        #
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        #
        GR = 0
        PE = 0
        poisson_mineralogical = 0
        for mineral in minerals_list:
            GR += mineral["GR"]*mineral_amounts[mineral["mineral"]]
            PE += mineral["PE"]*mineral_amounts[mineral["mineral"]]
            poisson_mineralogical += mineral["nu"]*mineral_amounts[mineral["mineral"]]
        #
        if dict_output == False:
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 3)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append(self.fluid)
            data.append([round(GR, 3), round(PE, 3)])
            data.append(amounts)
            #
            return data
        else:
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            return results
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

# class Sandstone():
#     """ Class that generates geophysical and geochemical data of magnetite"""
#     #
#     def __init__(self, traces_list=[], impurity="pure"):
#         self.traces_list = traces_list
#         self.impurity = impurity
#         if len(self.traces_list) > 0:
#             self.impurity = "impure"
#         if self.impurity == "random":
#             minors = ["Mg", "Zn", "Mn", "Ni", "Cr", "Ti", "V", "Al"]
#             n = rd.randint(1, len(minors))
#             while len(self.traces_list) < n:
#                 selection = rd.choice(minors)
#                 if selection not in self.traces_list:
#                     self.traces_list.append(selection)
#                 else:
#                     continue
#     #
#     def create_sandstone(self):
#         # Major elements
#         quartz = oxides.Quartz(impurity="random").create_quartz()

#####################
## SANDSTONE ROCKS ##
#####################
class Sandstone:
    #
    def __init__(self, fluid="water", actualThickness=100):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()
        self.data_calcite = Carbonates(impurity="pure", data_type=True).create_calcite()
        self.data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()
        self.data_water = Water.water("")
    #
    def create_sandstone(self, number, porosity=None):
        #
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()   # variable
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()         # variable
        data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()              # variable
        #
        assemblage = [self.data_quartz, data_alkalifeldspar, data_plagioclase,  self.data_calcite, self.data_hematite,
                      data_chlorite]
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
                if mineral == "Qz":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(1 - w_total, 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.85 <= value <= 1.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Cal":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Chl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.15), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
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
                phi_helper = round(rd.uniform(0.0, 0.05), 4)
            else:
                try:
                    phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
                except:
                    phi_helper = round(porosity, 4)
            #
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_chlorite = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Qz"][n] * self.data_quartz["chemistry"][element], 4)
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
                if element in data_chlorite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Chl"][n] * data_chlorite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_hematite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Hem"][n] * self.data_hematite["chemistry"][element], 4)
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
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_alkalifeldspar["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["PE"], 3)
                    elif mineral == "Chl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_chlorite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_chlorite["PE"], 3)
                    elif mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["PE"], 3)
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
        results["rock"] = "Sandstone"
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
    def create_conglomerate_alt(self, number=1, composition=None, porosity=None):
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        mineralogy = {"Qz": self.data_quartz, "Kfs": data_alkalifeldspar, "Pl": data_plagioclase, "Bt": data_biotite,
                      "Cal": self.data_calcite, "Hem": self.data_hematite}
        #
        condition = False
        #
        while condition == False:
            elements_list = []
            phi_minerals = {}
            w_minerals = {}
            w_elements = {}
            #
            if composition != None:
                phi_qz = composition["Qz"]
                phi_kfs = composition["Kfs"]
                phi_bt = composition["Bt"]
                phi_cal = composition["Cal"]
                phi_hem = composition["Hem"]
                #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
                phi_minerals["Cal"] = phi_cal
                phi_minerals["Hem"] = phi_hem
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_qz = round(rd.uniform(0.5, 0.8), 4)
                    phi_kfs = round(rd.uniform(0.0, (1.0 - phi_qz)), 4)
                    phi_pl = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs)), 4)
                    phi_bt = round(rd.uniform(0.0, (1.0 - phi_qz - phi_kfs - phi_pl)), 4)
                    phi_cal = round(rd.uniform(0.1, (1.0 - phi_qz - phi_kfs - phi_pl - phi_bt)), 4)
                    phi_hem = round(1 - phi_qz - phi_kfs - phi_pl - phi_bt - phi_cal, 4)
                    phi_total = phi_qz + phi_kfs + phi_pl + phi_bt + phi_cal + phi_hem
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.5 <= phi_qz <= 0.8 and 0.0 <= phi_kfs <= 0.2 and 0.0 <= phi_pl <= 0.2 \
                                and 0.0 <= phi_bt <= 0.05 and 0.1 <= phi_cal <= 0.2 and 0.0 <= phi_hem <= 0.05:
                            condition_2 = True
                    #
                phi_minerals["Qz"] = phi_qz
                phi_minerals["Kfs"] = phi_kfs
                phi_minerals["Pl"] = phi_pl
                phi_minerals["Bt"] = phi_bt
                phi_minerals["Cal"] = phi_cal
                phi_minerals["Hem"] = phi_hem
            #
            rho_s = 0
            for key, value in phi_minerals.items():
                rho_s += phi_minerals[key] * mineralogy[key]["rho"]
                for element, value in mineralogy[key]["chemistry"].items():
                    if element not in elements_list:
                        elements_list.append(element)
                        w_elements[element] = 0.0
            rho_s = round(rho_s, 3)
            for key, value in phi_minerals.items():
                w_minerals[key] = round((phi_minerals[key] * mineralogy[key]["rho"]) / rho_s, 4)
            #
            if porosity == None:
                porosity = round(rd.uniform(0.0, 0.1), 4)
            rho = round((1 - porosity) * rho_s + porosity * self.data_water[2] / 1000, 3)
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
        youngs_mod = round((9 * bulk_mod * shear_mod) / (3 * bulk_mod + shear_mod), 3)
        poisson_rat = round((3 * bulk_mod - 2 * shear_mod) / (6 * bulk_mod + 2 * shear_mod), 4)
        vP = round(((bulk_mod * 10 ** 9 + 4 / 3 * shear_mod * 10 ** 9) / (rho)) ** 0.5, 3)
        vS = round(((shear_mod * 10 ** 9) / (rho)) ** 0.5, 3)
        vPvS = round(vP / vS, 3)
        #
        results = {}
        results["rock"] = "Conglomerate"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = porosity
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
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in w_minerals.items():
                single_amounts_mineralogy[mineral] = value
            for element, value in w_elements.items():
                single_amounts_chemistry[element] = value
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = porosity
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
    def create_conglomerate(self, number, porosity=None):
        #
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
        #
        assemblage = [self.data_quartz, data_alkalifeldspar, data_plagioclase, data_biotite, self.data_calcite,
                      self.data_hematite]
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
                if mineral == "Qz":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(1 - w_total, 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.5 <= value <= 0.8:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kfs":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Pl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Cal":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.1, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.1 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Bt":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.05), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
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
                phi_helper = round(rd.uniform(0.0, 0.2), 4)
            else:
                phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
            #
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            data_biotite = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Qz"][n] * self.data_quartz["chemistry"][element], 4)
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
                if element in data_biotite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Bt"][n] * data_biotite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_hematite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Hem"][n] * self.data_hematite["chemistry"][element], 4)
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
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_alkalifeldspar["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["PE"], 3)
                    elif mineral == "Bt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_biotite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_biotite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_biotite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_biotite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_biotite["PE"], 3)
                    elif mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_hematite["PE"], 3)
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
        results["rock"] = "Conglomerate"
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
    def create_mudstone(self, number, porosity=None, dominance="Mnt", dominant_group="Clays"):
        #
        data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
        data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
        data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
        data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
        data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
        data_organics = Organics(data_type=True).create_organics_matter()
        #
        assemblage = [data_montmorillonite, data_illite, data_kaolinite, self.data_quartz, data_alkalifeldspar,
                      data_plagioclase, self.data_calcite]
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
        mineral_list_clays = ["Ilt", "Kln", "Mnt"]
        mineral_list_silic = ["Kfs", "Pl", "Qz"]
        mineral_list_carb = ["Cal"]
        elements.sort()
        #
        n = 0
        amounts_helper = []
        while n < number:
            if dominant_group == "Clays":
                w_clay = round(rd.uniform(0.4, 0.7), 4)
                w_silic = round(rd.uniform(0.0, (1 - w_clay)), 4)
                w_carb = round(1 - w_clay - w_silic, 4)
            elif dominant_group == "Siliciclastics":
                w_silic = round(rd.uniform(0.4, 0.7), 4)
                w_clay = round(rd.uniform(0.0, (1 - w_silic)), 4)
                w_carb = round(1 - w_silic - w_clay, 4)
            elif dominant_group == "Carbonates":
                w_carb = round(rd.uniform(0.4, 0.7), 4)
                w_clay = round(rd.uniform(0.0, (1 - w_carb)), 4)
                w_silic = round(1 - w_carb - w_clay, 4)
            mineral_fractions = {"Clays": w_clay, "Siliciclastics": w_silic, "Carbonates": w_carb}
            w_total = 0
            w_total_carb = 0
            w_total_silic = 0
            w_total_clay = 0
            n_minerals = 0
            n_carb = 0
            n_silic = 0
            n_clay = 0
            #
            for group, fraction in mineral_fractions.items():
                if group == "Clays":
                    for mineral in mineral_list_clays:
                        if dominance == "Mnt" and n_minerals < 3:
                            if mineral == "Mnt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.8, 1.0), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.6 <= value <= 1.0:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                            elif mineral == "Ilt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.0, 0.05), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.0 <= value <= 0.05:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                            elif mineral == "Kln":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.0, 0.05), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.0 <= value <= 0.05:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                        elif dominance == "Ilt" and n_minerals < 3:
                            if mineral == "Ilt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.8, 1.0), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.6 <= value <= 1.0:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                            elif mineral == "Mnt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.0, 0.05), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.0 <= value <= 0.05:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                            elif mineral == "Kln":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.0, 0.05), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.0 <= value <= 0.05:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                        elif dominance == "Kln" and n_minerals < 3:
                            if mineral == "Kln":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.8, 1.0), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.6 <= value <= 1.0:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                            elif mineral == "Mnt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.0, 0.05), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.0 <= value <= 0.05:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                            elif mineral == "Ilt":
                                if n_minerals < len(mineral_list) - 1:
                                    value = round(fraction * rd.uniform(0.0, 0.05), 4)
                                else:
                                    value = round(1 - w_total, 4)
                                if value >= 0.0 and 0.0 <= value <= 0.05:
                                    amounts_helper.append(value)
                                    w_total += value
                                    w_total_clay += value
                                    n_minerals += 1
                                    n_clay += 1
                #
                elif group == "Siliciclastics":
                    for mineral in mineral_list_silic:
                        if mineral == "Qz":
                            if n_minerals < len(mineral_list) - 1:
                                value = round(fraction * rd.uniform(0.0, 1.0), 4)
                            elif n_silic == 2:
                                value = round(w_silic - w_total_silic, 4)
                            else:
                                value = round(1 - w_total, 4)
                            if value >= 0.0 and 0.0 <= value <= fraction:
                                amounts_helper.append(value)
                                w_total += value
                                w_total_silic += value
                                n_minerals += 1
                                n_silic += 1
                        elif mineral == "Kfs":
                            if n_minerals < len(mineral_list) - 1:
                                value = round(fraction * rd.uniform(0.0, 1.0), 4)
                            else:
                                value = round(1 - w_total, 4)
                            if value >= 0.0 and 0.0 <= value <= fraction:
                                amounts_helper.append(value)
                                w_total += value
                                w_total_silic += value
                                n_minerals += 1
                                n_silic += 1
                        elif mineral == "Pl":
                            if n_minerals < len(mineral_list) - 1:
                                value = round(fraction * rd.uniform(0.0, 1.0), 4)
                            else:
                                value = round(1 - w_total, 4)
                            if value >= 0.0 and 0.0 <= value <= fraction:
                                amounts_helper.append(value)
                                w_total += value
                                w_total_silic += value
                                n_minerals += 1
                                n_silic += 1
                elif group == "Carbonates":
                    for mineral in mineral_list_carb:
                        if mineral == "Cal":
                            if n_minerals < len(mineral_list) - 1:
                                value = round(fraction * rd.uniform(0.0, 1.0), 4)
                            elif n_silic == 2:
                                value = round(w_carb - w_total_carb, 4)
                            else:
                                value = round(1 - w_total, 4)
                            if value >= 0.0 and 0.0 <= value <= fraction:
                                amounts_helper.append(value)
                                w_total += value
                                w_total_carb += value
                                n_minerals += 1
                                n_carb += 1
            #
            mineral_list_new = []
            mineral_list_new.extend(mineral_list_clays)
            mineral_list_new.extend(mineral_list_silic)
            mineral_list_new.extend(mineral_list_carb)
            if np.sum(amounts_helper) == 1.0 and n_minerals == len(mineral_list_new):
                for index, mineral in enumerate(mineral_list_new):
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
            data_montmorillonite = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
            data_illite = Phyllosilicates(impurity="pure", data_type=True).create_illite()
            data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
            data_alkalifeldspar = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            data_plagioclase = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            #
            old_index = elements.index("O")
            elements += [elements.pop(old_index)]
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in data_montmorillonite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Mnt"][n] * data_montmorillonite["chemistry"][element], 4)
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
                if element in data_kaolinite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Kln"][n] * data_kaolinite["chemistry"][element], 4)
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
                if element in self.data_calcite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Cal"][n] * self.data_calcite["chemistry"][element], 4)
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
                    if mineral == "Mnt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * data_montmorillonite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_montmorillonite["PE"], 3)
                    elif mineral == "Ilt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_illite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_illite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_illite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_illite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_illite["PE"], 3)
                    elif mineral == "Kln":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * data_kaolinite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["PE"], 3)
                    elif mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_quartz["PE"], 3)
                    elif mineral == "Kfs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * data_alkalifeldspar["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_alkalifeldspar["PE"], 3)
                    elif mineral == "Pl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * data_plagioclase["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_plagioclase["PE"], 3)
                #
                anisotropic_factor = round(rd.uniform(1.0, 2.0), 2)
                bulkmod_helper = bulkmod_helper/anisotropic_factor
                shearmod_helper = shearmod_helper / anisotropic_factor
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
        results["rock"] = "Mudstone"
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
#
## TEST
# print(Sandstone(fluid="water", actualThickness=0).create_sandstone(number=100))
# print(Sandstone(fluid="water", actualThickness=0).create_conglomerate(number=10))