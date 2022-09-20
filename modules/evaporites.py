#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		evaporites.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		20.09.2022

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
from modules import minerals, geochemistry
import random as rd
from modules import fluids
from modules.geophysics import Elasticity as elast
from modules.sulfates import Sulfates
from modules.halides import Halides
from modules.carbonates import Carbonates
from modules.fluids import Water

class evaporites:
    #
    def __init__(self,):
        pass
    #
    def createRocksalt(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemHalite = minerals.halides.halite("")
        chemSylvite = minerals.halides.sylvite("")
        chemAnhydrite = minerals.sulfates.anhydrite("")
        chemGypsum = minerals.sulfates.gypsum("")
        chemIllite = minerals.phyllosilicates.illite("")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = fluids.Water.water("")
        #
        rocksalt = []
        #
        cond = False
        composition = []
        while cond == False:
            xHalides = round(randint(90, 100)/100, 2)
            xHalite2 = round(randint(90, 100)/100, 2)
            xSylvite2 = 1 - xHalite2
            xHalite = round(xHalides * xHalite2, 2)
            xSylvite = round(xHalides * xSylvite2, 2)
            xSulfates = round(randint(0, 10) / 100, 2)
            xAnhydrite2 = round(randint(0, 100) / 100, 2)
            xGypsum2 = 1-xAnhydrite2
            xAnhydrite = round(xSulfates * xAnhydrite2, 2)
            xGypsum = round(xSulfates * xGypsum2, 2)
            xIllite = round(randint(0, 5)/100, 2)
            sumMin = round(xHalite + xAnhydrite + xGypsum + xSylvite + xIllite, 2)
            #print("sumMin:", sumMin)
            if sumMin == 1:
                cond = True
                composition.extend([["Hl", round(xHalite,2), round(chemHalite[1],2)], ["Anh", round(xAnhydrite,2), round(chemAnhydrite[1],2)], ["Gp", round(xGypsum,2), round(chemGypsum[1],2)], ["Syl", round(xSylvite, 2)], ["Ilt", round(xIllite,2), round(chemIllite[1],2)]])
            else:
                cond = False
        xHalite = composition[0][1]
        xAnhydrite = composition[1][1]
        xGypsum = composition[2][1]
        xSylvite = composition[3][1]
        xIllite = composition[4][1]
        rocksalt.append(composition)
        mineralogy = [chemHalite, chemAnhydrite, chemGypsum, chemSylvite, chemIllite]
        #
        rhoSolid = (xHalite*chemHalite[2] + xAnhydrite * chemAnhydrite[2] + xGypsum * chemGypsum[2] + xSylvite * chemSylvite[2] + xIllite * chemIllite[2]) / 1000
        X = [xHalite, xAnhydrite, xGypsum, xSylvite, xIllite]
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
        #phi = randint(0, 0)/100
        phi = 0.0
        rho = (1 - phi) * rhoSolid + phi * chemWater[2] / 1000
        vP = (1-phi)*vP_solid + phi*chemWater[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = xHalite*chemHalite[5][0] + xAnhydrite*chemAnhydrite[5][0] + xGypsum*chemGypsum[5][0] + xSylvite*chemSylvite[5][0] + xIllite*chemIllite[5][0]
        PE = xHalite*chemHalite[5][1] + xAnhydrite*chemAnhydrite[5][1] + xGypsum*chemGypsum[5][1] + xSylvite*chemSylvite[5][1] + xIllite*chemIllite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = xHalite*chemHalite[3][3] + xAnhydrite*chemAnhydrite[3][3] + xGypsum*chemGypsum[3][3] + xSylvite*chemSylvite[3][3] + xIllite*chemIllite[3][3]
        #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
        #
        rocksalt.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        rocksalt.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        rocksalt.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[4][0], 2)])
        rocksalt.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        rocksalt.append("water")
        rocksalt.append([GR, PE])
        #
        #  rocksalt = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return rocksalt
        #
    #
    def createAnhydrite(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemAnhydrite = minerals.sulfates.anhydrite("")
        chemCalcite = minerals.carbonates.calcite("")
        chemHalite = minerals.halides.halite("")
        chemGalena = minerals.sulfides.galena("")
        chemChalcopyrite = minerals.sulfides.chalcopyrite("")
        chemMolybdenite = minerals.sulfides.molybdenite("")
        chemPyrite = minerals.sulfides.pyrite("")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = fluids.Water.water("")
        #
        anhydrite = []
        #
        cond = False
        composition = []
        magicnumber = rd.randint(0,100)
        while cond == False:
            xAnhydrite = round(randint(90, 100)/100, 2)
            xCalcite = round(randint(0, 10)/100, 2)
            xHalite = round(randint(0, 10)/100, 2)
            xSulfides = round(randint(0, 5) / 100, 2)
            if magicnumber in range(0,5):
                xGalena2 = round(randint(0, 100) / 100, 2)
                xChalcopyrite2 = 0.00
                xMolybdenite2 = 1-xGalena2
                xPyrite2 = 0.00
            elif magicnumber in range(5,20):
                xGalena2 = 0.00
                xChalcopyrite2 = round(randint(0, 100) / 100, 2)
                xMolybdenite2 = 0.00
                xPyrite2 = 1-xChalcopyrite2
            else:
                xGalena2 = 0.00
                xChalcopyrite2 = 0.00
                xMolybdenite2 = 0.00
                xPyrite2 = 0.00
            xGalena = round(xSulfides * xGalena2, 2)
            xChalcopyrite = round(xSulfides * xChalcopyrite2, 2)
            xMolybdenite = round(xSulfides * xMolybdenite2, 2)
            xPyrite = round(xSulfides * xPyrite2, 2)
            sumMin = round(xAnhydrite + xCalcite + xHalite + xGalena + xChalcopyrite + xMolybdenite + xPyrite, 2)
            #print("sumMin:", sumMin)
            if sumMin == 1:
                cond = True
                composition.extend([["Anh", round(xAnhydrite,2), round(chemAnhydrite[1],2)], ["Hl", round(xHalite,2), round(chemHalite[1],2)], ["Cal", round(xCalcite,2), round(chemCalcite[1],2)], ["Gn", round(xGalena,2), round(chemGalena[1],2)], ["Ccp", round(xChalcopyrite,2), round(chemChalcopyrite[1],2)], ["Mol", round(xMolybdenite,2), round(chemMolybdenite[1],2)], ["Py", round(xPyrite,2), round(chemPyrite[1],2)]])
            else:
                cond = False
        xAnhydrite = composition[0][1]
        xHalite = composition[1][1]
        xCalcite = composition[2][1]
        xGalena = composition[3][1]
        xChalcopyrite = composition[4][1]
        xMolybdenite = composition[5][1]
        xPyrite = composition[6][1]
        anhydrite.append(composition)
        mineralogy = [chemAnhydrite, chemHalite, chemCalcite, chemGalena, chemChalcopyrite, chemMolybdenite, chemPyrite]
        #
        rhoSolid = (xAnhydrite * chemAnhydrite[2] + xHalite*chemHalite[2] + xCalcite*chemCalcite[2] + xGalena*chemGalena[2] + xChalcopyrite*chemChalcopyrite[2] + xMolybdenite*chemMolybdenite[2] + xPyrite*chemPyrite[2]) / 1000
        X = [xAnhydrite, xHalite, xCalcite, xGalena, xChalcopyrite, xMolybdenite, xPyrite]
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
        phi = randint(0, 2)/100
        rho = (1 - phi) * rhoSolid + phi * chemWater[2] / 1000
        vP = (1-phi)*vP_solid + phi*chemWater[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = xAnhydrite*chemAnhydrite[5][0] + xHalite*chemHalite[5][0] + xCalcite*chemCalcite[5][0] + xGalena*chemGalena[5][0] + xChalcopyrite*chemChalcopyrite[5][0] + xMolybdenite*chemMolybdenite[5][0] + xPyrite*chemPyrite[5][0]
        PE = xAnhydrite*chemAnhydrite[5][1] + xHalite*chemHalite[5][1] + xCalcite*chemCalcite[5][1] + xGalena*chemGalena[5][1] + xChalcopyrite*chemChalcopyrite[5][1] + xMolybdenite*chemMolybdenite[5][1] + xPyrite*chemPyrite[5][1]
        #poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        #poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = xAnhydrite*chemAnhydrite[3][3] + xHalite*chemHalite[3][3] + xCalcite*chemCalcite[3][3] + xGalena*chemGalena[3][3] + xChalcopyrite*chemChalcopyrite[3][3] + xMolybdenite*chemMolybdenite[3][3] + xPyrite*chemPyrite[3][3]
        #
        anhydrite.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        anhydrite.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        anhydrite.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[4][0], 2)])
        anhydrite.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        anhydrite.append("water")
        anhydrite.append([GR, PE])
        #
        #  anhydrite = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return anhydrite
        #
class Evaporites:
    #
    def __init__(self, fluid, actualThickness):
        self.fluid = fluid
        self.actualThickness = actualThickness
        #
        self.data_anhydrite = Sulfates(data_type=True).create_anhydrite()
        self.data_gypsum = Sulfates(data_type=True).create_gypsum()
        self.data_kainite = Sulfates(data_type=True).create_kainite()
        self.data_kieserite = Sulfates(data_type=True).create_kieserite()
        self.data_celestine = Sulfates(data_type=True).create_celestine()
        self.data_calcite = Carbonates(data_type=True).create_calcite()
        self.data_dolomite = Carbonates(data_type=True).create_dolomite()
        self.data_magnesite = Carbonates(data_type=True).create_magnesite()
        self.data_halite = Halides(dict=True).create_halite()
        self.data_sylvite = Halides(dict=True).create_sylvite()
        self.data_carnallite = Halides(dict=True).create_carnallite()
        self.data_water = Water.water("")
    #
    def create_simple_rocksalt(self, w_Na=None, w_Cl=None, amounts=None, porosity=None, dict=False):
        #
        results = {}
        results["rock"] = "Rock Salt"
        #
        self.w_Na = w_Na
        self.w_Cl = w_Cl
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        halite = minerals.halides.halite("")
        anhydrite = minerals.sulfates.anhydrite("")
        gypsum = minerals.sulfates.gypsum("")
        sylvite = minerals.halides.sylvite("")
        #
        mineralogy = [halite, anhydrite, gypsum, sylvite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Na == None and self.w_Cl == None and self.amounts == None:
                w_hl = round(rd.uniform(0.90, 1), 4)
                w_acc = round((1-w_hl), 4)
                w_anh = round(w_acc*rd.uniform(0, 1), 4)
                w_gp = round(w_acc*rd.uniform(0, (1-w_anh)), 4)
                w_syl = round((1-w_hl-w_anh-w_gp), 4)
            elif self.w_Na != None:
                w_hl = round((self.w_Na)/(halite[6][0]), 4)
                w_acc = round((1-w_hl), 4)
                w_anh = round(abs(w_acc*rd.uniform(0, 1)), 4)
                w_gp = round(abs(w_acc*rd.uniform(0, (1-w_anh))), 4)
                w_syl = round(abs(w_acc*(1-w_anh-w_gp)), 4)
            elif self.w_Cl != None:
                w_acc = round(abs(rd.uniform(0.0, 0.1)), 4)
                w_syl = round(abs(w_acc*rd.uniform(0, 1)), 4)
                w_anh = round(abs(w_acc*rd.uniform(0, (1-w_syl))), 4)
                w_gp = round(abs(w_acc*(1-w_syl-w_anh)), 4)
                w_hl = round((self.w_Cl-w_syl*sylvite[6][0])/(halite[6][0]), 4)
            elif type(self.amounts) is list:
                w_hl = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_anh = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_gp = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_syl = round(1-w_hl-w_anh-w_gp, 4)
            #
            if w_hl >= 0.0 and w_anh >= 0.0 and w_gp >= 0.0 and w_syl >= 0.0:
                sumMin = round(w_hl + w_anh + w_gp + w_syl, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_gp*gypsum[6][0], 4)
            w_O = round(w_anh*anhydrite[6][0] + w_gp*gypsum[6][1], 4)
            w_Na = round(w_hl*halite[6][0], 4)
            w_S = round(w_anh*anhydrite[6][1] + w_gp*gypsum[6][2], 4)
            w_Cl = round(w_hl*halite[6][1] + w_syl*sylvite[6][0], 4)
            w_K = round(w_syl*sylvite[6][1], 4)
            w_Ca = round(w_anh*anhydrite[6][2] + w_gp*gypsum[6][3], 4)
            sumConc = round(w_H + w_O + w_Na + w_S + w_Cl + w_K + w_Ca, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Hl", "Anh", "Gp", "Syl"]))
                concentrations = [w_H, w_O, w_Na, w_S, w_Cl, w_K, w_Ca]
                amounts = [w_hl, w_anh, w_gp, w_syl]
            else:
                cond = False
        #
        element_list = ["H", "O", "Na", "S", "Cl", "K", "Ca"]
        mineral_list = ["Hl", "Anh", "Gp", "Syl"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_hl*halite[2] + w_anh*anhydrite[2] + w_gp*gypsum[2] + w_syl*sylvite[2]) / 1000
        X = [w_hl, w_anh, w_gp, w_syl]
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
        GR = w_hl*halite[5][0] + w_anh*anhydrite[5][0] + w_gp*gypsum[5][0] + w_syl*sylvite[5][0]
        PE = w_hl*halite[5][1] + w_anh*anhydrite[5][1] + w_gp*gypsum[5][1] + w_syl*sylvite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_hl*halite[3][3] + w_anh*anhydrite[3][3] + w_gp*gypsum[3][3] + w_syl*sylvite[3][3]
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
    def create_simple_anhydrite(self, w_Ca=None, amounts=None, porosity=None, dict=False):
        #
        results = {}
        results["rock"] = "Anhydrite"
        #
        self.w_Ca = w_Ca
        self.amounts = amounts
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        anhydrite = minerals.sulfates.anhydrite("")
        calcite = minerals.carbonates.calcite("")
        dolomite = minerals.carbonates.dolomite("")
        gypsum = minerals.sulfates.gypsum("")
        halite = minerals.halides.halite("")
        #
        mineralogy = [anhydrite, calcite, dolomite, gypsum, halite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Ca == None and self.amounts == None:
                w_anh = round(rd.uniform(0.9, 1), 4)
                w_acc = round((1-w_anh), 4)
                w_cal = round(w_acc*rd.uniform(0, 1), 4)
                w_dol = round(w_acc*rd.uniform(0, (1-w_cal)), 4)
                w_gp = round(w_acc*rd.uniform(0, (1-w_cal-w_dol)), 4)
                w_hl = round(1-w_anh-w_cal-w_dol-w_gp, 4)
            elif self.w_Ca != None:
                w_acc = round(abs(rd.uniform(0, 0.2)), 4)
                w_cal = round(abs(w_acc*rd.uniform(0, 1)), 4)
                w_dol = round(abs(w_acc*rd.uniform(0, (1-w_cal))), 4)
                w_gp = round(abs(w_acc*rd.uniform(0, (1-w_cal-w_dol))), 4)
                w_hl = round(abs(w_acc*(1-w_cal-w_dol-w_gp)), 4)
                w_anh = round((self.w_Ca - w_cal*calcite[6][2] - w_dol*dolomite[6][3] - w_gp*gypsum[6][3])/(anhydrite[6][2]), 4)
            elif type(self.amounts) is list:
                w_anh = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_cal = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_dol = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_gp = round(abs(np.random.normal(self.amounts[3], 0.025)), 4)
                w_hl = round(1-w_anh-w_cal-w_dol-w_gp, 4)
            #
            if w_anh >= 0.0 and w_cal >= 0.0 and w_dol >= 0.0 and w_gp >= 0.0 and w_hl >= 0.0:
                sumMin = round(w_anh + w_cal + w_dol + w_gp + w_hl, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_gp*gypsum[6][0], 4)
            w_C = round(w_cal*calcite[6][0] + w_dol*dolomite[6][0], 4)
            w_O = round(w_anh*anhydrite[6][0] + w_cal*calcite[6][1] + w_dol*dolomite[6][1] + w_gp*gypsum[6][1], 4)
            w_Na = round(w_hl*halite[6][0], 4)
            w_Mg = round(w_dol*dolomite[6][2], 4)
            w_S = round(w_anh*anhydrite[6][1] + w_gp*gypsum[6][2], 4)
            w_Cl = round(w_hl*halite[6][1], 4)
            w_Ca = round(w_anh*anhydrite[6][2] + w_cal*calcite[6][2] + w_dol*dolomite[6][3] + w_gp*gypsum[6][3], 4)
            sumConc = round(w_H + w_C + w_O + w_Na + w_Mg + w_S + w_Cl + w_Ca, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Anh", w_anh, round(anhydrite[1], 2)], ["Cal", w_cal, round(calcite[1], 2)], ["Dol", w_dol, round(dolomite[1], 2)], ["Gp", w_gp, round(gypsum[1], 2)], ["Hl", w_hl, round(halite[1], 2)]))
                concentrations = [w_H, w_C, w_O, w_Na, w_Mg, w_S, w_Cl, w_Ca]
                amounts = [w_anh, w_cal, w_dol, w_gp, w_hl]
                # phi_V = geochemistry.Fractions.calculate_volume_fraction(self, mineralogy=mineralogy, w=amounts)
                # #print(np.around(phi_V, 4))
                # if 0.8 <= phi_V[0] <= 1.0:
                #     cond = True
                # else:
                #     composition = []
                #     cond = False
            else:
                cond = False
        #
        element_list = ["H", "C", "O", "Na", "Mg", "S", "Cl", "Ca"]
        mineral_list = ["Anh", "Cal", "Dol", "Gp", "Hl"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_anh*anhydrite[2] + w_cal*calcite[2] + w_dol*dolomite[2] + w_gp*gypsum[2] + w_hl*halite[2]) / 1000
        X = [w_anh, w_cal, w_dol, w_gp, w_hl]
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
        GR = w_anh*anhydrite[5][0] + w_cal*calcite[5][0] + w_dol*dolomite[5][0] + w_gp*gypsum[5][0] + w_hl*halite[5][0]
        PE = w_anh*anhydrite[5][1] + w_cal*calcite[5][1] + w_dol*dolomite[5][1] + w_gp*gypsum[5][1] + w_hl*halite[5][1]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_anh*anhydrite[3][3] + w_cal*calcite[3][3] + w_dol*dolomite[3][3] + w_gp*gypsum[3][3] + w_hl*halite[3][3]
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
    def create_anhydrite(self, number, porosity=None):
        #
        #
        assemblage = [self.data_anhydrite, self.data_calcite, self.data_dolomite, self.data_magnesite, self.data_halite,
                      self.data_sylvite]
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
                if mineral == "Anh":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.7, 1.0), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if 0.7 <= value <= 1.0:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Cal":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Dol":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.2), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if 0.0 <= value <= 0.2:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Mgs":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.1), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if 0.0 <= value <= 0.1:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.1), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if 0.0 <= value <= 0.1:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Syl":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.0, 0.1), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if 0.0 <= value <= 0.1:
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
                phi_helper = round(rd.uniform(porosity[0], porosity[1]), 4)
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in self.data_anhydrite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Anh"][n] * self.data_anhydrite["chemistry"][element], 4)
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
                if element in self.data_halite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Hl"][n] * self.data_halite["chemistry"][element], 4)
                    else:
                        value = round(1 - w_total, 4)
                    amounts_helper[element] += value
                    w_total += value
                if element in self.data_sylvite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = round(amounts_mineralogy["Syl"][n] * self.data_sylvite["chemistry"][element], 4)
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
                    if mineral == "Anh":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_anhydrite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_anhydrite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_anhydrite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_anhydrite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_anhydrite["PE"], 3)
                    elif mineral == "Cal":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_calcite["G"],
                            3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_calcite["PE"], 3)
                    elif mineral == "Dol":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_dolomite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_dolomite["PE"], 3)
                    elif mineral == "Mgs":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["K"], 3)
                        shearmod_helper += round(shear_factor * amounts_mineralogy[mineral][n] * self.data_magnesite["G"],
                                                 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_magnesite["PE"], 3)
                    elif mineral == "Hl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_halite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_halite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_halite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_halite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_halite["PE"], 3)
                    elif mineral == "Syl":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * self.data_sylvite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * self.data_sylvite["K"], 3)
                        shearmod_helper += round(
                            shear_factor * amounts_mineralogy[mineral][n] * self.data_sylvite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * self.data_sylvite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * self.data_sylvite["PE"], 3)
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
        results["rock"] = "Anhydrite"
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
    def create_anhydrite_rock(self, number=1, composition=None, porosity=None):
        #
        mineralogy = {"Anh": self.data_anhydrite, "Cal": self.data_calcite, "Dol": self.data_dolomite,
                      "Mgs": self.data_magnesite, "Hl": self.data_halite, "Syl": self.data_sylvite}
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
                phi_anh = composition["Anh"]
                phi_cal = composition["Cal"]
                phi_dol = composition["Dol"]
                phi_mgs = composition["Mgs"]
                phi_hl = composition["Hl"]
                phi_syl = composition["Syl"]
                #
                phi_minerals["Anh"] = phi_anh
                phi_minerals["Cal"] = phi_cal
                phi_minerals["Dol"] = phi_dol
                phi_minerals["Mgs"] = phi_mgs
                phi_minerals["Hl"] = phi_hl
                phi_minerals["Syl"] = phi_syl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_anh = round(rd.uniform(0.67, 1.0), 4)
                    phi_cal = round(rd.uniform(0.0, (1.0 - phi_anh)), 4)
                    phi_dol = round(rd.uniform(0.0, (1.0 - phi_anh - phi_cal)), 4)
                    phi_mgs = round(rd.uniform(0.0, (1.0 - phi_anh - phi_cal - phi_dol)), 4)
                    phi_hl = round(rd.uniform(0.0, (1.0 - phi_anh - phi_cal- phi_dol - phi_mgs)), 4)
                    phi_syl = round(1 - phi_anh - phi_cal - phi_dol - phi_mgs - phi_hl, 4)
                    phi_total = phi_anh + phi_cal + phi_dol + phi_mgs + phi_hl + phi_syl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.67 <= phi_anh <= 1.0 and 0.0 <= phi_cal <= 0.2 and 0.0 <= phi_dol <= 0.2 \
                                and 0.0 <= phi_mgs <= 0.1 and 0.0 <= phi_hl <= 0.05 and 0.0 <= phi_syl <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Anh"] = phi_anh
                phi_minerals["Cal"] = phi_cal
                phi_minerals["Dol"] = phi_dol
                phi_minerals["Mgs"] = phi_mgs
                phi_minerals["Hl"] = phi_hl
                phi_minerals["Syl"] = phi_syl
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
                phi = round(rd.uniform(0.1, 0.4), 4)
            else:
                phi = round(rd.uniform(porosity[0], porosity[1]), 4)
            rho = round((1 - phi) * rho_s + phi * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Anhydrite"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = phi
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
            results["phi"] = phi
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
    def create_rocksalt(self, number=1, composition=None, porosity=None):
        #
        mineralogy = {"Anh": self.data_anhydrite, "Cal": self.data_calcite, "Dol": self.data_dolomite,
                      "Mgs": self.data_magnesite, "Hl": self.data_halite, "Syl": self.data_sylvite}
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
                phi_anh = composition["Anh"]
                phi_cal = composition["Cal"]
                phi_dol = composition["Dol"]
                phi_mgs = composition["Mgs"]
                phi_hl = composition["Hl"]
                phi_syl = composition["Syl"]
                #
                phi_minerals["Anh"] = phi_anh
                phi_minerals["Cal"] = phi_cal
                phi_minerals["Dol"] = phi_dol
                phi_minerals["Mgs"] = phi_mgs
                phi_minerals["Hl"] = phi_hl
                phi_minerals["Syl"] = phi_syl
            else:
                condition_2 = False
                while condition_2 == False:
                    phi_hl = round(rd.uniform(0.75, 1.0), 4)
                    phi_cal = round(rd.uniform(0.0, (1.0 - phi_hl)), 4)
                    phi_dol = round(rd.uniform(0.0, (1.0 - phi_hl - phi_cal)), 4)
                    phi_mgs = round(rd.uniform(0.0, (1.0 - phi_hl - phi_cal - phi_dol)), 4)
                    phi_anh = round(rd.uniform(0.0, (1.0 - phi_hl - phi_cal - phi_dol - phi_mgs)), 4)
                    phi_syl = round(1 - phi_hl - phi_cal - phi_dol - phi_mgs - phi_anh, 4)
                    phi_total = phi_hl + phi_cal + phi_dol + phi_mgs + phi_anh + phi_syl
                    #
                    if np.isclose(phi_total, 1.0000) == True:
                        if 0.67 <= phi_hl <= 1.0 and 0.0 <= phi_cal <= 0.2 and 0.0 <= phi_dol <= 0.2 \
                                and 0.0 <= phi_mgs <= 0.1 and 0.0 <= phi_anh <= 0.05 and 0.0 <= phi_syl <= 0.05:
                            condition_2 = True
                #
                phi_minerals["Anh"] = phi_anh
                phi_minerals["Cal"] = phi_cal
                phi_minerals["Dol"] = phi_dol
                phi_minerals["Mgs"] = phi_mgs
                phi_minerals["Hl"] = phi_hl
                phi_minerals["Syl"] = phi_syl
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
                phi = round(rd.uniform(0.1, 0.4), 4)
            else:
                phi = round(rd.uniform(porosity[0], porosity[1]), 4)
            rho = round((1 - phi) * rho_s + phi * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Rock Salt"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = phi
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
            results["phi"] = phi
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
    def create_potash(self, number=1, composition=None, porosity=None, dominance=None):
        #
        mineralogy = {"Hl": self.data_halite, "Syl": self.data_sylvite, "Car": self.data_carnallite,
                      "Ka": self.data_kainite, "Kie": self.data_kieserite}
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
                phi_hl = composition["Hl"]
                phi_syl = composition["Syl"]
                phi_car = composition["Car"]
                phi_Ka = composition["Ka"]
                phi_Kie = composition["Kie"]
                #
                phi_minerals["Hl"] = phi_hl
                phi_minerals["Syl"] = phi_syl
                phi_minerals["Car"] = phi_car
                phi_minerals["Ka"] = phi_Ka
                phi_minerals["Kie"] = phi_Kie
            else:
                condition_2 = False
                while condition_2 == False:
                    if dominance == None:
                        phi_hl = round(rd.uniform(0.5, 0.8), 4)
                        phi_syl = round(rd.uniform(0.15, (1.0 - phi_hl)), 4)
                        phi_car = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl)), 4)
                        phi_Ka = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl - phi_car)), 4)
                        phi_Kie = round(1 - phi_hl - phi_syl - phi_car - phi_Ka, 4)
                        phi_total = phi_hl + phi_syl + phi_car + phi_Ka + phi_Kie
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.5 <= phi_hl <= 0.8 and 0.15 <= phi_syl <= 0.3 and 0.0 <= phi_car <= 0.2 \
                                    and 0.0 <= phi_Ka <= 0.2 and 0.0 <= phi_Kie <= 0.2:
                                condition_2 = True
                    elif dominance == "Syl":
                        phi_hl = round(rd.uniform(0.3, 0.6), 4)
                        phi_syl = round(rd.uniform(0.3, (1.0 - phi_hl)), 4)
                        phi_car = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl)), 4)
                        phi_Ka = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl - phi_car)), 4)
                        phi_Kie = round(1 - phi_hl - phi_syl - phi_car - phi_Ka, 4)
                        phi_total = phi_hl + phi_syl + phi_car + phi_Ka + phi_Kie
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.3 <= phi_hl <= 0.6 and 0.3 <= phi_syl <= 0.6 and 0.0 <= phi_car <= 0.2 \
                                    and 0.0 <= phi_Ka <= 0.2 and 0.0 <= phi_Kie <= 0.2:
                                condition_2 = True
                    elif dominance == "Kie":
                        phi_hl = round(rd.uniform(0.3, 0.6), 4)
                        phi_syl = round(rd.uniform(0.0, (1.0 - phi_hl)), 4)
                        phi_car = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl)), 4)
                        phi_Ka = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl - phi_car)), 4)
                        phi_Kie = round(1 - phi_hl - phi_syl - phi_car - phi_Ka, 4)
                        phi_total = phi_hl + phi_syl + phi_car + phi_Ka + phi_Kie
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.3 <= phi_hl <= 0.6 and 0.0 <= phi_syl <= 0.2 and 0.0 <= phi_car <= 0.2 \
                                    and 0.0 <= phi_Ka <= 0.2 and 0.3 <= phi_Kie <= 0.6:
                                condition_2 = True
                    elif dominance == "Car":
                        phi_hl = round(rd.uniform(0.3, 0.6), 4)
                        phi_syl = round(rd.uniform(0.0, (1.0 - phi_hl)), 4)
                        phi_car = round(rd.uniform(0.3, (1.0 - phi_hl - phi_syl)), 4)
                        phi_Ka = round(rd.uniform(0.0, (1.0 - phi_hl - phi_syl - phi_car)), 4)
                        phi_Kie = round(1 - phi_hl - phi_syl - phi_car - phi_Ka, 4)
                        phi_total = phi_hl + phi_syl + phi_car + phi_Ka + phi_Kie
                        #
                        if np.isclose(phi_total, 1.0000) == True:
                            if 0.3 <= phi_hl <= 0.6 and 0.0 <= phi_syl <= 0.2 and 0.3 <= phi_car <= 0.6 \
                                    and 0.0 <= phi_Ka <= 0.2 and 0.3 <= phi_Kie <= 0.2:
                                condition_2 = True
                #
                phi_minerals["Hl"] = phi_hl
                phi_minerals["Syl"] = phi_syl
                phi_minerals["Car"] = phi_car
                phi_minerals["Ka"] = phi_Ka
                phi_minerals["Kie"] = phi_Kie
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
                phi = round(rd.uniform(0.1, 0.4), 4)
            else:
                phi = round(rd.uniform(porosity[0], porosity[1]), 4)
            rho = round((1 - phi) * rho_s + phi * self.data_water[2] / 1000, 3)
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
        results["rock"] = "Potash"
        if number > 1:
            results["mineralogy"] = w_minerals
            results["chemistry"] = w_elements
            results["phi"] = phi
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
            results["phi"] = phi
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