#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		evaporites.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.04.2020

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
from modules import minerals
import random as rd
from modules import fluids
from modules.geophysics import Elasticity as elast

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