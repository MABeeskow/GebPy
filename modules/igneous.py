#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		igneous.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		01.05.2020

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
from random import *
import random as rd
from modules import minerals
from modules import fluids
from modules.geophysics import Elasticity as elast

class plutonic:
    #
    def __init__(self,):
        pass
    #
    def createGranite(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemQuartz = minerals.oxides.quartz("")
        chemAlkaliFeldspar = minerals.feldspars.alkalifeldspar(self, "K")
        chemPlagioclase = minerals.feldspars.plagioclase(self, "Na")
        chemBiotite = minerals.Biotites.biotite_group(self, "Biotite")
        chemMuscovite = minerals.phyllosilicates.muscovite("")
        chem_actinolite = minerals.inosilicates.actinolite("")
        chem_tremolite = minerals.inosilicates.tremolite("")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = minerals.oxides.water("")
        #
        granite = []
        #
        cond = False
        composition = []
        while cond == False:
            xQuartz = rd.randint(24,27)/100
            xAlkalifeldspar = rd.randint(14,51)/100
            xPlagioclase = rd.randint(10,32)/100
            xBiotite = rd.randint(3,8)/100
            xMuscovite = rd.randint(5,8)/100
            xAmphibole = rd.randint(6,14)/100
            xActinolite2 = rd.randint(0,100)/100
            xTremolite2 = 1-xActinolite2
            xActinolite = xAmphibole*xActinolite2
            xTremolite = xAmphibole*xTremolite2
            sumMin = round(xQuartz + xAlkalifeldspar + xPlagioclase + xBiotite + xMuscovite + xActinolite + xTremolite, 2)
            if sumMin == 1:
                cond = True
                composition.extend([["Qz", round(xQuartz,2), round(chemQuartz[1],2)], ["Kfs", round(xAlkalifeldspar,2), round(chemAlkaliFeldspar[1][0],2), chemAlkaliFeldspar[1][1]], ["Pl", round(xPlagioclase,2), round(chemPlagioclase[1][0],2), chemPlagioclase[1][1]], ["Bt", round(xBiotite,2), round(chemBiotite[1][0],2), chemBiotite[1][1], chemBiotite[1][2]], ["Ms", round(xMuscovite,2), round(chemMuscovite[1],2)], ["Act", round(xActinolite,2), round(chem_actinolite[1][0],2), chem_actinolite[1][1]], ["Tr", round(xTremolite,2), round(chem_tremolite[1],2)]])
            else:
                cond = False
        xQuartz = composition[0][1]
        xAlkalifeldspar = composition[1][1]
        xPlagioclase = composition[2][1]
        xBiotite = composition[3][1]
        xMuscovite = composition[4][1]
        xActinolite = composition[5][1]
        xTremolite = composition[6][1]
        granite.append(composition)
        mineralogy = [chemQuartz, chemAlkaliFeldspar, chemPlagioclase, chemBiotite, chemMuscovite, chem_actinolite, chem_tremolite]
        #
        rhoSolid = (xQuartz*chemQuartz[2] + xAlkalifeldspar*chemAlkaliFeldspar[2] + xPlagioclase*chemPlagioclase[2] + xBiotite*chemBiotite[2] + xMuscovite*chemMuscovite[2] + xActinolite*chem_actinolite[2] + xTremolite*chem_tremolite[2]) / 1000
        X = [xQuartz, xAlkalifeldspar, xPlagioclase, xBiotite, xMuscovite, xActinolite, xTremolite]
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
        vP = (1-phi)*vP_solid + phi*chemWater[5]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = xQuartz*chemQuartz[5][0] + xAlkalifeldspar*chemAlkaliFeldspar[5][0] + xPlagioclase*chemPlagioclase[5][0] + xBiotite*chemBiotite[5][0] + xMuscovite*chemMuscovite[5][0] + xActinolite*chem_actinolite[5][0] + xTremolite*chem_tremolite[5][0]
        PE = xQuartz*chemQuartz[5][1] + xAlkalifeldspar*chemAlkaliFeldspar[5][1] + xPlagioclase*chemPlagioclase[5][1] + xBiotite*chemBiotite[5][1] + xMuscovite*chemMuscovite[5][1] + xActinolite*chem_actinolite[5][1] + xTremolite*chem_tremolite[5][1]
        #poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        #poisson_elastic = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
        poisson_mineralogical = xQuartz*chemQuartz[3][3] + xAlkalifeldspar*chemAlkaliFeldspar[3][3] + xPlagioclase*chemPlagioclase[3][3] + xBiotite*chemBiotite[3][3] + xMuscovite*chemMuscovite[3][3] + xActinolite*chem_actinolite[3][3] + xTremolite*chem_tremolite[3][3]
        #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
        #
        granite.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        granite.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        granite.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[4], 2)])
        granite.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        granite.append("water")
        granite.append([GR, PE])
        #granite.append([["Qz", round(chemQuartz[1][0], 2)], ["Kfs", round(chemAlkaliFeldspar[1][0], 2), chemAlkaliFeldspar[1][1]], ["Pl", round(chemPlagioclase[1][0], 2), chemPlagioclase[1][1]], ["Bt", round(chemBiotite[1][0], 2), chemBiotite[1][1], chemBiotite[1][2]], ["Ms", round(chemMuscovite[1][0], 2)], ["Act", round(chem_actinolite[1][0], 2), chem_actinolite[1][1]], ["Tr", round(chem_tremolite[1][0], 2)]])
        #
        #  granite = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return granite
        #

    def createPeridotite(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemOlivine = minerals.nesosilicates.olivine("")
        chemQuartz = minerals.oxides.quartz("")
        chemOrthoclase = minerals.tectosilicates.orthoclase("")
        chemAlbite = minerals.tectosilicates.albite("")
        chemAnorthite = minerals.tectosilicates.anorthite("")
        chemBiotite = minerals.Biotites.biotite_group(self, "Biotite")
        chemMuscovite = minerals.phyllosilicates.muscovite("")
        chemPyrite = minerals.sulfides.pyrite("")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = minerals.oxides.water("")
        #
        peridotite = []
        #
        cond = False
        composition = []
        while cond == False:
            xQuartz = round(randint(30, 36) / 100, 2)
            xOrthoclase = round(randint(24, 32) / 100, 2)
            xPlagioclase = round(randint(22, 34) / 100, 2)
            xAlbite2 = round(randint(0, 100) / 100, 2)
            xAnorthite2 = 1 - xAlbite2
            xAlbite = round(xPlagioclase * xAlbite2, 2)
            xAnorthite = round(xPlagioclase * xAnorthite2, 2)
            xMica = round(randint(4, 9) / 100, 2)
            xBiotite2 = round(randint(0, 100) / 100, 2)
            xMuscovite2 = 1 - xBiotite2
            xBiotite = round(xMica * xBiotite2, 2)
            xMuscovite = round(xMica * xMuscovite2, 2)
            xPyrite = round(randint(0, 5) / 100, 2)
            sumMin = round(xQuartz + xOrthoclase + xAlbite + xAnorthite + xBiotite + xMuscovite + xPyrite, 2)
            if sumMin == 1:
                cond = True
                composition.extend([["Qz", round(xQuartz, 2)], ["Or", round(xOrthoclase, 2)], ["Ab", round(xAlbite, 2)],
                                    ["An", round(xAnorthite, 2)], ["Bt", round(xBiotite, 2)],
                                    ["Ms", round(xMuscovite, 2)], ["Py", round(xPyrite, 2)]])
            else:
                cond = False
        xQuartz = composition[0][1]
        xOrthoclase = composition[1][1]
        xAlbite = composition[2][1]
        xAnorthite = composition[3][1]
        xBiotite = composition[4][1]
        xMuscovite = composition[5][1]
        xPyrite = composition[6][1]
        peridotite.append(composition)
        #
        rhoSolid = (xQuartz*chemQuartz[2] + xOrthoclase*chemOrthoclase[2] + xAlbite * chemAlbite[1] + xAnorthite * chemAnorthite[1] + xBiotite * chemBiotite[2] + xMuscovite * chemMuscovite[1] + xPyrite * chemPyrite[1]) / 1000
        vPSolid = xQuartz*chemQuartz[4][0] + xOrthoclase*chemOrthoclase[4][0] + xAlbite * chemAlbite[4] + xAnorthite * chemAnorthite[4] + xBiotite * chemBiotite[4][0] + xMuscovite * chemMuscovite[4] + xPyrite * chemPyrite[4]
        vSSolid = xQuartz*chemQuartz[4][1] + xOrthoclase*chemOrthoclase[4][1] + xAlbite * chemAlbite[5] + xAnorthite * chemAnorthite[5] + xBiotite * chemBiotite[4][1] + xMuscovite * chemMuscovite[5] + xPyrite * chemPyrite[5]
        #
        phi = randint(0, 2) / 100
        rho = (1 - phi) * rhoSolid + phi * chemWater[2] / 1000
        vP = (1 - phi) * vPSolid + phi * chemWater[5]
        vS = (1 - phi) * vSSolid
        shearModulus = vS ** 2 * rho
        bulkModulus = vP ** 2 * rho - 4 / 3 * shearModulus
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        phiWyllie = (chemWater[4]) / (vP) * (vPSolid - vP) / (vPSolid - chemWater[4])
        phiRaymer = [np.real(chemWater[5] / (2 * vPSolid) - ((chemWater[5] ** 2) / (4 * vPSolid ** 2) - (vP - vPSolid) / (vPSolid)) ** 0.5), np.real(chemWater[5] / (2 * vPSolid) + ((chemWater[5] ** 2) / (4 * vPSolid ** 2) - (vP - vPSolid) / (vPSolid)) ** 0.5)]
        GR = xQuartz*chemQuartz[5][0] + xOrthoclase*chemOrthoclase[5][0] + xAlbite*chemAlbite[6] + xAnorthite*chemAnorthite[6] + xBiotite*chemBiotite[5][0] + xMuscovite*chemMuscovite[6]
        #
        peridotite.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        peridotite.append([round(bulkModulus, 3), round(shearModulus, 3), round(poisson_mineralogical, 3)])
        peridotite.append([round(vP, 2), round(vS, 2), round(vPSolid, 2), round(chemWater[4], 2)])
        peridotite.append([round(phi, 3), round(phiD, 3), round(phiN, 3), round(phiWyllie, 3), phiRaymer])
        peridotite.append("water")
        peridotite.append(GR)
        #
        #  peridotite = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return peridotite
#
class volcanic:
    #
    def __init__(self,):
        pass
    #
    def createBasalt(self):
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chem_plagioclase = minerals.feldspars.plagioclase(self, "Ca")
        chem_enstatite = minerals.inosilicates.enstatite("")
        chem_ferrosilite = minerals.inosilicates.ferrosilite("")
        chem_biotite = minerals.Biotites.biotite_group(self, "Biotite")
        chem_actinolite = minerals.inosilicates.actinolite("")
        chem_tremolite = minerals.inosilicates.tremolite("")
        chem_olivine = minerals.nesosilicates.olivine(self, "Olivine")
        #
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemWater = minerals.oxides.water("")
        #
        basalt = []
        #
        cond = False
        composition = []
        while cond == False:
            xPlagioclase = rd.randint(13,54)/100
            xPyroxene = rd.randint(9,74)/100
            xEnstatite2 = rd.randint(0,100)/100
            xFerrosilite2 = 1-xEnstatite2
            xEnstatite = xPyroxene*xEnstatite2
            xFerrosilite = xPyroxene*xFerrosilite2
            xBiotite = rd.randint(0,4)/100
            xAmphibole = rd.randint(0,23)/100
            xActinolite2 = rd.randint(0,100)/100
            xTremolite2 = 1-xActinolite2
            xActinolite = xAmphibole*xActinolite2
            xTremolite = xAmphibole*xTremolite2
            xOlivine = rd.randint(0,14)/100
            sumMin = round(xPlagioclase + (xEnstatite+xFerrosilite) + xBiotite + (xActinolite+xTremolite) + xOlivine, 2)
            if sumMin == 1:
                cond = True
                composition.extend([["Pl", round(xPlagioclase,2), round(chem_plagioclase[1][0],2), chem_plagioclase[1][1]], ["En", round(xEnstatite,2), round(chem_enstatite[1],2)], ["Fs", round(xFerrosilite,2), round(chem_ferrosilite[1],2)], ["Bt", round(xBiotite,2), round(chem_biotite[1][0],2), round(chem_biotite[1][1],2), round(chem_biotite[1][2],2)], ["Act", round(xActinolite,2), round(chem_actinolite[1][0],2), chem_actinolite[1][1]], ["Tr", round(xTremolite,2), round(chem_tremolite[1],2)], ["Ol", round(xOlivine,2), round(chem_olivine[1][0],2), round(chem_olivine[1][1],2), round(chem_olivine[1][2],2), round(chem_olivine[1][3],2)]])
            else:
                cond = False
        xPlagioclase = composition[0][1]
        xEnstatite = composition[1][1]
        xFerrosilite = composition[2][1]
        xBiotite = composition[3][1]
        xActinolite = composition[4][1]
        xTremolite = composition[5][1]
        xOlivine = composition[6][1]
        mineralogy = [chem_plagioclase, chem_enstatite, chem_ferrosilite, chem_biotite, chem_actinolite, chem_tremolite, chem_olivine]
        #
        basalt.append(composition)
        #
        rhoSolid = 1.15*(xPlagioclase*chem_plagioclase[2] + xEnstatite*chem_enstatite[2] + xFerrosilite*chem_ferrosilite[2] + xBiotite*chem_biotite[2] + xActinolite*chem_actinolite[2] + xTremolite*chem_tremolite[2] + xOlivine*chem_olivine[2]) / 1000
        X = [xPlagioclase, xEnstatite, xFerrosilite, xBiotite, xActinolite, xTremolite, xOlivine]
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
        phi = randint(0, 15)/100
        rho = (1 - phi) * rhoSolid + phi * chemWater[2] / 1000
        vP = (1-phi)*vP_solid + phi*chemWater[5]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - chemWater[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = xPlagioclase*chem_plagioclase[5][0] + xEnstatite*chem_enstatite[5][0] + xFerrosilite*chem_ferrosilite[5][0] + xBiotite*chem_biotite[5][0] + xActinolite*chem_actinolite[5][0] + xTremolite*chem_tremolite[5][0] + xOlivine*chem_olivine[5][0]
        PE = xPlagioclase*chem_plagioclase[5][1] + xEnstatite*chem_enstatite[5][1] + xFerrosilite*chem_ferrosilite[5][1] + xBiotite*chem_biotite[5][1] + xActinolite*chem_actinolite[5][1] + xTremolite*chem_tremolite[5][1] + xOlivine*chem_olivine[5][1]
        #poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        #poisson_elastic = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
        poisson_mineralogical = xPlagioclase*chem_plagioclase[3][3] + xEnstatite*chem_enstatite[3][3] + xFerrosilite*chem_ferrosilite[3][3] + xBiotite*chem_biotite[3][3] + xActinolite*chem_actinolite[3][3] + xTremolite*chem_tremolite[3][3] + xOlivine*chem_olivine[3][3]
        #print("Poisson:", round(poisson_seismic,3), round(poisson_elastic,3), round(poisson_mineralogical,3))
        #
        basalt.append([round(rho, 3), round(rhoSolid, 3), round(chemWater[2] / 1000, 6)])
        basalt.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
        basalt.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(chemWater[4], 2)])
        basalt.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
        basalt.append("water")
        basalt.append([GR, PE])
        #
        #  basalt = [[mineralogical compositon], [densities], [elastic properties], [seismic velocities], [porosities], fluid name, GR]
        #
        return basalt