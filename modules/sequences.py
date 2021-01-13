#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sequences.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		21.02.2020

#-----------------------------------------------

## MODULES
import numpy as np
from numpy import round
import random as rd
from random import randint
from modules.carbonates import limestone, dolomite
from modules.siliciclastics import sandstone, shale, ore, Soil
from modules.igneous import plutonic, volcanic
from modules.evaporites import evaporites
from modules import minerals
from modules.elements import elements
from modules import fluids

class surface:
    #
    def __init__(self, surface, actualThickness):
        self.surface = surface
        self.actualThickness = actualThickness
    #
    def createSurface(self):
        # [symbol, atomic number, atomic mass, oxidation states, melting point, boiling point, density, electronegativity]
        chemH = ["H", 1, 1.0078, 1, 13.99, 20.271, 0.084, 2.2]
        chemC = elements.C(self)
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        chemSeawater = fluids.Water.seawater("")
        chemWater = [18.0146, 997, 2.08, 1444.0]
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemDolomite = minerals.carbonates.dolomite("")
        chemCalcite = minerals.carbonates.calcite("")
        chemQuartz = minerals.oxides.quartz("")
        chemOrthoclase = minerals.tectosilicates.orthoclase("")
        chemBiotite = minerals.Biotites.biotite_group(self, "Biotite")
        #
        # sequence = [lithology, thickness, top depth, bottom depth, density, vP, GR, neutron porosity]
        sequence = []
        #
        if self.surface == "water":
            thicknessUnit = randint(10, 50)
            newThickness = self.actualThickness + thicknessUnit
            #
            cond = False
            composition = []
            while cond == False:
                xWater = 1
                sumMin = xWater
                if sumMin == 1:
                    cond = True
                    composition.extend([["Seawater", round(xWater, 2)]])
                else:
                    cond = False
            xWater = composition[0][1]
            #
            phi = 1
            rhoSolid = (xWater * chemSeawater[2]) / 1000
            vPSolid = xWater * chemSeawater[4][0]
            vSSolid = xWater * chemSeawater[4][1]
            rho = (1 - phi) * rhoSolid + phi * chemSeawater[2] / 1000
            vP = (1-phi)*vPSolid + phi*chemSeawater[4][0]
            vS = (1 - phi) * vSSolid
            shearModulus = vS**2 * rho
            bulkModulus = vP**2 * rho - 4/3*shearModulus
            poisson = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
            velocities = [round(vP,2), round(vS,2), round(vPSolid,2), round(chemSeawater[4][0],2)]
            GR = chemSeawater[5][0]
            PE = 0.81
            phiD = (rhoSolid - rho) / (rhoSolid - chemSeawater[1]/1000 / 1000)
            #phiN = ((2 * phi ** 2 - phiD ** 2) ** (0.5))*100
            phiN = 100
            sequence.extend(["water", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR,1), round(phiN,1), "water", composition, round(poisson,2), round(PE,2)])
        elif self.surface == "dry sand":
            thicknessUnit = randint(10, 50)
            newThickness = self.actualThickness + thicknessUnit
            #
            cond = False
            composition = []
            while cond == False:
                xQuartz = round(randint(65, 100)/100, 2)
                xCalcite = round(randint(0, 10)/100, 2)
                xDolomite = round(randint(0, 10)/100, 2)
                xOrthoclase = round(randint(5, 30)/100, 2)
                xBiotite = round(randint(0, 5)/100, 2)
                sumMin = xQuartz + xCalcite + xDolomite + xOrthoclase + xBiotite
                if sumMin == 1:
                    cond = True
                    composition.extend([["Qz", round(xQuartz, 2)], ["Cal", round(xCalcite, 2)], ["Dol", round(xDolomite, 2)], ["Or", round(xOrthoclase, 2)], ["Bt", round(xBiotite,2), round(chemBiotite[1][0],2), round(chemBiotite[1][1],2), round(chemBiotite[1][2],2)]])
                else:
                    cond = False
            xQuartz = composition[0][1]
            xCalcite = composition[1][1]
            xDolomite = composition[2][1]
            xOrthoclase = composition[3][1]
            xBiotite = composition[4][1]
            #
            phi = randint(35, 50) / 100
            rhoSolid = (xQuartz*chemQuartz[2] + xCalcite*chemCalcite[2] + xDolomite*chemDolomite[2] + xOrthoclase*chemOrthoclase[2] + xBiotite *chemBiotite[2]) / 1000
            vPSolid = xQuartz*chemQuartz[4][0] + xCalcite*chemCalcite[4][0] + xDolomite*chemDolomite[4][0] + xOrthoclase*chemOrthoclase[4][0] + xBiotite * chemBiotite[4][0]
            vSSolid = xQuartz*chemQuartz[4][1] + xCalcite*chemCalcite[4][1] + xDolomite*chemDolomite[4][1] + xOrthoclase*chemOrthoclase[4][1] + xBiotite * chemBiotite[4][1]
            rho = (1 - phi) * rhoSolid + phi * chemWater[1] / 1000
            vP = (1-phi)*vPSolid + phi*chemWater[3]
            vS = (1 - phi) * vSSolid
            shearModulus = vS**2 * rho
            bulkModulus = vP**2 * rho - 4/3*shearModulus
            poisson = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
            velocities = [round(vP,2), round(vS,2), round(vPSolid,2), round(chemWater[3],2)]
            GR = xQuartz*chemQuartz[5][0] + xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xOrthoclase*chemOrthoclase[5][0] + xBiotite * chemBiotite[5][0]
            PE = xQuartz*chemQuartz[5][1] + xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xOrthoclase*chemOrthoclase[5][1] + xBiotite * chemBiotite[5][1]
            phiD = (rhoSolid - rho) / (rhoSolid - chemWater[1] / 1000)
            phiN = ((2 * phi ** 2 - phiD ** 2) ** (0.5))*100
            sequence.extend(["dry sand", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR,1), round(phiN,1), "water", composition, round(poisson,2), round(PE,2)])
        elif self.surface == "soil":
            thicknessUnit = randint(10, 50)
            newThickness = self.actualThickness + thicknessUnit
            #
            cond = False
            composition = []
            while cond == False:
                xOrganic = round(randint(3, 7)/100, 2)
                xQuartz = round(randint(65, 100)/100, 2)
                xCalcite = round(randint(0, 10)/100, 2)
                xOrthoclase = round(randint(5, 15)/100, 2)
                xBiotite = round(randint(0, 5)/100, 2)
                sumMin = xQuartz + xCalcite + xOrganic + xOrthoclase + xBiotite
                if sumMin == 1:
                    cond = True
                    composition.extend([["Qz", round(xQuartz, 2)], ["Cal", round(xCalcite, 2)], ["Org", round(xOrganic, 2)], ["Or", round(xOrthoclase, 2)], ["Bt", round(xBiotite,2), round(chemBiotite[1][0],2), round(chemBiotite[1][1],2), round(chemBiotite[1][2],2)]])
                else:
                    cond = False
            xQuartz = composition[0][1]
            xCalcite = composition[1][1]
            xOrganic = composition[2][1]
            xOrthoclase = composition[3][1]
            xBiotite = composition[4][1]
            #
            phi = randint(30, 50) / 100
            rhoSolid = (xQuartz*chemQuartz[2] + xCalcite*chemCalcite[2] + xOrganic * chemC[4] + xOrthoclase*chemOrthoclase[2] + xBiotite *chemBiotite[2]) / 1000
            vPSolid = xQuartz*chemQuartz[4][0] + xCalcite*chemCalcite[4][0] + xOrganic * chemC[8] + xOrthoclase*chemOrthoclase[4][0] + xBiotite * chemBiotite[4][0]
            vSSolid = xQuartz*chemQuartz[4][1] + xCalcite*chemCalcite[4][1] + xOrganic * chemC[9] + xOrthoclase*chemOrthoclase[4][1] + xBiotite * chemBiotite[4][1]
            rho = (1 - phi) * rhoSolid + phi * chemWater[1] / 1000
            vP = (1-phi)*vPSolid + phi*chemWater[3]
            vS = (1 - phi) * vSSolid
            shearModulus = vS**2 * rho
            bulkModulus = vP**2 * rho - 4/3*shearModulus
            poisson = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
            velocities = [round(vP,2), round(vS,2), round(vPSolid,2), round(chemWater[3],2)]
            GR = xQuartz*chemQuartz[5][0] + xCalcite*chemCalcite[5][0] + xOrthoclase*chemOrthoclase[5][0] + xBiotite * chemBiotite[5][0]
            PE = xQuartz*chemQuartz[5][1] + xCalcite*chemCalcite[5][1] + xOrthoclase*chemOrthoclase[5][1] + xBiotite * chemBiotite[5][1]
            phiD = (rhoSolid - rho)/(rhoSolid - (chemWater[1] + 1.2920)/1000)
            phiN = np.real((2*phi**2 - phiD**2)**(0.5))*100
            sequence.extend(["soil", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR,1), round(phiN,1), "water", composition, round(poisson,2), round(PE,2)])
        #
        return sequence
    #
class sand:
    #
    def __init__(self, actualThickness):
        self.actualThickness = actualThickness
    #
    def createWetSand(self):
        # [molar mass, density, bulk modulus, vP]
        chemWater = [18.0146, 997, 2.08, 1444]
        # [chemical formula, molar mass, density, bulk modulus, shear modulus, vP, vS]
        chemDolomite = minerals.carbonates.dolomite("")
        chemCalcite = minerals.carbonates.calcite("")
        chemQuartz = minerals.oxides.quartz("")
        chemOrthoclase = minerals.tectosilicates.orthoclase("")
        chemBiotite = minerals.Biotites.biotite_group(self, "Biotite")
        #
        # sequence = [lithology, thickness, top depth, bottom depth, density, vP, GR, neutron porosity]
        sequence = []
        #
        thicknessUnit = randint(10, 25)
        newThickness = self.actualThickness + thicknessUnit
        #
        cond = False
        composition = []
        while cond == False:
            xQuartz = round(randint(65, 100)/100, 2)
            xCalcite = round(randint(0, 10)/100, 2)
            xDolomite = round(randint(0, 10)/100, 2)
            xOrthoclase = round(randint(5, 30)/100, 2)
            xBiotite = round(randint(0, 5)/100, 2)
            sumMin = xQuartz + xCalcite + xDolomite + xOrthoclase + xBiotite
            if sumMin == 1:
                cond = True
                composition.extend([["Qz", round(xQuartz,2), round(chemQuartz[1],2)], ["Cal", round(xCalcite,2), round(chemCalcite[1],2)], ["Dol", round(xDolomite,2), round(chemDolomite[1],2)], ["Or", round(xOrthoclase,2), round(chemOrthoclase[1],2)], ["Bt", round(xBiotite,2), round(chemBiotite[1][0],2), round(chemBiotite[1][1],2), round(chemBiotite[1][2],2)]])
            else:
                cond = False
        xQuartz = composition[0][1]
        xCalcite = composition[1][1]
        xDolomite = composition[2][1]
        xOrthoclase = composition[3][1]
        xBiotite = composition[4][1]
        #
        phi = randint(25, 40) / 100
        rhoSolid = (xQuartz*chemQuartz[2] + xCalcite*chemCalcite[2] + xDolomite*chemDolomite[2] + xOrthoclase*chemOrthoclase[2] + xBiotite *chemBiotite[2]) / 1000
        vPSolid = xQuartz*chemQuartz[4][0] + xCalcite*chemCalcite[4][0] + xDolomite*chemDolomite[4][0] + xOrthoclase*chemOrthoclase[4][0] + xBiotite * chemBiotite[4][0]
        vSSolid = xQuartz*chemQuartz[4][1] + xCalcite*chemCalcite[4][1] + xDolomite*chemDolomite[4][1] + xOrthoclase*chemOrthoclase[4][1] + xBiotite * chemBiotite[4][1]
        rho = (1 - phi) * rhoSolid + phi * chemWater[1] / 1000
        vP = (1-phi)*vPSolid + phi*chemWater[3]
        vS = (1 - phi) * vSSolid
        shearModulus = vS**2 * rho
        bulkModulus = vP**2 * rho - 4/3*shearModulus
        poisson = (3*bulkModulus - 2*shearModulus)/(6*bulkModulus + 2*shearModulus)
        velocities = [round(vP,2), round(vS,2), round(vPSolid,2), round(chemWater[3],2)]
        GR = xQuartz*chemQuartz[5][0] + xCalcite*chemCalcite[5][0] + xDolomite*chemDolomite[5][0] + xOrthoclase*chemOrthoclase[5][0] + xBiotite*chemBiotite[5][0]
        PE = xQuartz*chemQuartz[5][1] + xCalcite*chemCalcite[5][1] + xDolomite*chemDolomite[5][1] + xOrthoclase*chemOrthoclase[5][1] + xBiotite*chemBiotite[5][1]
        phiD = (rhoSolid - rho)/(rhoSolid - chemWater[1]/1000)
        phiN = np.real((2*phi**2 - phiD**2)**(0.5))*100
        sequence.extend(["wet sand", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR,1), round(phiN,1), "water", composition, round(poisson,2), round(PE,2)])
        #
        return sequence
#
class subsurface:
    #
    def __init__(self, actualThickness, rockType, parts):
        self.actualThickness = actualThickness
        self.rockType = rockType
        self.parts = parts
    #
    def createSiliciclastics(self):
        #
        sequence = []
        minThickness = 5
        maxThickness = 50
        #
        if self.rockType == "shale" or self.rockType == "halite" or self.rockType == "anhydrite":
            magicnumber = randint(0, 2)
            #########
            # water #
            #########
            if magicnumber == 0:
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = sandstone("water", top)
                    dataSandstone = data.createSandstone()
                    rho = round(dataSandstone[1][0], 3)
                    vP = round(dataSandstone[3][0], 1)
                    vS = round(dataSandstone[3][1], 1)
                    velocities = dataSandstone[3]
                    phiN = round(dataSandstone[4][2] * 100, 1)
                    GR = round(dataSandstone[6][0], 1)
                    PE = round(dataSandstone[6][1], 1)
                    fluid = dataSandstone[5]
                    composition = dataSandstone[0]
                    bulkModulus = dataSandstone[2][0]
                    shearModulus = dataSandstone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataSandstone[2][3]
                    sequence.append(["sandstone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
            ###############
            # oil - water #
            ###############
            elif magicnumber == 1:
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = sandstone("oil", top)
                    dataSandstone = data.createSandstone()
                    rho = round(dataSandstone[1][0], 3)
                    vP = round(dataSandstone[3][0], 1)
                    vS = round(dataSandstone[3][1], 1)
                    velocities = dataSandstone[3]
                    phiN = round(dataSandstone[4][2] * 100, 1)
                    GR = round(dataSandstone[6][0], 1)
                    PE = round(dataSandstone[6][1], 1)
                    fluid = dataSandstone[5]
                    composition = dataSandstone[0]
                    bulkModulus = dataSandstone[2][0]
                    shearModulus = dataSandstone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataSandstone[2][3]
                    sequence.append(["sandstone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
                #
                thicknessUnit2 = randint(minThickness, maxThickness)
                N = self.parts
                d2 = float(thicknessUnit2/N)
                for i in range(0, N):
                    top2 = bottom + i*d2
                    bottom2 = top2 + d2
                    data = sandstone("water", top2)
                    dataSandstone = data.createSandstone()
                    rho = round(dataSandstone[1][0], 3)
                    vP = round(dataSandstone[3][0], 1)
                    vS = round(dataSandstone[3][1], 1)
                    velocities = dataSandstone[3]
                    phiN = round(dataSandstone[4][2] * 100, 1)
                    GR = round(dataSandstone[6][0], 1)
                    PE = round(dataSandstone[6][1], 1)
                    fluid = dataSandstone[5]
                    composition = dataSandstone[0]
                    bulkModulus = dataSandstone[2][0]
                    shearModulus = dataSandstone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataSandstone[2][3]
                    sequence.append(["sandstone", round(d2,1), round(top2,1), round(bottom2,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
            #####################
            # gas - oil - water #
            #####################
            elif magicnumber == 2:
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = sandstone("gas", top)
                    dataSandstone = data.createSandstone()
                    rho = round(dataSandstone[1][0], 3)
                    vP = round(dataSandstone[3][0], 1)
                    vS = round(dataSandstone[3][1], 1)
                    velocities = dataSandstone[3]
                    phiN = round(dataSandstone[4][2] * 100, 1)
                    GR = round(dataSandstone[6][0], 1)
                    PE = round(dataSandstone[6][1], 1)
                    fluid = dataSandstone[5]
                    composition = dataSandstone[0]
                    bulkModulus = dataSandstone[2][0]
                    shearModulus = dataSandstone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataSandstone[2][3]
                    sequence.append(["sandstone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
                #
                thicknessUnit2 = randint(minThickness, maxThickness)
                N = self.parts
                d2 = float(thicknessUnit2/N)
                for i in range(0, N):
                    top2 = bottom + i*d2
                    bottom2 = top2 + d2
                    data = sandstone("oil", top2)
                    dataSandstone = data.createSandstone()
                    rho = round(dataSandstone[1][0], 3)
                    vP = round(dataSandstone[3][0], 1)
                    vS = round(dataSandstone[3][1], 1)
                    velocities = dataSandstone[3]
                    phiN = round(dataSandstone[4][2] * 100, 1)
                    GR = round(dataSandstone[6][0], 1)
                    PE = round(dataSandstone[6][1], 1)
                    fluid = dataSandstone[5]
                    composition = dataSandstone[0]
                    bulkModulus = dataSandstone[2][0]
                    shearModulus = dataSandstone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataSandstone[2][3]
                    sequence.append(["sandstone", round(d2,1), round(top2,1), round(bottom2,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
                    #
                thicknessUnit3 = randint(minThickness, maxThickness)
                N = self.parts
                d3 = float(thicknessUnit3/N)
                for i in range(0, N):
                    top3 = bottom2 + i*d3
                    bottom3 = top3 + d3
                    data = sandstone("water", top3)
                    dataSandstone = data.createSandstone()
                    rho = round(dataSandstone[1][0], 3)
                    vP = round(dataSandstone[3][0], 1)
                    vS = round(dataSandstone[3][1], 1)
                    velocities = dataSandstone[3]
                    phiN = round(dataSandstone[4][2] * 100, 1)
                    GR = round(dataSandstone[6][0], 1)
                    PE = round(dataSandstone[6][1], 1)
                    fluid = dataSandstone[5]
                    composition = dataSandstone[0]
                    bulkModulus = dataSandstone[2][0]
                    shearModulus = dataSandstone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataSandstone[2][3]
                    sequence.append(["sandstone", round(d3,1), round(top3,1), round(bottom3,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
        else:
            #########
            # water #
            #########
            thicknessUnit = randint(minThickness, maxThickness)
            #
            N = self.parts
            d = float(thicknessUnit/N)
            for i in range(0, N):
                top = self.actualThickness + i*d
                bottom = top + d
                data = sandstone("water", top)
                dataSandstone = data.createSandstone()
                rho = round(dataSandstone[1][0], 3)
                vP = round(dataSandstone[3][0], 1)
                vS = round(dataSandstone[3][1], 1)
                velocities = dataSandstone[3]
                phiN = round(dataSandstone[4][2] * 100, 1)
                GR = round(dataSandstone[6][0], 1)
                PE = round(dataSandstone[6][1], 1)
                fluid = dataSandstone[5]
                composition = dataSandstone[0]
                bulkModulus = dataSandstone[2][0]
                shearModulus = dataSandstone[2][1]
                #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                poisson = dataSandstone[2][3]
                sequence.append(["sandstone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
        #
        return sequence
    #
    def createCarbonates(self, keyword="carbonates"):
        #
        sequence = []
        minThickness = 15
        maxThickness = 60
        self.keyword = keyword
        #
        magicnumber = randint(0, 1)
        if magicnumber == 0 or self.keyword == "limestone":
            if self.rockType == "shale" or self.rockType == "halite" or self.rockType == "anhydrite":
                magicnumber = randint(0, 2)
                #########
                # water #
                #########
                if magicnumber == 0:
                    thicknessUnit = randint(minThickness, maxThickness)
                    N = self.parts
                    d = float(thicknessUnit/N)
                    for i in range(0, N):
                        top = self.actualThickness + i*d
                        bottom = top + d
                        data = limestone("water", top)
                        dataLimestone = data.createLimestone()
                        rho = round(dataLimestone[1][0], 3)
                        vP = round(dataLimestone[3][0], 1)
                        vS = round(dataLimestone[3][1], 1)
                        velocities = dataLimestone[3]
                        phiN = round(dataLimestone[4][2] * 100, 1)
                        GR = round(dataLimestone[6][0], 1)
                        PE = round(dataLimestone[6][1], 1)
                        fluid = dataLimestone[5]
                        composition = dataLimestone[0]
                        bulkModulus = dataLimestone[2][0]
                        shearModulus = dataLimestone[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataLimestone[2][3]
                        sequence.append(["limestone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                ###############
                # oil - water #
                ###############
                elif magicnumber == 1:
                    thicknessUnit = randint(minThickness, maxThickness)
                    N = self.parts
                    d = float(thicknessUnit/N)
                    for i in range(0, N):
                        top = self.actualThickness + i*d
                        bottom = top + d
                        data = limestone("oil", top)
                        dataLimestone = data.createLimestone()
                        rho = round(dataLimestone[1][0], 3)
                        vP = round(dataLimestone[3][0], 1)
                        vS = round(dataLimestone[3][1], 1)
                        velocities = dataLimestone[3]
                        phiN = round(dataLimestone[4][2] * 100, 1)
                        GR = round(dataLimestone[6][0], 1)
                        PE = round(dataLimestone[6][1], 1)
                        fluid = dataLimestone[5]
                        composition = dataLimestone[0]
                        bulkModulus = dataLimestone[2][0]
                        shearModulus = dataLimestone[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataLimestone[2][3]
                        sequence.append(["limestone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                    #
                    thicknessUnit2 = randint(minThickness, maxThickness)
                    N = self.parts
                    d2 = float(thicknessUnit2/N)
                    for i in range(0, N):
                        top2 = bottom + i*d2
                        bottom2 = top2 + d2
                        data = limestone("water", top2)
                        dataLimestone = data.createLimestone()
                        rho = round(dataLimestone[1][0], 3)
                        vP = round(dataLimestone[3][0], 1)
                        vS = round(dataLimestone[3][1], 1)
                        velocities = dataLimestone[3]
                        phiN = round(dataLimestone[4][2] * 100, 1)
                        GR = round(dataLimestone[6][0], 1)
                        PE = round(dataLimestone[6][1], 1)
                        fluid = dataLimestone[5]
                        composition = dataLimestone[0]
                        bulkModulus = dataLimestone[2][0]
                        shearModulus = dataLimestone[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataLimestone[2][3]
                        sequence.append(["limestone", round(d2,1), round(top2,1), round(bottom2,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                #####################
                # gas - oil - water #
                #####################
                elif magicnumber == 2:
                    thicknessUnit = randint(minThickness, maxThickness)
                    N = self.parts
                    d = float(thicknessUnit/N)
                    for i in range(0, N):
                        top = self.actualThickness + i*d
                        bottom = top + d
                        data = limestone("gas", top)
                        dataLimestone = data.createLimestone()
                        rho = round(dataLimestone[1][0], 3)
                        vP = round(dataLimestone[3][0], 1)
                        vS = round(dataLimestone[3][1], 1)
                        velocities = dataLimestone[3]
                        phiN = round(dataLimestone[4][2] * 100, 1)
                        GR = round(dataLimestone[6][0], 1)
                        PE = round(dataLimestone[6][1], 1)
                        fluid = dataLimestone[5]
                        composition = dataLimestone[0]
                        bulkModulus = dataLimestone[2][0]
                        shearModulus = dataLimestone[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataLimestone[2][3]
                        sequence.append(["limestone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                    #
                    thicknessUnit2 = randint(minThickness, maxThickness)
                    N = self.parts
                    d2 = float(thicknessUnit2/N)
                    for i in range(0, N):
                        top2 = bottom + i*d2
                        bottom2 = top2 + d2
                        data = limestone("oil", top2)
                        dataLimestone = data.createLimestone()
                        rho = round(dataLimestone[1][0], 3)
                        vP = round(dataLimestone[3][0], 1)
                        vS = round(dataLimestone[3][1], 1)
                        velocities = dataLimestone[3]
                        phiN = round(dataLimestone[4][2] * 100, 1)
                        GR = round(dataLimestone[6][0], 1)
                        PE = round(dataLimestone[6][1], 1)
                        fluid = dataLimestone[5]
                        composition = dataLimestone[0]
                        bulkModulus = dataLimestone[2][0]
                        shearModulus = dataLimestone[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataLimestone[2][3]
                        sequence.append(["limestone", round(d2,1), round(top2,1), round(bottom2,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                        #
                    thicknessUnit3 = randint(minThickness, maxThickness)
                    N = self.parts
                    d3 = float(thicknessUnit3/N)
                    for i in range(0, N):
                        top3 = bottom2 + i*d3
                        bottom3 = top3 + d3
                        data = limestone("water", top3)
                        dataLimestone = data.createLimestone()
                        rho = round(dataLimestone[1][0], 3)
                        vP = round(dataLimestone[3][0], 1)
                        vS = round(dataLimestone[3][1], 1)
                        velocities = dataLimestone[3]
                        phiN = round(dataLimestone[4][2] * 100, 1)
                        GR = round(dataLimestone[6][0], 1)
                        PE = round(dataLimestone[6][1], 1)
                        fluid = dataLimestone[5]
                        composition = dataLimestone[0]
                        bulkModulus = dataLimestone[2][0]
                        shearModulus = dataLimestone[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataLimestone[2][3]
                        sequence.append(["limestone", round(d3,1), round(top3,1), round(bottom3,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
            else:
                #########
                # water #
                #########
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = limestone("water", top)
                    dataLimestone = data.createLimestone()
                    rho = round(dataLimestone[1][0], 3)
                    vP = round(dataLimestone[3][0], 1)
                    vS = round(dataLimestone[3][1], 1)
                    velocities = dataLimestone[3]
                    phiN = round(dataLimestone[4][2] * 100, 1)
                    GR = round(dataLimestone[6][0], 1)
                    PE = round(dataLimestone[6][1], 1)
                    fluid = dataLimestone[5]
                    composition = dataLimestone[0]
                    bulkModulus = dataLimestone[2][0]
                    shearModulus = dataLimestone[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataLimestone[2][3]
                    sequence.append(["limestone", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                #
        elif magicnumber == 1 or self.keyword == "dolomite":
            if self.rockType == "shale" or self.rockType == "halite" or self.rockType == "anhydrite":
                magicnumber = randint(0, 2)
                #########
                # water #
                #########
                if magicnumber == 0:
                    thicknessUnit = randint(minThickness, maxThickness)
                    N = self.parts
                    d = float(thicknessUnit/N)
                    for i in range(0, N):
                        top = self.actualThickness + i*d
                        bottom = top + d
                        data = dolomite("water", top)
                        dataDolomite = data.createDolomite()
                        rho = round(dataDolomite[1][0], 3)
                        vP = round(dataDolomite[3][0], 1)
                        vS = round(dataDolomite[3][1], 1)
                        velocities = dataDolomite[3]
                        phiN = round(dataDolomite[4][2] * 100, 1)
                        GR = round(dataDolomite[6][0], 1)
                        PE = round(dataDolomite[6][1], 1)
                        fluid = dataDolomite[5]
                        composition = dataDolomite[0]
                        bulkModulus = dataDolomite[2][0]
                        shearModulus = dataDolomite[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataDolomite[2][3]
                        sequence.append(["dolomite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                ###############
                # oil - water #
                ###############
                elif magicnumber == 1:
                    thicknessUnit = randint(minThickness, maxThickness)
                    N = self.parts
                    d = float(thicknessUnit/N)
                    for i in range(0, N):
                        top = self.actualThickness + i*d
                        bottom = top + d
                        data = dolomite("oil", top)
                        dataDolomite = data.createDolomite()
                        rho = round(dataDolomite[1][0], 3)
                        vP = round(dataDolomite[3][0], 1)
                        vS = round(dataDolomite[3][1], 1)
                        velocities = dataDolomite[3]
                        phiN = round(dataDolomite[4][2] * 100, 1)
                        GR = round(dataDolomite[6][0], 1)
                        PE = round(dataDolomite[6][1], 1)
                        fluid = dataDolomite[5]
                        composition = dataDolomite[0]
                        bulkModulus = dataDolomite[2][0]
                        shearModulus = dataDolomite[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataDolomite[2][3]
                        sequence.append(["dolomite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                    #
                    thicknessUnit2 = randint(minThickness, maxThickness)
                    N = self.parts
                    d2 = float(thicknessUnit2/N)
                    for i in range(0, N):
                        top2 = bottom + i*d2
                        bottom2 = top2 + d2
                        data = dolomite("water", top2)
                        dataDolomite = data.createDolomite()
                        rho = round(dataDolomite[1][0], 3)
                        vP = round(dataDolomite[3][0], 1)
                        vS = round(dataDolomite[3][1], 1)
                        velocities = dataDolomite[3]
                        phiN = round(dataDolomite[4][2] * 100, 1)
                        GR = round(dataDolomite[6][0], 1)
                        PE = round(dataDolomite[6][1], 1)
                        fluid = dataDolomite[5]
                        composition = dataDolomite[0]
                        bulkModulus = dataDolomite[2][0]
                        shearModulus = dataDolomite[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataDolomite[2][3]
                        sequence.append(["dolomite", round(d2,1), round(top2,1), round(bottom2,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                #####################
                # gas - oil - water #
                #####################
                elif magicnumber == 2:
                    thicknessUnit = randint(minThickness, maxThickness)
                    N = self.parts
                    d = float(thicknessUnit/N)
                    for i in range(0, N):
                        top = self.actualThickness + i*d
                        bottom = top + d
                        data = dolomite("gas", top)
                        dataDolomite = data.createDolomite()
                        rho = round(dataDolomite[1][0], 3)
                        vP = round(dataDolomite[3][0], 1)
                        vS = round(dataDolomite[3][1], 1)
                        velocities = dataDolomite[3]
                        phiN = round(dataDolomite[4][2] * 100, 1)
                        GR = round(dataDolomite[6][0], 1)
                        PE = round(dataDolomite[6][1], 1)
                        fluid = dataDolomite[5]
                        composition = dataDolomite[0]
                        bulkModulus = dataDolomite[2][0]
                        shearModulus = dataDolomite[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataDolomite[2][3]
                        sequence.append(["dolomite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                    #
                    thicknessUnit2 = randint(minThickness, maxThickness)
                    N = self.parts
                    d2 = float(thicknessUnit2/N)
                    for i in range(0, N):
                        top2 = bottom + i*d2
                        bottom2 = top2 + d2
                        data = dolomite("oil", top2)
                        dataDolomite = data.createDolomite()
                        rho = round(dataDolomite[1][0], 3)
                        vP = round(dataDolomite[3][0], 1)
                        vS = round(dataDolomite[3][1], 1)
                        velocities = dataDolomite[3]
                        phiN = round(dataDolomite[4][2] * 100, 1)
                        GR = round(dataDolomite[6][0], 1)
                        PE = round(dataDolomite[6][1], 1)
                        fluid = dataDolomite[5]
                        composition = dataDolomite[0]
                        bulkModulus = dataDolomite[2][0]
                        shearModulus = dataDolomite[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataDolomite[2][3]
                        sequence.append(["dolomite", round(d2,1), round(top2,1), round(bottom2,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
                    #
                    thicknessUnit3 = randint(minThickness, maxThickness)
                    N = self.parts
                    d3 = float(thicknessUnit3/N)
                    for i in range(0, N):
                        top3 = bottom2 + i*d3
                        bottom3 = top3 + d3
                        data = dolomite("water", top3)
                        dataDolomite = data.createDolomite()
                        rho = round(dataDolomite[1][0], 3)
                        vP = round(dataDolomite[3][0], 1)
                        vS = round(dataDolomite[3][1], 1)
                        velocities = dataDolomite[3]
                        phiN = round(dataDolomite[4][2] * 100, 1)
                        GR = round(dataDolomite[6][0], 1)
                        PE = round(dataDolomite[6][1], 1)
                        fluid = dataDolomite[5]
                        composition = dataDolomite[0]
                        bulkModulus = dataDolomite[2][0]
                        shearModulus = dataDolomite[2][1]
                        #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                        poisson = dataDolomite[2][3]
                        sequence.append(["dolomite", round(d3,1), round(top3,1), round(bottom3,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
            else:
                #########
                # water #
                #########
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = dolomite("water", top)
                    dataDolomite = data.createDolomite()
                    rho = round(dataDolomite[1][0], 3)
                    vP = round(dataDolomite[3][0], 1)
                    vS = round(dataDolomite[3][1], 1)
                    velocities = dataDolomite[3]
                    phiN = round(dataDolomite[4][2] * 100, 1)
                    GR = round(dataDolomite[6][0], 1)
                    PE = round(dataDolomite[6][1], 1)
                    fluid = dataDolomite[5]
                    composition = dataDolomite[0]
                    bulkModulus = dataDolomite[2][0]
                    shearModulus = dataDolomite[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataDolomite[2][3]
                    sequence.append(["dolomite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        #
        return sequence
    #
    def createShales(self):
        #
        sequence = []
        minThickness = 10
        maxThickness = 100
        #
        thicknessUnit = randint(minThickness, maxThickness)
        N = self.parts
        d = float(thicknessUnit/N)
        for i in range(0, N):
            top = self.actualThickness + i*d
            bottom = top + d
            data = shale()
            dataShale = data.createShale()
            rho = round(dataShale[1][0], 3)
            vP = round(dataShale[3][0], 1)
            vS = round(dataShale[3][1], 1)
            velocities = dataShale[3]
            phiN = round(dataShale[4][2] * 100, 1)
            GR = round(dataShale[6][0], 1)
            PE = round(dataShale[6][1], 1)
            fluid = dataShale[5]
            composition = dataShale[0]
            bulkModulus = dataShale[2][0]
            shearModulus = dataShale[2][1]
            #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
            poisson = dataShale[2][3]
            sequence.append(["shale", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        #
        return sequence
    #
    def createEvaporites(self, keyword):
        #
        sequence = []
        minThickness = 25
        maxThickness = 150
        #
        magicnumber = rd.randint(0, 1)
        #
        if magicnumber == 0 and keyword != "Anhydrite":
            thicknessUnit = randint(minThickness, maxThickness)
            N = self.parts
            d = float(thicknessUnit/N)
            for i in range(0, N):
                top = self.actualThickness + i*d
                bottom = top + d
                data = evaporites()
                dataHalite = data.createRocksalt()
                rho = round(dataHalite[1][0], 3)
                vP = round(dataHalite[3][0], 1)
                vS = round(dataHalite[3][1], 1)
                velocities = dataHalite[3]
                phiN = round(dataHalite[4][2] * 100, 1)
                GR = round(dataHalite[6][0], 1)
                PE = round(dataHalite[6][1], 1)
                fluid = dataHalite[5]
                composition = dataHalite[0]
                bulkModulus = dataHalite[2][0]
                shearModulus = dataHalite[2][1]
                #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                poisson = dataHalite[2][3]
                sequence.append(["halite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        elif magicnumber == 1 or keyword == "Anhydrite":
            thicknessUnit = randint(minThickness, maxThickness)
            N = self.parts
            d = float(thicknessUnit/N)
            for i in range(0, N):
                top = self.actualThickness + i*d
                bottom = top + d
                data = evaporites()
                dataAnhydrite = data.createAnhydrite()
                rho = round(dataAnhydrite[1][0], 3)
                vP = round(dataAnhydrite[3][0], 1)
                vS = round(dataAnhydrite[3][1], 1)
                velocities = dataAnhydrite[3]
                phiN = round(dataAnhydrite[4][2] * 100, 1)
                GR = round(dataAnhydrite[6][0], 1)
                PE = round(dataAnhydrite[6][1], 1)
                fluid = dataAnhydrite[5]
                composition = dataAnhydrite[0]
                bulkModulus = dataAnhydrite[2][0]
                shearModulus = dataAnhydrite[2][1]
                #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                poisson = dataAnhydrite[2][3]
                sequence.append(["anhydrite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        #
        return sequence
    #
    def createIgneousRocks(self, rock="None"):
        #
        sequence = []
        minThickness = 15
        maxThickness = 60
        #
        if rock == "None":
            magicnumber = rd.randint(0, 1)
            #
            if magicnumber == 0:
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = plutonic()
                    dataGranite = data.createGranite()
                    rho = round(dataGranite[1][0], 3)
                    vP = round(dataGranite[3][0], 1)
                    vS = round(dataGranite[3][1], 1)
                    velocities = dataGranite[3]
                    phiN = round(dataGranite[4][2] * 100, 1)
                    GR = round(dataGranite[6][0], 1)
                    PE = round(dataGranite[6][1], 1)
                    fluid = dataGranite[5]
                    composition = dataGranite[0]
                    bulkModulus = dataGranite[2][0]
                    shearModulus = dataGranite[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataGranite[2][3]
                    sequence.append(["granite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
            elif magicnumber == 1:
                thicknessUnit = randint(minThickness, maxThickness)
                N = self.parts
                d = float(thicknessUnit/N)
                for i in range(0, N):
                    top = self.actualThickness + i*d
                    bottom = top + d
                    data = volcanic()
                    dataBasalt = data.createBasalt()
                    rho = round(dataBasalt[1][0], 3)
                    vP = round(dataBasalt[3][0], 1)
                    vS = round(dataBasalt[3][1], 1)
                    velocities = dataBasalt[3]
                    phiN = round(dataBasalt[4][2] * 100, 1)
                    GR = round(dataBasalt[6][0], 1)
                    PE = round(dataBasalt[6][1], 1)
                    fluid = dataBasalt[5]
                    composition = dataBasalt[0]
                    bulkModulus = dataBasalt[2][0]
                    shearModulus = dataBasalt[2][1]
                    #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                    poisson = dataBasalt[2][3]
                    sequence.append(["basalt", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        elif rock == "granite":
            thicknessUnit = randint(minThickness, maxThickness)
            N = self.parts
            d = float(thicknessUnit/N)
            for i in range(0, N):
                top = self.actualThickness + i*d
                bottom = top + d
                data = plutonic()
                dataGranite = data.createGranite()
                rho = round(dataGranite[1][0], 3)
                vP = round(dataGranite[3][0], 1)
                vS = round(dataGranite[3][1], 1)
                velocities = dataGranite[3]
                phiN = round(dataGranite[4][2] * 100, 1)
                GR = round(dataGranite[6][0], 1)
                PE = round(dataGranite[6][1], 1)
                fluid = dataGranite[5]
                composition = dataGranite[0]
                bulkModulus = dataGranite[2][0]
                shearModulus = dataGranite[2][1]
                #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                poisson = dataGranite[2][3]
                sequence.append(["granite", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        elif rock == "basalt":
            thicknessUnit = randint(minThickness, maxThickness)
            N = self.parts
            d = float(thicknessUnit/N)
            for i in range(0, N):
                top = self.actualThickness + i*d
                bottom = top + d
                data = volcanic()
                dataBasalt = data.createBasalt()
                rho = round(dataBasalt[1][0], 3)
                vP = round(dataBasalt[3][0], 1)
                vS = round(dataBasalt[3][1], 1)
                velocities = dataBasalt[3]
                phiN = round(dataBasalt[4][2] * 100, 1)
                GR = round(dataBasalt[6][0], 1)
                PE = round(dataBasalt[6][1], 1)
                fluid = dataBasalt[5]
                composition = dataBasalt[0]
                bulkModulus = dataBasalt[2][0]
                shearModulus = dataBasalt[2][1]
                #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
                poisson = dataBasalt[2][3]
                sequence.append(["basalt", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        #
        return sequence
        #
    def create_ores(self):
        #
        sequence = []
        minThickness = 15
        maxThickness = 45
        #
        ## Fe-bearing ore
        thicknessUnit = randint(minThickness, maxThickness)
        N = self.parts
        d = float(thicknessUnit/N)
        for i in range(0, N):
            top = self.actualThickness + i*d
            bottom = top + d
            data = ore()
            data_ore_Fe = data.create_ore_Fe()
            rho = round(data_ore_Fe[1][0], 3)
            vP = round(data_ore_Fe[3][0], 1)
            vS = round(data_ore_Fe[3][1], 1)
            velocities = data_ore_Fe[3]
            phiN = round(data_ore_Fe[4][2] * 100, 1)
            GR = round(data_ore_Fe[6][0], 1)
            PE = round(data_ore_Fe[6][1], 1)
            fluid = data_ore_Fe[5]
            composition = data_ore_Fe[0]
            bulkModulus = data_ore_Fe[2][0]
            shearModulus = data_ore_Fe[2][1]
            #poisson = (vP**2 - 2*vS**2)/(2*(vP**2 - vS**2))
            poisson = data_ore_Fe[2][3]
            sequence.append(["ore", round(d,1), round(top,1), round(bottom,1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        #
        return sequence
    #
class MineralDeposits:
    #
    def __init__(self):
        pass
    #
    def create_vms_deposit(self):
        s = 0
#
class SedimentaryBasin:
    #
    def __init__(self, actualThickness, parts):
        self.actualThickness = actualThickness
        self.parts = parts
    #
    def create_soil(self):
        sequence = []
        thicknessUnit = rd.randint(5, 15)
        newThickness = self.actualThickness + thicknessUnit
        #
        N = self.parts
        d = float(thicknessUnit/N)
        #
        for i in range(N):
            data = Soil()
            data_soil = data.create_simple_soil()
            top = self.actualThickness + i*d
            bottom = top + d
            sequence.append(["soil", round(d,1), round(top,1), round(bottom,1), data_soil[1:]])