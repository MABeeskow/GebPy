#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sequences.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.05.2023

#-----------------------------------------------

## MODULES
import sys
import numpy as np
from numpy import round
import random as rd
from random import randint
from modules.carbonates import limestone, dolomite
from modules.siliciclastics import sandstone, shale, Soil
from modules.igneous import plutonic, volcanic, Plutonic, Volcanic
from modules.evaporites import evaporites, Evaporites
from modules import minerals
from modules.elements import elements
from modules import fluids

class DataProcessing:
    #
    def __init__(self, dataset):
        self.dataset = dataset
    #
    def extract_lithology(self, type="sequence"):
        """Returns a list that contains the lithology of a previously generated sequence or the name of a rock.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of rock names
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][0])
        else:
            data.append(self.dataset[0][0])
        #
        return data
    #
    def extract_thickness(self, type="sequence"):
        """Returns a list that contains the thickness of previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of thicknesses
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][1])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][1])
        #
        return np.array(data)
    #
    def extract_top(self, type="sequence"):
        """Returns a list that contains the top depths of previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of top depths
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][2])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][2])
        #
        return np.array(data)
    #
    def extract_bottom(self, type="sequence"):
        """Returns a list that contains the bottom depths of previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of bottom depths
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][3])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][3])
        #
        return np.array(data)
    #
    def extract_mineralogy(self, type="sequence"):
        """Returns a list that contains the mineralogical assemblage of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of incorporated minerals
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][0][1])
        else:
            data = self.dataset[0][4][0][1]
        #
        return data
    #
    def extract_elements(self, type="sequence"):
        """Returns a list that contains the chemical elements incorporated in the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of incorporated elements
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][0][0])
        else:
            data = self.dataset[0][4][0][0]
        #
        return data
    #
    def extract_molar_mass(self):
        """Returns a list that contains the molar mass of the previously generated minerals.
        **Arguments**:
        **Outputs**:
            data: list of molar masses
        """
        data = []
        for i in range(len(self.dataset)):
            data.append(self.dataset[i][1])
        #
        return np.array(data)
    #
    def extract_densities(self, type="sequence", keyword="bulk", dict=False):
        """Returns a list that contains the densities of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
            keyword: all, bulk, solid, fluid
        **Outputs**:
            data: list of densities
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    if keyword == "all":
                        data.append(self.dataset[i][j][4][1])
                    elif keyword == "bulk":
                        data.append(self.dataset[i][j][4][1][0])
                    elif keyword == "solid":
                        data.append(self.dataset[i][j][4][1][1])
                    elif keyword == "fluid":
                        data.append(self.dataset[i][j][4][1][2])
        elif type == "random":
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][1])
                elif keyword == "bulk":
                    data.append(self.dataset[i][1][0])
                elif keyword == "solid":
                    data.append(self.dataset[i][1][1])
                elif keyword == "fluid":
                    data.append(self.dataset[i][1][2])
        elif type == "mineral":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][2])
        elif type == "rock":
            for i in range(len(self.dataset)):
                if dict == False:
                    data.append(self.dataset[i][2])
                else:
                    data.append(self.dataset[i]["rho"])
        else:
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][4][1])
                elif keyword == "bulk":
                    data.append(self.dataset[i][4][1][0])
                elif keyword == "solid":
                    data.append(self.dataset[i][4][1][1])
                elif keyword == "fluid":
                    data.append(self.dataset[i][4][1][2])
        #
        return np.array(data)
    #
    def extract_elastic_moduli(self, type="sequence", keyword="bulk", dict=False):
        """Returns a list that contains the elastic moduli of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of elastic moduli
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    if keyword == "all":
                        data.append(self.dataset[i][j][4][2])
                    elif keyword == "bulk" or keyword == "K":
                        data.append(self.dataset[i][j][4][2][0])
                    elif keyword == "shear" or keyword == "G":
                        data.append(self.dataset[i][j][4][2][1])
                    elif keyword == "young" or keyword == "E":
                        data.append(self.dataset[i][j][4][2][2])
                    elif keyword == "poisson" or keyword == "nu":
                        data.append(self.dataset[i][j][4][2][3])
        elif type == "random":
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][2])
                elif keyword == "bulk" or keyword == "K":
                    data.append(self.dataset[i][2][0])
                elif keyword == "shear" or keyword == "G":
                    data.append(self.dataset[i][2][1])
                elif keyword == "young" or keyword == "E":
                    data.append(self.dataset[i][2][2])
                elif keyword == "poisson" or keyword == "nu":
                    data.append(self.dataset[i][2][3])
        elif type == "mineral":
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][3])
                elif keyword == "bulk" or keyword == "K":
                    data.append(self.dataset[i][3][0])
                elif keyword == "shear" or keyword == "G":
                    data.append(self.dataset[i][3][1])
                elif keyword == "young" or keyword == "E":
                    data.append(self.dataset[i][3][2])
                elif keyword == "poisson" or keyword == "nu":
                    data.append(self.dataset[i][3][3])
        elif type == "rock":
            for i in range(len(self.dataset)):
                if dict == False:
                    if keyword == "all":
                        data.append(self.dataset[i][2])
                    elif keyword == "bulk" or keyword == "K":
                        data.append(self.dataset[i][2][0])
                    elif keyword == "shear" or keyword == "G":
                        data.append(self.dataset[i][2][1])
                    elif keyword == "young" or keyword == "E":
                        data.append(self.dataset[i][2][2])
                    elif keyword == "poisson" or keyword == "nu":
                        data.append(self.dataset[i][2][3])
                else:
                    data.append(self.dataset[i][keyword])
        else:
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][4][2])
                elif keyword == "bulk" or keyword == "K":
                    data.append(self.dataset[i][4][2][0])
                elif keyword == "shear" or keyword == "G":
                    data.append(self.dataset[i][4][2][1])
                elif keyword == "young" or keyword == "E":
                    data.append(self.dataset[i][4][2][2])
                elif keyword == "poisson" or keyword == "nu":
                    data.append(self.dataset[i][4][2][3])
        #
        return np.array(data)
    #
    def extract_seismic_velocities(self, type="sequence", keyword="p-wave", dict=False):
        """Returns a list that contains the seismic velocities of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of seismic velocities
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    if keyword == "all":
                        data.append(self.dataset[i][j][4][3])
                    elif keyword in ["p-wave", "compressional", "vP"]:
                        data.append(self.dataset[i][j][4][3][0])
                    elif keyword in ["s-wave", "shear", "vS"]:
                        data.append(self.dataset[i][j][4][3][1])
        elif type == "random":
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][3])
                elif keyword in ["p-wave", "compressional", "vP"]:
                    data.append(self.dataset[i][3][0])
                elif keyword in ["s-wave", "shear", "vS"]:
                    data.append(self.dataset[i][3][1])
        elif type == "mineral":
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][4])
                elif keyword in ["p-wave", "compressional", "vP"]:
                    data.append(self.dataset[i][4][0])
                elif keyword in ["s-wave", "shear", "vS"]:
                    data.append(self.dataset[i][4][1])
                elif keyword in ["vPvS"]:
                    data.append(self.dataset[i][4][2])
        elif type == "rock":
            for i in range(len(self.dataset)):
                if dict == False:
                    if keyword == "all":
                        data.append(self.dataset[i][3])
                    elif keyword in ["p-wave", "compressional", "vP"]:
                        data.append(self.dataset[i][3][0])
                    elif keyword in ["s-wave", "shear", "vS"]:
                        data.append(self.dataset[i][3][1])
                else:
                    data.append(self.dataset[i][keyword])
        else:
            for i in range(len(self.dataset)):
                if keyword == "all":
                    data.append(self.dataset[i][4][3])
                elif keyword in ["p-wave", "compressional", "vP"]:
                    data.append(self.dataset[i][4][3][0])
                elif keyword in ["s-wave", "shear", "vS"]:
                    data.append(self.dataset[i][4][3][1])
        #
        return np.array(data)
    #
    def extract_porosity(self, type="sequence", keyword="phi", dict=False):
        """Returns a list that contains the porosities of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of porosities
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][4][0])
        elif type == "random":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][4][0])
        elif type == "rock":
            for i in range(len(self.dataset)):
                if dict == False:
                    data.append(self.dataset[i][4][0])
                else:
                    data.append(self.dataset[i][keyword])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][4][4][0])
        #
        return np.array(data)
    #
    def extract_fluid(self, type="sequence"):
        """Returns a list that contains the fluid within the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of fluids
        """
        data = []
        if type == "sequence":
            data.append(self.dataset[0][0][4][5])
        else:
            data.append(self.dataset[0][4][5])
        #
        return data
    #
    def extract_gamma_ray(self, type="sequence", keyword="GR", dict=False):
        """Returns a list that contains the natural gamma ray values of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of natural gamma ray values
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][6][0])
        elif type == "random":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][6][0])
        elif type == "mineral":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][5][0])
        elif type == "rock":
            for i in range(len(self.dataset)):
                if dict == False:
                    data.append(self.dataset[i][6][0])
                else:
                    data.append(self.dataset[i][keyword])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][4][6][0])
        #
        return np.array(data)
    #
    def extract_photoelectricity(self, type="sequence", keyword="PE", dict=False):
        """Returns a list that contains the photoelectricity values of the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of photoelectricity values
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][6][1])
        elif type == "random":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][6][1])
        elif type == "mineral":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][5][2])
        elif type == "rock":
            for i in range(len(self.dataset)):
                if dict == False:
                    data.append(self.dataset[i][6][1])
                else:
                    data.append(self.dataset[i][keyword])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][4][6][1])
        #
        return np.array(data)
    #
    def extract_data(self, keyword="rho"):
        """Returns a list that contains the photoelectricity values of the previously generated rock units.
        **Arguments**:
            type: keyword
        **Outputs**:
            data: list of photoelectricity values
        """
        data = []
        for index, item in enumerate(self.dataset, start=0):
            data.append(item[keyword])
        #
        return np.array(data)
    #
    def extract_element_amounts(self, type="sequence", pos=None, element=None):
        """Returns a list that contains the amounts of the elements within the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of element amounts
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][7])
        elif type == "mineral":
            for i in range(len(self.dataset)):
                if pos != None:
                    data.append(self.dataset[i][6][pos])
                if element != None:
                    for item in self.dataset[i][6]:
                        if element in item:
                            data.append(item[2])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][4][7])
        #
        return np.array(data)
    #
    def extract_mineral_amounts(self, type="sequence"):
        """Returns a list that contains the amounts of the minerals within the previously generated rock units.
        **Arguments**:
            type: sequence, rock
        **Outputs**:
            data: list of mineral amounts
        """
        data = []
        if type == "sequence":
            for i in range(len(self.dataset)):
                for j in range(len(self.dataset[i])):
                    data.append(self.dataset[i][j][4][8])
        elif type == "rock":
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][0][1])
        else:
            for i in range(len(self.dataset)):
                data.append(self.dataset[i][4][8])
        #
        return np.array(data)

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
            sequence.extend(["water", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR, 1), round(phiN, 1), "water", composition, round(poisson,2), round(PE,2)])
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
            sequence.extend(["dry sand", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR, 1), round(phiN, 1), "water", composition, round(poisson,2), round(PE,2)])
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
            sequence.extend(["soil", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR, 1), round(phiN, 1), "water", composition, round(poisson,2), round(PE,2)])
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
        sequence.extend(["wet sand", thicknessUnit, self.actualThickness, newThickness, round(rho,3), velocities, round(GR, 1), round(phiN, 1), "water", composition, round(poisson,2), round(PE,2)])
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
                    sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                    sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                    sequence.append(["sandstone", round(d2, 1), round(top2, 1), round(bottom2, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                    sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                    sequence.append(["sandstone", round(d2, 1), round(top2, 1), round(bottom2, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                    sequence.append(["sandstone", round(d3, 1), round(top3, 1), round(bottom3, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2), dataSandstone[2]])
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
                        sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["limestone", round(d2, 1), round(top2, 1), round(bottom2, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["limestone", round(d2, 1), round(top2, 1), round(bottom2, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["limestone", round(d3, 1), round(top3, 1), round(bottom3, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                    sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["dolomite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["dolomite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["dolomite", round(d2, 1), round(top2, 1), round(bottom2, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["dolomite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["dolomite", round(d2, 1), round(top2, 1), round(bottom2, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                        sequence.append(["dolomite", round(d3, 1), round(top3, 1), round(bottom3, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                    sequence.append(["dolomite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
            sequence.append(["shale", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                sequence.append(["halite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                sequence.append(["anhydrite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                    sequence.append(["granite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                    sequence.append(["basalt", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                sequence.append(["granite", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
                sequence.append(["basalt", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
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
            sequence.append(["ore", round(d, 1), round(top, 1), round(bottom, 1), rho, velocities, GR, phiN, fluid, composition, round(poisson,2), round(PE,2)])
        #
        return sequence
#
class SedimentaryBasin:
    #
    def __init__(self, actualThickness=0, parts=None, maximum_thickness=None):
        self.actualThickness = actualThickness
        self.parts = parts
        self.maximum_thickness = maximum_thickness
    #
    def create_soil(self, thickness=None, grainsize=False):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(5, 10)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Soil()
                data_soil = data.create_simple_soil(grainsize_list=grainsize)
                top = self.actualThickness
                bottom = top + d
                sequence.append(["soil", round(d, 1), round(top, 1), round(bottom, 1), data_soil])
            else:
                data = Soil()
                data_soil = data.create_simple_soil(amounts=sequence[-1][4][-2], grainsize_list=grainsize)
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["soil", round(d, 1), round(top, 1), round(bottom, 1), data_soil])
        #
        return sequence
    #
    def create_sand(self, thickness=None, grainsize=False):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(5, 15)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Soil()
                data_sand = data.create_simple_sand(grainsize_list=grainsize)
                top = self.actualThickness
                bottom = top + d
                sequence.append(["sand", round(d, 1), round(top, 1), round(bottom, 1), data_sand])
            else:
                data = Soil()
                data_sand = data.create_simple_sand(amounts=sequence[-1][4][-2], grainsize_list=grainsize)
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["sand", round(d, 1), round(top, 1), round(bottom, 1), data_sand])
        #
        return sequence
    #
    def create_sandstone(self, thickness=None, fluid="water", phi=None, Qz_rich=False, keyword=None):
        self.fluid = fluid
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(5, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0 and self.fluid == "water":
                data = sandstone("water", actualThickness=self.actualThickness)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(porosity=phi)
                top = self.actualThickness
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i == 0 and self.fluid == "gas":
                data = sandstone("gas", actualThickness=self.actualThickness)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(porosity=phi)
                top = self.actualThickness
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i == 0 and self.fluid == "oil":
                data = sandstone("oil", actualThickness=self.actualThickness)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(porosity=phi)
                top = self.actualThickness
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i == 0 and self.fluid == "air":
                data = sandstone("air", actualThickness=self.actualThickness)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(porosity=phi)
                top = self.actualThickness
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i > 0 and self.fluid == "water":
                data = sandstone("water", actualThickness=bottom)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(amounts=sequence[-1][-1][8], porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(amounts=sequence[-1][-1][8], porosity=phi)
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i > 0 and self.fluid == "gas":
                data = sandstone("gas", actualThickness=bottom)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(amounts=sequence[-1][-1][8], porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(amounts=sequence[-1][-1][8], porosity=phi)
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i > 0 and self.fluid == "oil":
                data = sandstone("oil", actualThickness=bottom)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(amounts=sequence[-1][-1][8], porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(amounts=sequence[-1][-1][8], porosity=phi)
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
            elif i > 0 and self.fluid == "air":
                data = sandstone("air", actualThickness=bottom)
                if keyword == None:
                    data_sandstone = data.create_simple_sandstone(amounts=sequence[-1][-1][8], porosity=phi, pure=Qz_rich)
                else:
                    data_sandstone = data.create_feldspathic_sandstone(amounts=sequence[-1][-1][8], porosity=phi)
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["sandstone", round(d, 1), round(top, 1), round(bottom, 1), data_sandstone])
        #
        return sequence
    #
    def create_limestone(self, thickness=None, fluid="water"):
        self.fluid = fluid
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(5, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0 and self.fluid == "water":
                data = limestone("water", actualThickness=self.actualThickness)
                data_limestone = data.create_simple_limestone()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), data_limestone])
            elif i == 0 and self.fluid == "gas":
                data = limestone("gas", actualThickness=self.actualThickness)
                data_limestone = data.create_simple_limestone()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), data_limestone])
            elif i == 0 and self.fluid == "oil":
                data = limestone("oil", actualThickness=self.actualThickness)
                data_limestone = data.create_simple_limestone()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), data_limestone])
            elif i > 0 and self.fluid == "water":
                data = limestone("water", actualThickness=bottom)
                data_limestone = data.create_simple_limestone(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), data_limestone])
            elif i > 0 and self.fluid == "gas":
                data = limestone("gas", actualThickness=bottom)
                data_limestone = data.create_simple_limestone(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), data_limestone])
            elif i > 0 and self.fluid == "oil":
                data = limestone("oil", actualThickness=bottom)
                data_limestone = data.create_simple_limestone(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["limestone", round(d, 1), round(top, 1), round(bottom, 1), data_limestone])
        #
        return sequence
    #
    def create_shale(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(5, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = shale()
                data_shale = data.create_shale()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["shale", round(d, 1), round(top, 1), round(bottom, 1), data_shale])
            else:
                data = shale()
                data_shale = data.create_shale(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["shale", round(d, 1), round(top, 1), round(bottom, 1), data_shale])
        #
        return sequence
    #
    def create_shale_complex(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(5, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = shale()
                data_shale = data.create_simple_shale()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["shale", round(d, 1), round(top, 1), round(bottom, 1), data_shale])
            else:
                data = shale()
                data_shale = data.create_simple_shale(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["shale", round(d, 1), round(top, 1), round(bottom, 1), data_shale])
        #
        return sequence
    #
    def create_rocksalt(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Evaporites("water", actualThickness=self.actualThickness)
                data_halite = data.create_simple_rocksalt()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["rock salt", round(d, 1), round(top, 1), round(bottom, 1), data_halite])
            else:
                data = Evaporites("water", actualThickness=self.actualThickness)
                data_halite = data.create_simple_rocksalt(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["rock salt", round(d, 1), round(top, 1), round(bottom, 1), data_halite])
        #
        return sequence
    #
    def create_granite(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_granite = data.create_simple_granite()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["granite", round(d, 1), round(top, 1), round(bottom, 1), data_granite])
            else:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_granite = data.create_simple_granite(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["granite", round(d, 1), round(top, 1), round(bottom, 1), data_granite])
        #
        return sequence
    #
    def create_basalt(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Volcanic("water", actualThickness=self.actualThickness)
                data_basalt = data.create_simple_basalt()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["basalt", round(d, 1), round(top, 1), round(bottom, 1), data_basalt])
            else:
                data = Volcanic("water", actualThickness=self.actualThickness)
                data_basalt = data.create_simple_basalt(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["basalt", round(d, 1), round(top, 1), round(bottom, 1), data_basalt])
        #
        return sequence
    #
    def create_sedimentary_basin(self, maximum_thickness=500, n_units = 15, csv_stratigraphy=False, csv_lithology=False,
                                 excludeRocksalt=False, excludeLimestone=False):
        self.maximum_thickness = maximum_thickness
        sequence = []
        thickness_list = [0, rd.randint(5, 10), rd.randint(5, 15)]
        depth_list = [thickness_list[1]]
        #
        for i in range(2, n_units):
            if len(thickness_list) < n_units:
                thickness_list.append(rd.randint(int(0.75*(maximum_thickness-sum(thickness_list))/(n_units-len(thickness_list)+1)), int(1.25*(maximum_thickness-sum(thickness_list))/(n_units-len(thickness_list)+1))))
                depth_list.append(sum(thickness_list[:i+1]))
            else:
                thickness_list.append(maximum_thickness-sum(thickness_list))
                depth_list.append(sum(thickness_list[:i+1]))
        depth_list.append(sum(thickness_list))
        del thickness_list[0]
        #
        condition = False
        while condition == False:
            for i in range(n_units):
                if len(sequence) == 0:
                    data = SedimentaryBasin(parts=self.parts)
                    data_soil = data.create_soil(thickness=thickness_list[i])
                    sequence.append(data_soil)
                    actual_depth = depth_list[i]
                elif len(sequence) == 1:
                    data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                    data_sand = data.create_sand(thickness=thickness_list[i])
                    sequence.append(data_sand)
                    actual_depth = depth_list[i]
                elif actual_depth < 0.85*self.maximum_thickness and sequence[-1][-1][0] not in ["shale", "rock salt"] and sequence[-1][-1][-1][5] not in ["gas", "oil"]:
                    condition_2 = False
                    while condition_2 == False:
                        magicnumber = rd.random()
                        if magicnumber < 0.4:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_sandstone = data.create_sandstone(thickness=thickness_list[i])
                            sequence.append(data_sandstone)
                            actual_depth = depth_list[i]
                            condition_2 = True
                        elif 0.4 <= magicnumber < 0.8:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_shale = data.create_shale(thickness=thickness_list[i])
                            sequence.append(data_shale)
                            actual_depth = depth_list[i]
                            condition_2 = True
                        elif 0.8 <= magicnumber < 0.95 and excludeLimestone == False:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_limestone = data.create_limestone(thickness=thickness_list[i])
                            sequence.append(data_limestone)
                            actual_depth = depth_list[i]
                            condition_2 = True
                        elif magicnumber >= 0.95 and excludeRocksalt == False:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_rocksalt = data.create_rocksalt(thickness=thickness_list[i])
                            sequence.append(data_rocksalt)
                            actual_depth = depth_list[i]
                            condition_2 = True
                        elif magicnumber >= 0.95 and excludeRocksalt == True:
                            condition_2 = False
                elif actual_depth < 0.85*self.maximum_thickness and sequence[-1][-1][0] in ["shale", "rock salt"]:
                    condition_3 = False
                    while condition_3 == False:
                        magicnumber = rd.random()
                        if magicnumber < 0.1 and excludeLimestone == False:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_limestone = data.create_limestone(thickness=thickness_list[i], fluid="gas")
                            sequence.append(data_limestone)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif 0.1 <= magicnumber < 0.2 and excludeLimestone == False:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_limestone = data.create_limestone(thickness=thickness_list[i], fluid="oil")
                            sequence.append(data_limestone)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif 0.2 <= magicnumber < 0.3:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_sandstone = data.create_sandstone(thickness=thickness_list[i], fluid="gas")
                            sequence.append(data_sandstone)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif 0.3 <= magicnumber < 0.4:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_sandstone = data.create_sandstone(thickness=thickness_list[i], fluid="oil")
                            sequence.append(data_sandstone)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif 0.4 <= magicnumber < 0.6:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_sandstone = data.create_sandstone(thickness=thickness_list[i], fluid="water")
                            sequence.append(data_sandstone)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif 0.6 <= magicnumber < 0.8 and excludeLimestone == False:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_limestone = data.create_limestone(thickness=thickness_list[i], fluid="water")
                            sequence.append(data_limestone)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif magicnumber >= 0.8:
                            data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                            data_shale = data.create_shale(thickness=thickness_list[i])
                            sequence.append(data_shale)
                            actual_depth = depth_list[i]
                            condition_3 = True
                        elif excludeLimestone == True:
                            condition_3 = False
                elif actual_depth < 0.85*self.maximum_thickness and sequence[-1][-1][4][5] == "gas" and sequence[-1][-1][0] == "sandstone":
                    data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                    data_sandstone = data.create_sandstone(thickness=thickness_list[i], fluid="oil")
                    sequence.append(data_sandstone)
                    actual_depth = depth_list[i]
                elif actual_depth < 0.85*self.maximum_thickness and sequence[-1][-1][4][5] == "oil" and sequence[-1][-1][0] == "sandstone":
                    data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                    data_sandstone = data.create_sandstone(thickness=thickness_list[i], fluid="water")
                    sequence.append(data_sandstone)
                    actual_depth = depth_list[i]
                elif actual_depth < 0.85*self.maximum_thickness and sequence[-1][-1][4][5] == "gas" and sequence[-1][-1][0] == "limestone" and excludeLimestone == False:
                    data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                    data_limestone = data.create_limestone(thickness=thickness_list[i], fluid="oil")
                    sequence.append(data_limestone)
                    actual_depth = depth_list[i]
                elif actual_depth < 0.85*self.maximum_thickness and sequence[-1][-1][4][5] == "oil" and sequence[-1][-1][0] == "limestone" and excludeLimestone == False:
                    data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                    data_limestone = data.create_limestone(thickness=thickness_list[i], fluid="water")
                    sequence.append(data_limestone)
                    actual_depth = depth_list[i]
                elif actual_depth >= 0.85*self.maximum_thickness and actual_depth < maximum_thickness:
                    if sequence[-1][-1][0] not in ["granite", "basalt"]:
                        magicnumber = rd.randint(0, 1)
                    elif sequence[-1][-1][0] == "granite":
                        magicnumber = 0
                    elif sequence[-1][-1][0] == "basalt":
                        magicnumber = 1
                    if magicnumber == 0 or sequence[-1][-1][0] == "granite":
                        data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                        data_granite = data.create_granite(thickness=thickness_list[i])
                        sequence.append(data_granite)
                        actual_depth = depth_list[i]
                    elif magicnumber == 1 or sequence[-1][-1][0] == "basalt":
                        data = SedimentaryBasin(parts=self.parts, actualThickness=actual_depth)
                        data_basalt = data.create_basalt(thickness=thickness_list[i])
                        sequence.append(data_basalt)
                        actual_depth = depth_list[i]
                elif actual_depth >= self.maximum_thickness:
                    condition = True
                #print("Actual depth for", len(sequence), "units:", actual_depth, "m. -->", sequence[-1][-1][0], sequence[-1][-1][4][5])
        #
        units_list = []
        for i in range(len(sequence)):
            if sequence[i][0][0] not in units_list:
                units_list.append(sequence[i][0][0])
        #
        if csv_stratigraphy == True:
            try:
                file_sb = open("Data_SB_Stratigraphy.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            file_sb.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("POISSON") + "," + str("FLUID") + "\n")
            for i in range(len(sequence)):
                for j in range(len(sequence[i])):
                    if sequence[i][j][0] == "soil":
                        file_sb.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5][0]) + "\n")
                    else:
                        file_sb.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "\n")
        #
        if csv_lithology == True:
            try:
                file_soil = open("Data_SB_Soil.csv", "w")
                file_sand = open("Data_SB_Sand.csv", "w")
                file_sandstone = open("Data_SB_Sandstone.csv", "w")
                file_limestone = open("Data_SB_Limestone.csv", "w")
                file_shale = open("Data_SB_Shale.csv", "w")
                #file_shale_complex = open("Data_SB_Shale_complex.csv", "w")
                if "rock salt" in units_list:
                    file_rocksalt = open("Data_SB_Rocksalt.csv", "w")
                    file_rocksalt.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("HL") + "," + str("ANH") + "," + str("GP") + "," + str("SYL") + "\n")
                if "granite" in units_list:
                    file_granite = open("Data_SB_Granite.csv", "w")
                    file_granite.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("QZ") + "," + str("KFS") + "," + str("PL") + "," + str("BT") + "," + str("MS") + "," + str("ACT") + "," + str("TR") + "," + str("AUG") + "\n")
                if "basalt" in units_list:
                    file_basalt = open("Data_SB_Basalt.csv", "w")
                    file_basalt.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("QZ") + "," + str("KFS") + "," + str("PL") + "," + str("BT") + "," + str("MS") + "," + str("ACT") + "," + str("TR") + "," + str("AUG") + "\n")
            except:
                print("Error")
                sys.exit(0)
            #
            file_soil.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL1") + "," + str("RHOFL2") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("QZ") + "," + str("ILT") + "," + str("KLN") + "," + str("ORG") + "\n")
            file_sand.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("QZ") + "," + str("ILT") + "," + str("KLN") + "," + str("ORG") + "\n")
            file_sandstone.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("QZ") + "," + str("KFS") + "," + str("PL") + "," + str("CAL") + "," + str("CHL") + "," + str("MS") + "," + str("HEM") + "\n")
            file_limestone.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("CAL") + "," + str("ARG") + "," + str("DOL") + "," + str("SD") + "," + str("QZ") + "," + str("KFS") + "," + str("PL") + "," + str("MNT") + "," + str("KLN") + "," + str("CHL") + "," + str("PY") + "\n")
            file_shale.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("ORG") + "," + str("QZ") + "," + str("CAL") + "," + str("ILT") + "," + str("KLN") + "," + str("MNT") + "\n")
            #file_shale_complex.write(str("LITHOLOGY") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("RHOMA") + "," + str("RHOFL") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PE") + "," + str("PHIN") + "," + str("KMOD") + "," + str("GMOD") + "," + str("EMOD") + "," + str("POISSON") + "," + str("FLUID") + "," + str("ORG") + "," + str("QZ") + "," + str("CAL") + "," + str("PY") + "," + str("ILT") + "," + str("KLN") + "," + str("MNT") + "," + str("BT") + "," + str("MS") + "," + str("URN") + "\n")
            for i in range(len(sequence)):
                for j in range(len(sequence[i])):
                    if sequence[i][j][0] == "soil":
                        file_soil.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][1][3]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5][0]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "\n")
                    elif sequence[i][j][0] == "sand":
                        file_sand.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "\n")
                    elif sequence[i][j][0] == "sandstone":
                        file_sandstone.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "," + str(sequence[i][j][4][8][4]) + "," + str(sequence[i][j][4][8][5]) + "," + str(sequence[i][j][4][8][6]) + "\n")
                    elif sequence[i][j][0] == "limestone":
                        file_limestone.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "," + str(sequence[i][j][4][8][4]) + "," + str(sequence[i][j][4][8][5]) + "," + str(sequence[i][j][4][8][6]) + "," + str(sequence[i][j][4][8][7]) + "," + str(sequence[i][j][4][8][8]) + "," + str(sequence[i][j][4][8][9]) + "," + str(sequence[i][j][4][8][10]) + "\n")
                    elif sequence[i][j][0] == "shale":
                        file_shale.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + str(sequence[i][j][4][8][3]) + "," + str(sequence[i][j][4][8][4]) + "," + str(sequence[i][j][4][8][5]) + "\n")
                        #file_shale_complex.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "," + str(sequence[i][j][4][8][4]) + "," + str(sequence[i][j][4][8][5]) + "," + str(sequence[i][j][4][8][6]) + "," + str(sequence[i][j][4][8][7]) + "," + str(sequence[i][j][4][8][8]) + "," + str(sequence[i][j][4][8][9]) + "\n")
                    elif sequence[i][j][0] == "rock salt":
                        file_rocksalt.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "\n")
                    elif sequence[i][j][0] == "granite":
                        file_granite.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "," + str(sequence[i][j][4][8][4]) + "," + str(sequence[i][j][4][8][5]) + "," + str(sequence[i][j][4][8][6]) + "," + str(sequence[i][j][4][8][7]) + "\n")
                    elif sequence[i][j][0] == "basalt":
                        file_basalt.write(str(sequence[i][j][0]) + "," + str(sequence[i][j][1]) + "," + str(sequence[i][j][2]) + "," + str(sequence[i][j][3]) + "," + str(sequence[i][j][4][1][0]) + "," + str(sequence[i][j][4][1][1]) + "," + str(sequence[i][j][4][1][2]) + "," + str(sequence[i][j][4][3][0]) + "," + str(sequence[i][j][4][3][1]) + "," + str(sequence[i][j][4][6][0]) + "," + str(sequence[i][j][4][6][1]) + "," + str(sequence[i][j][4][4][0]) + "," + str(sequence[i][j][4][2][0]) + "," + str(sequence[i][j][4][2][1]) + "," + str(sequence[i][j][4][2][2]) + "," + str(sequence[i][j][4][2][3]) + "," + str(sequence[i][j][4][5]) + "," + str(sequence[i][j][4][8][0]) + "," + str(sequence[i][j][4][8][1]) + "," + str(sequence[i][j][4][8][2]) + "," + str(sequence[i][j][4][8][3]) + "," + str(sequence[i][j][4][8][4]) + "," + str(sequence[i][j][4][8][5]) + "," + str(sequence[i][j][4][8][6]) + "," + str(sequence[i][j][4][8][7]) + "\n")
        #
        return sequence
#
class Plutonite:
    #
    def __init__(self, actualThickness=0, parts=None, maximum_thickness=None):
        self.actualThickness = actualThickness
        self.parts = parts
        self.maximum_thickness = maximum_thickness
    #
    def create_gabbro(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_simple_gabbro()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["gabbro", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
            else:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_simple_gabbro(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["gabbro", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
        #
        return sequence
    #
    def create_granite(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_simple_granite()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["granite", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
            else:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_simple_granite(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["granite", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
        #
        return sequence
    #
    def create_felsic(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_felsic()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["felsic rock", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
            else:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_felsic(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["felsic rock", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
        #
        return sequence
    #
    def create_intermediate(self, thickness=None):
        sequence = []
        if thickness == None:
            thickness_max = rd.randint(10, 20)
        else:
            thickness_max = thickness
        #
        if self.parts == None:
            d = 1.0
            self.parts = thickness_max
        else:
            d = float(thickness_max/self.parts)
        #
        for i in range(self.parts):
            if i == 0:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_intermediate()
                top = self.actualThickness
                bottom = top + d
                sequence.append(["intermediate rock", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
            else:
                data = Plutonic("water", actualThickness=self.actualThickness)
                data_rock = data.create_intermediate(amounts=sequence[-1][-1][8])
                top = self.actualThickness + i*d
                bottom = top + d
                sequence.append(["intermediate rock", round(d, 1), round(top, 1), round(bottom, 1), data_rock])
        #
        return sequence