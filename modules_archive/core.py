#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		core.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		13.01.2021

#-----------------------------------------------

## MODULES
import numpy as np
from scipy import constants
from numpy import round
import sys
from random import *
from modules import minerals
from modules.elements import elements
import matplotlib.pyplot as plt
from modules.geochemistry import elementanalysis
from modules import siliciclastics, carbonates, igneous, evaporites
from modules import fluids

## CONSTANTS / PROPERTIES
pi = constants.pi
avogadro = constants.Avogadro

class sample:
    #
    def __init__(self, sequences, radius):
        self.sequences = sequences
        self.radius = radius
    #
    def createCore(self):
        volume = []
        for i in range(0, len(self.sequences)):
            volume.append(pi*self.radius**2*self.sequences[i][1])
        #
        mass = []
        for i in range(0, len(volume)):
            mass.append(volume[i]*self.sequences[i][4]*1000)
        #
        data = []
        for i in range(0, len(volume)):
            data.append([self.sequences[i][0], volume[i], mass[i], self.sequences[i][4], self.sequences[i][9]])
            #
            # WRITING A CSV FILE
            #
            try:
                SequencyStratigraphy = open("./outputs/SequenceStratigraphy.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            SequencyStratigraphy.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("PE") + "\n")
            for i in range(0, len(self.sequences)):
                SequencyStratigraphy.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][11]) + "\n")
            #
            try:
                RockProperties_Sandstone = open("./outputs/RockProperties_Sandstone.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Sandstone.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Qz") + "," + str("Cal") + "," + str("Dol") + "," + str("Kfs") + "," + str("Pl") + "," + str("Bt") + "," + str("Glt") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("w(Qz)") + "," + str("w(Cal)") + "," + str("w(Dol)") + "," + str("w(Kfs)") + "," + str("w(Pl)") + "," + str("w(Bt)") + "," + str("w(Glt)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "sandstone":
                    VCore = pi*self.radius**2*self.sequences[i][1]
                    mCore = pi*self.radius**2*self.sequences[i][1]*self.sequences[i][4]*1000
                    RockProperties_Sandstone.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "\n")
            #
            try:
                RockProperties_Limestone = open("./outputs/RockProperties_Limestone.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Limestone.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Cal") + "," + str("Dol") + "," + str("Sd") + "," + str("Qz") + "," + str("Kfs") + "," + str("Pl") + "," + str("Kln") + "," + str("Chl") + "," + str("Ilt") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Cal)") + "," + str("m(Dol)") + "," + str("m(Sd)") + "," + str("m(Qz)") + "," + str("m(Kfs)") + "," + str("m(Pl)") + "," + str("m(Kln)") + "," + str("m(Chl)") + "," + str("m(Ilt)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "limestone":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Limestone.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.sequences[i][9][7][1]) + "," + str(self.sequences[i][9][8][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "," + str(mCore*self.sequences[i][9][6][1]) + "," + str(mCore*self.sequences[i][9][7][1]) + "," + str(mCore*self.sequences[i][9][8][1]) + "\n")
            #
            try:
                RockProperties_Dolomite = open("./outputs/RockProperties_Dolomite.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Dolomite.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Cal") + "," + str("Dol") + "," + str("Mgs") + "," + str("Qz") + "," + str("Kfs") + "," + str("Pl") + "," + str("Kln") + "," + str("Chl") + "," + str("Ilt") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Cal)") + "," + str("m(Dol)") + "," + str("m(Mgs)") + "," + str("m(Qz)") + "," + str("m(Kfs)") + "," + str("m(Pl)") + "," + str("m(Kln)") + "," + str("m(Chl)") + "," + str("m(Ilt)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "dolomite":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Dolomite.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.sequences[i][9][7][1]) + "," + str(self.sequences[i][9][8][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "," + str(mCore*self.sequences[i][9][6][1]) + "," + str(mCore*self.sequences[i][9][7][1]) + "," + str(mCore*self.sequences[i][9][8][1]) + "\n")
            #
            try:
                RockProperties_Shale = open("./outputs/RockProperties_Shale.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Shale.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Qz") + "," + str("Kln") + "," + str("Chl") + "," + str("Ilt") + "," + str("Cal") + "," + str("Dol") + "," + str("Kfs") + "," + str("Pl") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Qz)") + "," + str("m(Kln)") + "," + str("m(Chl)") + "," + str("m(Ilt)") + "," + str("m(Cal)") + "," + str("m(Dol)") + "," + str("m(Kfs)") + "," + str("m(Pl)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "shale":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Shale.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.sequences[i][9][7][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "," + str(mCore*self.sequences[i][9][6][1]) + "," + str(mCore*self.sequences[i][9][7][1]) + "\n")
            #
            try:
                RockProperties_Halite = open("./outputs/RockProperties_Halite.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Halite.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Hl") + "," + str("Anh") + "," + str("Gp") + "," + str("Syl") + "," + str("Ilt") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Hl)") + "," + str("m(Anh)") + "," + str("m(Gp)") + "," + str("m(Syl)") + "," + str("m(Ilt)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "halite":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Halite.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "\n")
            #
            try:
                RockProperties_Anhydrite = open("./outputs/RockProperties_Anhydrite.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Anhydrite.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Anh") + "," + str("Hl") + "," + str("Cal") + "," + str("Gn") + "," + str("Ccp") + "," + str("Mol") + "," + str("Py") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Anh)") + "," + str("m(Hl)") + "," + str("m(Cal)") + "," + str("m(Gn)") + "," + str("m(Ccp)") + "," + str("m(Mol)") + "," + str("m(Py)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "anhydrite":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Anhydrite.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "," + str(mCore*self.sequences[i][9][6][1]) + "\n")
            #
            try:
                RockProperties_Granite = open("./outputs/RockProperties_Granite.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Granite.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Qz") + "," + str("Or") + "," + str("Ab") + "," + str("An") + "," + str("Bt") + "," + str("Ms") + "," + str("Py") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Qz)") + "," + str("m(Or)") + "," + str("m(Ab)") + "," + str("m(An)") + "," + str("m(Bt)") + "," + str("m(Ms)") + "," + str("m(Py)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "granite":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Granite.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "," + str(mCore*self.sequences[i][9][6][1]) + "\n")
            #
            try:
                RockProperties_Basalt = open("./outputs/RockProperties_Basalt.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Basalt.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Pl") + "," + str("En") + "," + str("Fs") + "," + str("Bt") + "," + str("Act") + "," + str("Tr") + "," + str("Ol") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Pl)") + "," + str("m(En)") + "," + str("m(Fs)") + "," + str("m(Bt)") + "," + str("m(Act)") + "," + str("m(Tr)") + "," + str("m(Ol)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "basalt":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Basalt.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.sequences[i][9][6][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "," + str(mCore*self.sequences[i][9][6][1]) + "\n")
            #
            try:
                RockProperties_Ore_Fe = open("./outputs/RockProperties_Ore_Fe.csv", "w")
            except:
                print("Error")
                sys.exit(0)
            RockProperties_Ore_Fe.write(str("ROCK") + "," + str("THICKNESS") + "," + str("TOP") + "," + str("BOTTOM") + "," + str("RHOB") + "," + str("VP") + "," + str("VS") + "," + str("GR") + "," + str("PHIN") + "," + str("FLUID") + "," + str("POISSON") + "," + str("Mag") + "," + str("Hem") + "," + str("Qz") + "," + str("Cal") + "," + str("Kfs") + "," + str("Pl") + "," + str("R(Core)") + "," + str("V(Core)") + "," + str("m(Core)") + "," + str("m(Mag)") + "," + str("m(Hem)") + "," + str("m(Qz)") + "," + str("m(Cal)") + "," + str("m(Kfs)") + "," + str("m(Pl)") + "\n")
            for i in range(0, len(self.sequences)):
                if self.sequences[i][0] == "ore":
                    VCore = pi * self.radius ** 2 * self.sequences[i][1]
                    mCore = pi * self.radius ** 2 * self.sequences[i][1] * self.sequences[i][4] * 1000
                    RockProperties_Ore_Fe.write(str(self.sequences[i][0]) + "," + str(self.sequences[i][1]) + "," + str(self.sequences[i][2]) + "," + str(self.sequences[i][3]) + "," + str(self.sequences[i][4]) + "," + str(self.sequences[i][5][0]) + "," + str(self.sequences[i][5][1]) + "," + str(self.sequences[i][6]) + "," + str(self.sequences[i][7]) + "," + str(self.sequences[i][8]) + "," + str(self.sequences[i][10]) + "," + str(self.sequences[i][9][0][1]) + "," + str(self.sequences[i][9][1][1]) + "," + str(self.sequences[i][9][2][1]) + "," + str(self.sequences[i][9][3][1]) + "," + str(self.sequences[i][9][4][1]) + "," + str(self.sequences[i][9][5][1]) + "," + str(self.radius) + "," + str(VCore) + "," + str(mCore) + "," + str(mCore*self.sequences[i][9][0][1]) + "," + str(mCore*self.sequences[i][9][1][1]) + "," + str(mCore*self.sequences[i][9][2][1]) + "," + str(mCore*self.sequences[i][9][3][1]) + "," + str(mCore*self.sequences[i][9][4][1]) + "," + str(mCore*self.sequences[i][9][5][1]) + "\n")
            #
        return data
    #
class geochemistry:
    #
    def __init__(self, input, sequences):
        self.input = input
        self.sequences = sequences
    #
    def calculateConcentration(self):
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemH = elements.H(self)
        chemO = elements.O(self)
        chemSi = elements.Si(self)
        chemCa = elements.Ca(self)
        chemC = elements.C(self)
        chemMg = elements.Mg(self)
        chemK = elements.K(self)
        chemAl = elements.Al(self)
        chemMn = elements.Mn(self)
        chemFe = elements.Fe(self)
        chemNa = elements.Na(self)
        chemCl = elements.Cl(self)
        chemS = elements.S(self)
        chemF = elements.F(self)
        chemPb = elements.Pb(self)
        chemMo = elements.Mo(self)
        chemCu = elements.Cu(self)
        # [molar mass, density, bulk modulus, shear modulus, vP, vS, GR]
        # Oxides
        chemQuartz = minerals.oxides.quartz("")
        chem_chromite = minerals.oxides.chromite("")
        chem_spinel = minerals.oxides.spinel("")
        chem_magnetite = minerals.oxides.magnetite("")
        chem_hematite = minerals.oxides.hematite("")
        chem_boehmite = minerals.oxides.boehmite("")
        chem_diaspore = minerals.oxides.diaspore("")
        chem_gibbsite = minerals.oxides.gibbsite("")
        chem_uraninite = minerals.oxides.uraninite("")
        chem_cuprite = minerals.oxides.cuprite("")
        # Carbonates
        chemDolomite = minerals.carbonates.dolomite("")
        chemCalcite = minerals.carbonates.calcite("")
        chemSiderite = minerals.carbonates.siderite("")
        chemMagnesite = minerals.carbonates.magnesite("")
        chem_ankerite = minerals.carbonates.ankerite("")
        chem_malachite = minerals.carbonates.malachite("")
        chem_aragonite = minerals.carbonates.aragonite("")
        # Tectosilicates
        chemOrthoclase = minerals.tectosilicates.orthoclase("")
        chemAnorthite = minerals.tectosilicates.anorthite("")
        chemAlbite = minerals.tectosilicates.albite("")
        chemAlkalifeldspar = minerals.feldspars.alkalifeldspar(self, "Alkalifeldspar")
        chemPlagioclase = minerals.feldspars.plagioclase(self, "Plagioclase")
        chem_alkalifeldspar_granite = minerals.feldspars.alkalifeldspar(self, "K")
        chem_plagioclase_granite = minerals.feldspars.plagioclase(self, "Na")
        chem_scapolite = minerals.tectosilicates.scapolite(self, "Scapolite")
        # Inosilicates
        chem_actinolite = minerals.inosilicates.actinolite("")
        chem_tremolite = minerals.inosilicates.tremolite("")
        chem_enstatite = minerals.inosilicates.enstatite("")
        chem_ferrosilite = minerals.inosilicates.ferrosilite("")
        chem_augite = minerals.inosilicates.augite("")
        chem_diopside = minerals.inosilicates.diopside("")
        chem_aegirine = minerals.inosilicates.aegirine("")
        chem_hedenbergite = minerals.inosilicates.hedenbergite("")
        # Phyllosilicates
        chemBiotite = minerals.phyllosilicates.biotite("")
        chemGlauconite = minerals.phyllosilicates.glauconite("")
        chemChlorite = minerals.phyllosilicates.chamosite("")
        chemIllite = minerals.phyllosilicates.illite("")
        chemMuscovite = minerals.phyllosilicates.muscovite("")
        chemKaolinite = minerals.phyllosilicates.kaolinite("")
        chem_biotite = minerals.Biotites.biotite_group(self, "Biotite")
        chem_clinochlore = minerals.phyllosilicates.clinochlore("")
        chem_montmorillonite = minerals.phyllosilicates.montmorillonite("")
        # Nesosilicates
        chem_olivine = minerals.nesosilicates.olivine(self, "Olivine")
        chem_garnet_pyralspite = minerals.nesosilicates.garnet_pyralspite(self)
        chem_garnet_ugrandite = minerals.nesosilicates.garnet_ugrandite(self)
        # Halides
        chemHalite = minerals.halides.halite("")
        chem_fluorite = minerals.halides.fluorite("")
        chemSylvite = minerals.halides.sylvite("")
        chem_carobbiite = minerals.halides.carobbiite("")
        # Sulfates
        chemAnhydrite = minerals.sulfates.anhydrite("")
        chemGypsum = minerals.sulfates.gypsum("")
        chem_scheelite = minerals.sulfates.scheelite("")
        chem_barite = minerals.sulfates.barite("")
        # Sulfides
        chemPyrite = minerals.sulfides.pyrite("")
        chem_bornite = minerals.sulfides.bornite("")
        chemGalena = minerals.sulfides.galena("")
        chemChalcopyrite = minerals.sulfides.chalcopyrite("")
        chemMolybdenite = minerals.sulfides.molybdenite("")
        chem_stibnite = minerals.sulfides.stibnite("")
        chem_arsenopyrite = minerals.sulfides.arsenopyrite("")
        chem_acanthite = minerals.sulfides.acanthite("")
        chem_argentite = minerals.sulfides.argentite("")
        chem_alabandite = minerals.sulfides.alabandite("")
        chem_berthierite = minerals.sulfides.berthierite("")
        chem_pyrrhotite = minerals.sulfides.pyrrhotite("")
        chem_cobaltite = minerals.sulfides.cobaltite("")
        chem_carrollite = minerals.sulfides.carrollite("")
        chem_chalcocite = minerals.sulfides.chalcocite("")
        chem_digenite = minerals.sulfides.digenite("")
        # Antimonides
        chem_allargentum = minerals.antimonides.allargentum("")
        chem_dyscrasite = minerals.antimonides.dyscrasite("")
        # Arsenides
        # Tourmalines
        chem_schorl = minerals.Tourmalines.schorl("")
        chem_dravite = minerals.Tourmalines.dravite("")
        chem_elbaite = minerals.Tourmalines.elbaite("")
        chem_liddicoatite = minerals.Tourmalines.liddicoatite("")
        # Natives
        chem_arsenic = minerals.natives.arsenic("")
        chem_bismuth = minerals.natives.bismuth("")
        # Coltans
        chem_columbite = minerals.Coltans.columbite("")
        chem_tantalite = minerals.Coltans.tantalite("")
        chem_coltan = minerals.Coltans.coltan("")

        #print(chemCalcite)
        #print(chem_aragonite)

        # print("Evaporites")
        # data = evaporites.Evaporites("water", 2500)
        # data_rocksalt_01 = data.create_simple_rocksalt()
        # print(data_rocksalt_01)
        # data_anhydrite_01 = data.create_simple_anhydrite()
        # print(data_anhydrite_01)

        # print("Soils")
        # data = siliciclastics.Soil()
        # data_soil_01 = data.create_simple_soil()
        # print(data_soil_01)
        # data_sand_01 = data.create_simple_sand()
        # print(data_sand_01)

        # print("Igneous rocks: Plutonics")
        # data = igneous.Plutonic("water", 4000)
        # data_granite_01 = data.create_simple_granite()
        # print(data_granite_01)
        # data_syenite_01 = data.create_simple_syenite()
        # print(data_syenite_01)
        # data_monzonite_01 = data.create_simple_monzonite()
        # print(data_monzonite_01)
        # data_gabbro_01 = data.create_simple_gabbro()
        # print(data_gabbro_01)
        # data_diorite_01 = data.create_simple_diorite()
        # print(data_diorite_01)
        # data_granodiorite_01 = data.create_simple_granodiorite()
        # print(data_granodiorite_01)
        # data_tonalite_01 = data.create_simple_tonalite()
        # print(data_tonalite_01)
        # data_quartzrich_granitoid_01 = data.create_simple_quartzrich_granitoid()
        # print(data_quartzrich_granitoid_01)
        # data_quartzolite_01 = data.create_simple_quartzolite()
        # print(data_quartzolite_01)

        # print("Igneous rocks: Volcanics")
        # data = igneous.Volcanic("water", 4000)
        # data_basalt_01 = data.create_simple_basalt()
        # print(data_basalt_01)
        # data_andesite_01 = data.create_simple_andesite()
        # print(data_andesite_01)
        # data_latite_01 = data.create_simple_latite()
        # print(data_latite_01)
        # data_trachyte_01 = data.create_simple_trachyte()
        # print(data_trachyte_01)
        # data_rhyolite_01 = data.create_simple_rhyolite()
        # print(data_rhyolite_01)
        # data_dacite_01 = data.create_simple_dacite()
        # print(data_dacite_01)

        # data = igneous.Plutonic("water", 4000)
        # data_granite_01 = data.create_simple_granite()
        # print(data_granite_01)
        # data_granite_01a = data.create_simple_granite(amounts=data_granite_01[8])
        # print(data_granite_01a)
        # data_granite_01b = data.create_simple_granite(amounts=data_granite_01a[8])
        # print(data_granite_01b)
        # data_granite_01c = data.create_simple_granite(amounts=data_granite_01b[8])
        # print(data_granite_01c)
        # data_granite_01d = data.create_simple_granite(amounts=data_granite_01c[8])
        # print(data_granite_01d)
        #
        # data = siliciclastics.sandstone("water", 100)
        # data_sandstone1 = data.create_simple_sandstone()
        # print(data_sandstone1)
        # data_sandstone1a = data.create_simple_sandstone(amounts=data_sandstone1[8])
        # print(data_sandstone1a)
        # data_sandstone1b = data.create_simple_sandstone(amounts=data_sandstone1a[8])
        # print(data_sandstone1b)
        # data_sandstone1c = data.create_simple_sandstone(amounts=data_sandstone1b[8])
        # print(data_sandstone1c)
        # data_sandstone1d = data.create_simple_sandstone(amounts=data_sandstone1c[8])
        # print(data_sandstone1d)
        # data_sandstone2 = data.create_simple_sandstone(w_Fe=0.1)
        # data_sandstone3 = data.create_simple_sandstone(w_Fe=0.2)
        # data_sandstone4 = data.create_simple_sandstone(w_Fe=0.3)
        # data_sandstone5 = data.create_simple_sandstone(w_Fe=0.4)
        # data_sandstone6 = data.create_simple_sandstone(w_Fe=0.5)
        # print(data_sandstone1)
        # print(data_sandstone2)
        # print(data_sandstone3)
        # print(data_sandstone4)
        # print(data_sandstone5)
        # print(data_sandstone6)
        #
        #data = siliciclastics.shale()
        #data_shale1 = data.create_simple_shale()
        #print(data_shale1)
        #data_shale1a = data.create_simple_shale(amounts=data_shale1[8])
        #print(data_shale1a)
        #data_shale1b = data.create_simple_shale(amounts=data_shale1a[8])
        #print(data_shale1b)
        #data_shale1c = data.create_simple_shale(amounts=data_shale1b[8])
        #print(data_shale1c)
        #data_shale1d = data.create_simple_shale(amounts=data_shale1c[8])
        #print(data_shale1d)
        # data_shale2 = data.create_simple_shale()
        # data_shale3 = data.create_simple_shale()
        # data_shale4 = data.create_simple_shale(w_C=0.05)
        # print(data_shale1)
        # print(data_shale2)
        # print(data_shale3)
        # print(data_shale4)
        #
        # data = carbonates.limestone("water", 100)
        # data_limestone1 = data.create_simple_limestone()
        # print(data_limestone1)
        # data_limestone1a = data.create_simple_limestone(amounts=data_limestone1[8])
        # print(data_limestone1a)
        # data_limestone1b = data.create_simple_limestone(amounts=data_limestone1a[8])
        # print(data_limestone1b)
        # data_limestone1c = data.create_simple_limestone(amounts=data_limestone1b[8])
        # print(data_limestone1c)
        # data_limestone1d = data.create_simple_limestone(amounts=data_limestone1c[8])
        # print(data_limestone1d)
        #data_limestone2 = data.create_simple_limestone(w_K=0.01)
        #print(data_limestone2)
        #data_limestone3 = data.create_simple_limestone(w_Mg=0.05)
        #data_limestone4 = data.create_simple_limestone(w_K=0.05)
        #data_limestone5 = data.create_simple_limestone(w_Ca=0.3)
        #data_limestone6 = data.create_simple_limestone(w_Fe=0.05)
        #print(data_limestone3)
        #print(data_limestone4)
        #print(data_limestone5)
        #print(data_limestone6)
        #
        # data = carbonates.dolomite("water", 100)
        # data_dolomite_01 = data.create_simple_dolomite()
        # print(data_dolomite_01)
        # data_dolomite_01a = data.create_simple_dolomite(amounts=data_dolomite_01[8])
        # print(data_dolomite_01a)
        # data_dolomite_01b = data.create_simple_dolomite(amounts=data_dolomite_01a[8])
        # print(data_dolomite_01b)
        # data_dolomite_01c = data.create_simple_dolomite(amounts=data_dolomite_01b[8])
        # print(data_dolomite_01c)
        # data_dolomite_01d = data.create_simple_dolomite(amounts=data_dolomite_01c[8])
        # print(data_dolomite_01d)
        #data_dolomite_02 = data.create_simple_dolomite(w_Mg=0.1)
        #print(data_dolomite_02)
        #data_dolomite_03 = data.create_simple_dolomite(w_Ca=0.15)
        #print(data_dolomite_03)
        #data_dolomite_04 = data.create_simple_dolomite(w_Fe=0.025)
        #print(data_dolomite_04)
        #
        #data = fluids.Hydrocarbons()
        #data_oil_01 = data.oil()
        #print(data_oil_01)
        #data_gas_01 = data.natural_gas()
        #print(data_gas_01)
        #
        #data = fluids.Water()
        #data_water_01 = data.water()
        #print(data_water_01)
        #
        mineralsSandstone = [chemQuartz, chemCalcite, chemDolomite, chemAlkalifeldspar, chemPlagioclase, chemBiotite, chemGlauconite]
        mineralsShale = [chemQuartz, chemKaolinite, chemChlorite, chemIllite, chemCalcite, chemDolomite, chemAlkalifeldspar, chemPlagioclase]
        mineralsLimestone = [chemCalcite, chemDolomite, chemSiderite, chemQuartz, chemAlkalifeldspar, chemPlagioclase, chemKaolinite, chemChlorite, chemIllite]
        mineralsDolomite = [chemCalcite, chemDolomite, chemMagnesite, chemQuartz, chemAlkalifeldspar, chemPlagioclase, chemKaolinite, chemChlorite, chemIllite]
        mineralsGranite = [chemQuartz, chem_alkalifeldspar_granite, chem_plagioclase_granite, chemBiotite, chemMuscovite, chem_actinolite, chem_tremolite]
        mineralsHalite = [chemHalite, chemAnhydrite, chemGypsum, chemSylvite, chemIllite]
        mineralsAnhydrite = [chemAnhydrite, chemHalite, chemCalcite, chemGalena, chemChalcopyrite, chemMolybdenite, chemPyrite]
        minerals_ore_Fe = [chem_magnetite, chem_hematite, chemQuartz, chemCalcite, chemAlkalifeldspar, chemPlagioclase]
        minerals_basalt = [chemPlagioclase, chem_enstatite, chem_ferrosilite, chem_biotite, chem_actinolite, chem_tremolite, chem_olivine]
        #
        # input = core --> createCore
        lithologies = []
        mineralogy = []
        massMinerals = []
        for i in range(2, len(self.input)):
            massMinerals.append([self.input[i][0], self.input[i][2]])
            if self.input[i][0] not in lithologies:
                lithologies.append(self.input[i][0])
        for i in range(0, len(massMinerals)):
            if massMinerals[i][0] == "sandstone":
                for j in range(0, len(mineralsSandstone)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "shale":
                for j in range(0, len(mineralsShale)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "limestone":
                for j in range(0, len(mineralsLimestone)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "dolomite":
                for j in range(0, len(mineralsDolomite)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "granite":
                for j in range(0, len(mineralsGranite)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "halite":
                for j in range(0, len(mineralsHalite)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "anhydrite":
                for j in range(0, len(mineralsAnhydrite)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "ore":
                for j in range(0, len(minerals_ore_Fe)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
            elif massMinerals[i][0] == "basalt":
                for j in range(0, len(minerals_basalt)):
                    massMinerals[i].extend([[self.input[i+2][4][j][0], self.input[i+2][4][j][1]*self.input[i+2][2]]])
                    if self.input[i+2][4][j][0] not in mineralogy:
                        mineralogy.append(self.input[i+2][4][j][0])
        #for i in range(0, len(massMinerals)):
        #    if massMinerals[i][0] == "sandstone":
        #        print(massMinerals[i])
        lithologies = sorted(lithologies, key=str.lower)
        #print(lithologies)
        mineralogy = sorted(mineralogy, key=str.lower)
        #print(mineralogy)
        #
        n = []
        for i in range(0, len(massMinerals)):
            n.append([massMinerals[i][0]])
        for i in range(0, len(n)):
            if n[i][0] == "sandstone":
                for j in range(0, len(mineralsSandstone)):
                    if mineralsSandstone[j][0] in ["Cal", "Dol", "Qz", "Or", "Hl", "Kln", "Bt", "Glt", "Chl", "Ilt", "Sd", "Mgs", "Anh", "Gp", "Hl", "Tr", "Gn", "Py", "Ccp", "Mol", "Syl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsSandstone[j][1])]])
                    elif mineralsSandstone[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsSandstone[j][1][0]), mineralsSandstone[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsSandstone[j][0])]])
            elif n[i][0] == "shale":
                for j in range(0, len(mineralsShale)):
                    if mineralsShale[j][0] in ["Cal", "Dol", "Qz", "Or", "Hl", "Kln", "Bt", "Glt", "Chl", "Ilt", "Sd", "Mgs", "Anh", "Gp", "Hl", "Tr", "Gn", "Py", "Ccp", "Mol", "Syl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsShale[j][1])]])
                    elif mineralsShale[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsShale[j][1][0]), mineralsShale[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsShale[j][0])]])
            elif n[i][0] == "limestone":
                for j in range(0, len(mineralsLimestone)):
                    if mineralsLimestone[j][0] in ["Cal", "Dol", "Qz", "Or", "Hl", "Kln", "Bt", "Glt", "Chl", "Ilt", "Sd", "Mgs", "Anh", "Gp", "Hl", "Tr", "Gn", "Py", "Ccp", "Mol", "Syl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsLimestone[j][1])]])
                    elif mineralsLimestone[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsLimestone[j][1][0]), mineralsLimestone[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsLimestone[j][0])]])
            elif n[i][0] == "dolomite":
                for j in range(0, len(mineralsDolomite)):
                    if mineralsDolomite[j][0] in ["Cal", "Dol", "Qz", "Or", "Hl", "Kln", "Bt", "Glt", "Chl", "Ilt", "Sd", "Mgs", "Anh", "Gp", "Hl", "Tr", "Gn", "Py", "Ccp", "Mol", "Syl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsDolomite[j][1])]])
                    elif mineralsDolomite[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsDolomite[j][1][0]), mineralsDolomite[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsDolomite[j][0])]])
            elif n[i][0] == "granite":
                for j in range(0, len(mineralsGranite)):
                    if mineralsGranite[j][0] in ["Qz", "Ms", "Tr"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsGranite[j][1])]])
                    elif mineralsGranite[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsGranite[j][1][0]), mineralsGranite[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsGranite[j][0])]])
            elif n[i][0] == "halite":
                for j in range(0, len(mineralsHalite)):
                    if mineralsHalite[j][0] in ["Cal", "Dol", "Qz", "Or", "Hl", "Kln", "Bt", "Glt", "Chl", "Ilt", "Sd", "Mgs", "Anh", "Gp", "Hl", "Tr", "Gn", "Py", "Ccp", "Mol", "Syl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsHalite[j][1])]])
                    elif mineralsHalite[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsHalite[j][1][0]), mineralsHalite[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsHalite[j][0])]])
            elif n[i][0] == "anhydrite":
                for j in range(0, len(mineralsAnhydrite)):
                    if mineralsAnhydrite[j][0] in ["Cal", "Dol", "Qz", "Or", "Hl", "Kln", "Bt", "Glt", "Chl", "Ilt", "Sd", "Mgs", "Anh", "Gp", "Hl", "Tr", "Gn", "Py", "Ccp", "Mol", "Syl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsAnhydrite[j][1])]])
                    elif mineralsAnhydrite[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsAnhydrite[j][1][0]), mineralsAnhydrite[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*mineralsAnhydrite[j][0])]])
            elif n[i][0] == "ore":
                for j in range(0, len(minerals_ore_Fe)):
                    if minerals_ore_Fe[j][0] in ["Mag", "Hem", "Qz", "Cal"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_ore_Fe[j][1])]])
                    elif minerals_ore_Fe[j][0] in ["Kfs", "Afs", "Pl"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_ore_Fe[j][1][0]), minerals_ore_Fe[j][1][1]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_ore_Fe[j][0])]])
            elif n[i][0] == "basalt":
                for j in range(0, len(minerals_basalt)):
                    if minerals_basalt[j][0] in ["En", "Fs", "Tr"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_basalt[j][1])]])
                    elif minerals_basalt[j][0] in ["Kfs", "Afs", "Pl", "Act"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_basalt[j][1][0]), minerals_basalt[j][1][1]]])
                    elif minerals_basalt[j][0] in ["Bt", "Ol"]:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_basalt[j][1][0]), minerals_basalt[j][1][1], minerals_basalt[j][1][2], minerals_basalt[j][1][3]]])
                    else:
                        n[i].extend([[massMinerals[i][j+2][0], massMinerals[i][j+2][1]/(0.001*minerals_basalt[j][0])]])
        #
        data = elementanalysis(n, massMinerals)
        if "Qz" in mineralogy:
            dataQz = data.analyzeQuartz()
        if "Cal" in mineralogy:
            dataCal = data.analyzeCalcite()
        if "Dol" in mineralogy:
            dataDol = data.analyzeDolomite()
        if "Bt" in mineralogy:
            data_Bt = data.analyze_biotite(self.sequences)
        if "Ilt" in mineralogy:
            dataIlt = data.analyzeIllite()
        if "Chl" in mineralogy:
            dataChl = data.analyzeChlorite()
        if "Glt" in mineralogy:
            dataGlt = data.analyzeGlauconite()
        if "Ms" in mineralogy:
            dataMs = data.analyzeMuscovite()
        if "Py" in mineralogy:
            dataPy = data.analyzePyrite()
        if "Hl" in mineralogy:
            dataHl = data.analyzeHalite()
        if "Anh" in mineralogy:
            dataAnh = data.analyzeAnhydrite()
        if "Gp" in mineralogy:
            dataGp = data.analyzeGypsum()
        if "Syl" in mineralogy:
            dataSyl = data.analyzeSylvite()
        if "Gn" in mineralogy:
            dataGn = data.analyzeGalena()
        if "Ccp" in mineralogy:
            dataCcp = data.analyzeChalcopyrite()
        if "Mol" in mineralogy:
            dataMol = data.analyzeMolybdenite()
        if "Kln" in mineralogy:
            dataKln = data.analyzeKaolinite()
        if "Afs" in mineralogy or "Kfs" in mineralogy:
            dataKfs = data.analyzeAlkalifeldspar(self.sequences)
        if "Pl" in mineralogy:
            dataPl = data.analyzePlagioclase(self.sequences)
        if "Sd" in mineralogy:
            dataSd = data.analyzeSiderite()
        if "Mgs" in mineralogy:
            dataMgs = data.analyzeMagnesite()
        if "Act" in mineralogy:
            dataAct = data.analyzeActinolite(self.sequences)
        if "Tr" in mineralogy:
            dataTr = data.analyzeTremolite()
        if "Mag" in mineralogy:
            data_Mag = data.analyzeMagnetite()
        if "Hem" in mineralogy:
            data_Hem = data.analyzeHematite()
        if "En" in mineralogy:
            data_En = data.analyze_enstatite()
        if "Fs" in mineralogy:
            data_Fs = data.analyze_ferrosilite()
        if "Ol" in mineralogy:
            data_Ol = data.analyze_olivine(self.sequences)
        #
        sandstone = []
        limestone = []
        dolomite = []
        shale = []
        granite = []
        halite = []
        anhydrite = []
        ore_Fe = []
        basalt = []
        for i in range(0, len(n)):
            if n[i][0] == "sandstone":
                sandstone.append([])
            elif n[i][0] == "shale":
                shale.append([])
            elif n[i][0] == "limestone":
                limestone.append([])
            elif n[i][0] == "dolomite":
                dolomite.append([])
            elif n[i][0] == "granite":
                granite.append([])
            elif n[i][0] == "halite":
                halite.append([])
            elif n[i][0] == "anhydrite":
                anhydrite.append([])
            elif n[i][0] == "ore":
                ore_Fe.append([])
            elif n[i][0] == "basalt":
                basalt.append([])
        #
        results = []
        #
        if "sandstone" in lithologies:
            ## SANDSTONE
            #
            QzSandstone = []
            CalSandstone = []
            DolSandstone = []
            KfsSandstone = []
            PlSandstone = []
            BtSandstone = []
            GltSandstone = []
            for i in range(0, len(dataQz)):
                if dataQz[i][0] == "sandstone":
                    QzSandstone.append(dataQz[i][1])
            for i in range(0, len(dataCal)):
                if dataCal[i][0] == "sandstone":
                    CalSandstone.append(dataCal[i][1])
            for i in range(0, len(dataDol)):
                if dataDol[i][0] == "sandstone":
                    DolSandstone.append(dataDol[i][1])
            for i in range(0, len(dataKfs)):
                if dataKfs[i][0] == "sandstone":
                    KfsSandstone.append(dataKfs[i][1])
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "sandstone":
                    PlSandstone.append(dataPl[i][1])
            for i in range(0, len(data_Bt)):
                if data_Bt[i][0] == "sandstone":
                    BtSandstone.append(data_Bt[i][1])
            for i in range(0, len(dataGlt)):
                if dataGlt[i][0] == "sandstone":
                    GltSandstone.append(dataGlt[i][1])
            for i in range(0, len(sandstone)):
                sandstone[i].extend(QzSandstone[i])
                sandstone[i].extend(CalSandstone[i])
                sandstone[i].extend(DolSandstone[i])
                sandstone[i].extend(KfsSandstone[i])
                sandstone[i].extend(PlSandstone[i])
                sandstone[i].extend(BtSandstone[i])
                sandstone[i].extend(GltSandstone[i])
            elementsSandstone = []
            for i in range(0, len(sandstone)):
                elementsSandstone.append([])
            for i in range(0, len(sandstone)):
                for j in range(0, len(sandstone[i])):
                    if [sandstone[i][j][0],0] in elementsSandstone[i]:
                        pass
                    else:
                        elementsSandstone[i].append([sandstone[i][j][0], 0])
            #print("Elements Sandstone:", elementsSandstone[0])
            for i in range(0, len(sandstone)):
                for j in range(0, len(sandstone[i])):
                    for k in range(0, len(elementsSandstone[i])):
                        if elementsSandstone[i][k][0] == sandstone[i][j][0]:
                            elementsSandstone[i][k][1] += sandstone[i][j][1]
                        else:
                            pass
            mSandstone = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "sandstone":
                    mSandstone.append(round(self.input[i][2],3))
            #print("m:", mSandstone)
            dataM = []
            for i in range(0, len(elementsSandstone)):
                m = 0
                for j in range(0, len(elementsSandstone[i])):
                    m += round(elementsSandstone[i][j][1], 3)
                dataM.append(round(m, 3))
            #print("dataM", dataM)
            for i in range(0, len(elementsSandstone)):
                for j in range(0, len(elementsSandstone[i])):
                    if elementsSandstone[i][j][0] == "H":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemH[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "C":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemC[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "O":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemO[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "F":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemF[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "Na":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemNa[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "Mg":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemMg[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "Al":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemAl[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "Si":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemSi[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "K":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemK[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "Ca":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemCa[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
                    elif elementsSandstone[i][j][0] == "Fe":
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/(0.001*chemFe[2]))
                        elementsSandstone[i][j].append(elementsSandstone[i][j][1]/dataM[i]*100)
            resultsSandstone = []
            for i in range(0, len(elementsSandstone)):
                resultsSandstone.append(["sandstone", round(dataM[i],3)])
            for i in range(0, len(elementsSandstone)):
                for j in range(0, len(elementsSandstone[i])):
                    resultsSandstone[i].append([elementsSandstone[i][j][0], round(elementsSandstone[i][j][1],3), round(elementsSandstone[i][j][2],1), round(elementsSandstone[i][j][3],2)])
            #
            results.append(resultsSandstone)
            #print("Results Sandstone:")
            #for j in range(0, len(resultsSandstone)):
            #    print(j+1, ":", resultsSandstone[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Sandstone = open("./outputs/Geochemistry_Sandstone.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Sandstone.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("C")+","+str("O")+","+str("F")+","+str("Na")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(resultsSandstone)):
                Geochemistry_Sandstone.write(str(resultsSandstone[i][0])+","+str(resultsSandstone[i][1])+","+str(resultsSandstone[i][11][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][5][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][3][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][12][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][7][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][6][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][9][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][2][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][8][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][4][1]/resultsSandstone[i][1])+","+str(resultsSandstone[i][10][1]/resultsSandstone[i][1])+"\n")
        #
        if "limestone" in lithologies:
            ## LIMESTONE
            #
            CalLimestone = []
            DolLimestone = []
            QzLimestone = []
            SdLimestone = []
            KfsLimestone = []
            PlLimestone = []
            KlnLimestone = []
            ChlLimestone = []
            IltLimestone = []
            for i in range(0, len(dataQz)):
                if dataQz[i][0] == "limestone":
                    QzLimestone.append(dataQz[i][1])
            for i in range(0, len(dataCal)):
                if dataCal[i][0] == "limestone":
                    CalLimestone.append(dataCal[i][1])
            for i in range(0, len(dataDol)):
                if dataDol[i][0] == "limestone":
                    DolLimestone.append(dataDol[i][1])
            for i in range(0, len(dataSd)):
                if dataSd[i][0] == "limestone":
                    SdLimestone.append(dataSd[i][1])
            for i in range(0, len(dataKfs)):
                if dataKfs[i][0] == "limestone":
                    KfsLimestone.append(dataKfs[i][1])
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "limestone":
                    PlLimestone.append(dataPl[i][1])
            for i in range(0, len(dataKln)):
                if dataKln[i][0] == "limestone":
                    KlnLimestone.append(dataKln[i][1])
            for i in range(0, len(dataChl)):
                if dataChl[i][0] == "limestone":
                    ChlLimestone.append(dataChl[i][1])
            for i in range(0, len(dataIlt)):
                if dataIlt[i][0] == "limestone":
                    IltLimestone.append(dataIlt[i][1])
            for i in range(0, len(limestone)):
                limestone[i].extend(CalLimestone[i])
                limestone[i].extend(DolLimestone[i])
                limestone[i].extend(QzLimestone[i])
                limestone[i].extend(SdLimestone[i])
                limestone[i].extend(KfsLimestone[i])
                limestone[i].extend(PlLimestone[i])
                limestone[i].extend(KlnLimestone[i])
                limestone[i].extend(ChlLimestone[i])
                limestone[i].extend(IltLimestone[i])
            elementsLimestone = []
            for i in range(0, len(limestone)):
                elementsLimestone.append([])
            for i in range(0, len(limestone)):
                for j in range(0, len(limestone[i])):
                    if [limestone[i][j][0],0] in elementsLimestone[i]:
                        pass
                    else:
                        elementsLimestone[i].append([limestone[i][j][0], 0])
            #print("Elements Limestone:", elementsLimestone[0])
            for i in range(0, len(limestone)):
                for j in range(0, len(limestone[i])):
                    for k in range(0, len(elementsLimestone[i])):
                        if elementsLimestone[i][k][0] == limestone[i][j][0]:
                            elementsLimestone[i][k][1] += limestone[i][j][1]
                        else:
                            pass
            mLimestone = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "limestone":
                    mLimestone.append(round(self.input[i][2],3))
            #print("m:", mLimestone)
            dataM = []
            for i in range(0, len(elementsLimestone)):
                m = 0
                for j in range(0, len(elementsLimestone[i])):
                    m += round(elementsLimestone[i][j][1], 3)
                dataM.append(round(m, 3))
            # print("dataM", dataM)
            for i in range(0, len(elementsLimestone)):
                for j in range(0, len(elementsLimestone[i])):
                    if elementsLimestone[i][j][0] == "Ca":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemCa[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "C":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemC[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "O":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemO[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "Mg":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemMg[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "Si":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemSi[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "K":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemK[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "Al":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemAl[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "Fe":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemFe[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "Na":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemNa[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
                    elif elementsLimestone[i][j][0] == "H":
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/(0.001*chemH[2]))
                        elementsLimestone[i][j].append(elementsLimestone[i][j][1]/dataM[i]*100)
            resultsLimestone = []
            for i in range(0, len(elementsLimestone)):
                resultsLimestone.append(["limestone", round(dataM[i],3)])
            for i in range(0, len(elementsLimestone)):
                for j in range(0, len(elementsLimestone[i])):
                    resultsLimestone[i].append([elementsLimestone[i][j][0], round(elementsLimestone[i][j][1],3), round(elementsLimestone[i][j][2],1), round(elementsLimestone[i][j][3],2)])
            #
            results.append(resultsLimestone)
            #print("Results Limestone:")
            #for j in range(0, len(resultsLimestone)):
            #    print(j+1, ":", resultsLimestone[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Limestone = open("./outputs/Geochemistry_Limestone.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Limestone.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("C")+","+str("O")+","+str("Na")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(resultsLimestone)):
                Geochemistry_Limestone.write(str(resultsLimestone[i][0])+","+str(resultsLimestone[i][1])+","+str(resultsLimestone[i][11][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][3][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][4][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][8][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][5][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][10][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][6][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][9][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][2][1]/resultsLimestone[i][1])+","+str(resultsLimestone[i][7][1]/resultsLimestone[i][1])+"\n")
        #
        if "dolomite" in lithologies:
            ## DOLOMITE
            #
            CalDolomite = []
            DolDolomite = []
            QzDolomite = []
            MgsDolomite = []
            KfsDolomite = []
            PlDolomite = []
            KlnDolomite = []
            ChlDolomite = []
            IltDolomite = []
            for i in range(0, len(dataQz)):
                if dataQz[i][0] == "dolomite":
                    QzDolomite.append(dataQz[i][1])
            for i in range(0, len(dataCal)):
                if dataCal[i][0] == "dolomite":
                    CalDolomite.append(dataCal[i][1])
            for i in range(0, len(dataDol)):
                if dataDol[i][0] == "dolomite":
                    DolDolomite.append(dataDol[i][1])
            for i in range(0, len(dataMgs)):
                if dataMgs[i][0] == "dolomite":
                    MgsDolomite.append(dataMgs[i][1])
            for i in range(0, len(dataKfs)):
                if dataKfs[i][0] == "dolomite":
                    KfsDolomite.append(dataKfs[i][1])
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "dolomite":
                    PlDolomite.append(dataPl[i][1])
            for i in range(0, len(dataKln)):
                if dataKln[i][0] == "dolomite":
                    KlnDolomite.append(dataKln[i][1])
            for i in range(0, len(dataChl)):
                if dataChl[i][0] == "dolomite":
                    ChlDolomite.append(dataChl[i][1])
            for i in range(0, len(dataIlt)):
                if dataIlt[i][0] == "dolomite":
                    IltDolomite.append(dataIlt[i][1])
            for i in range(0, len(dolomite)):
                dolomite[i].extend(CalDolomite[i])
                dolomite[i].extend(DolDolomite[i])
                dolomite[i].extend(QzDolomite[i])
                dolomite[i].extend(MgsDolomite[i])
                dolomite[i].extend(KfsDolomite[i])
                dolomite[i].extend(PlDolomite[i])
                dolomite[i].extend(KlnDolomite[i])
                dolomite[i].extend(ChlDolomite[i])
                dolomite[i].extend(IltDolomite[i])
            elementsDolomite = []
            for i in range(0, len(dolomite)):
                elementsDolomite.append([])
            for i in range(0, len(dolomite)):
                for j in range(0, len(dolomite[i])):
                    if [dolomite[i][j][0],0] in elementsDolomite[i]:
                        pass
                    else:
                        elementsDolomite[i].append([dolomite[i][j][0], 0])
            #print("Elements Dolomite:", elementsDolomite[0])
            for i in range(0, len(dolomite)):
                for j in range(0, len(dolomite[i])):
                    for k in range(0, len(elementsDolomite[i])):
                        if elementsDolomite[i][k][0] == dolomite[i][j][0]:
                            elementsDolomite[i][k][1] += dolomite[i][j][1]
                        else:
                            pass
            mDolomite = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "dolomite":
                    mDolomite.append(round(self.input[i][2],3))
            #print("m:", mDolomite)
            dataM = []
            for i in range(0, len(elementsDolomite)):
                m = 0
                for j in range(0, len(elementsDolomite[i])):
                    m += round(elementsDolomite[i][j][1], 3)
                dataM.append(round(m, 3))
            # print("dataM", dataM)
            for i in range(0, len(elementsDolomite)):
                for j in range(0, len(elementsDolomite[i])):
                    if elementsDolomite[i][j][0] == "Ca":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemCa[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "C":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemC[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "O":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemO[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "Mg":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemMg[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "Si":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemSi[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "K":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemK[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "Al":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemAl[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "Fe":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemFe[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "Na":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemNa[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
                    elif elementsDolomite[i][j][0] == "H":
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/(0.001*chemH[2]))
                        elementsDolomite[i][j].append(elementsDolomite[i][j][1]/dataM[i]*100)
            resultsDolomite = []
            for i in range(0, len(elementsDolomite)):
                resultsDolomite.append(["dolomite", round(dataM[i],3)])
            for i in range(0, len(elementsDolomite)):
                for j in range(0, len(elementsDolomite[i])):
                    resultsDolomite[i].append([elementsDolomite[i][j][0], round(elementsDolomite[i][j][1],3), round(elementsDolomite[i][j][2],1), round(elementsDolomite[i][j][3],2)])
            #
            results.append(resultsDolomite)
            #print("Results Dolomite:")
            #for j in range(0, len(resultsDolomite)):
            #    print(j+1, ":", resultsDolomite[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Dolomite = open("./outputs/Geochemistry_Dolomite.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Dolomite.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("C")+","+str("O")+","+str("Na")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(resultsDolomite)):
                Geochemistry_Dolomite.write(str(resultsDolomite[i][0])+","+str(resultsDolomite[i][1])+","+str(resultsDolomite[i][10][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][3][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][4][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][7][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][5][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][9][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][6][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][8][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][2][1]/resultsDolomite[i][1])+","+str(resultsDolomite[i][11][1]/resultsDolomite[i][1])+"\n")
        #
        if "shale" in lithologies:
            ## SHALE
            #
            QzShale = []
            KlnShale = []
            ChlShale = []
            IltShale = []
            CalShale = []
            DolShale = []
            KfsShale = []
            PlShale = []
            for i in range(0, len(dataQz)):
                if dataQz[i][0] == "shale":
                    QzShale.append(dataQz[i][1])
            for i in range(0, len(dataKln)):
                if dataKln[i][0] == "shale":
                    KlnShale.append(dataKln[i][1])
            for i in range(0, len(dataChl)):
                if dataChl[i][0] == "shale":
                    ChlShale.append(dataChl[i][1])
            for i in range(0, len(dataIlt)):
                if dataIlt[i][0] == "shale":
                    IltShale.append(dataIlt[i][1])
            for i in range(0, len(dataCal)):
                if dataCal[i][0] == "shale":
                    CalShale.append(dataCal[i][1])
            for i in range(0, len(dataDol)):
                if dataDol[i][0] == "shale":
                    DolShale.append(dataDol[i][1])
            for i in range(0, len(dataKfs)):
                if dataKfs[i][0] == "shale":
                    KfsShale.append(dataKfs[i][1])
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "shale":
                    PlShale.append(dataPl[i][1])
            for i in range(0, len(shale)):
                shale[i].extend(QzShale[i])
                shale[i].extend(KlnShale[i])
                shale[i].extend(ChlShale[i])
                shale[i].extend(IltShale[i])
                shale[i].extend(CalShale[i])
                shale[i].extend(DolShale[i])
                shale[i].extend(KfsShale[i])
                shale[i].extend(PlShale[i])
            elementsShale = []
            for i in range(0, len(shale)):
                elementsShale.append([])
            for i in range(0, len(shale)):
                for j in range(0, len(shale[i])):
                    if [shale[i][j][0],0] in elementsShale[i]:
                        pass
                    else:
                        elementsShale[i].append([shale[i][j][0], 0])
            #print("Elements Shale:", elementsShale[0])
            for i in range(0, len(shale)):
                for j in range(0, len(shale[i])):
                    for k in range(0, len(elementsShale[i])):
                        if elementsShale[i][k][0] == shale[i][j][0]:
                            elementsShale[i][k][1] += shale[i][j][1]
                        else:
                            pass
            mShale = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "shale":
                    mShale.append(round(self.input[i][2],3))
            #print("m:", mShale)
            dataM = []
            for i in range(0, len(elementsShale)):
                m = 0
                for j in range(0, len(elementsShale[i])):
                    m += round(elementsShale[i][j][1], 3)
                dataM.append(round(m, 3))
            # print("dataM", dataM)
            for i in range(0, len(elementsShale)):
                for j in range(0, len(elementsShale[i])):
                    if elementsShale[i][j][0] == "Si":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemSi[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "O":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemO[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "K":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemK[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "Mg":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemMg[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "Fe":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemFe[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "Al":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemAl[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "H":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemH[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "Ca":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemCa[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "Na":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemNa[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
                    elif elementsShale[i][j][0] == "C":
                        elementsShale[i][j].append(elementsShale[i][j][1]/(0.001*chemC[2]))
                        elementsShale[i][j].append(elementsShale[i][j][1]/dataM[i]*100)
            resultsShale = []
            for i in range(0, len(elementsShale)):
                resultsShale.append(["shale", round(dataM[i],3)])
            for i in range(0, len(elementsShale)):
                for j in range(0, len(elementsShale[i])):
                    resultsShale[i].append([elementsShale[i][j][0], round(elementsShale[i][j][1],3), round(elementsShale[i][j][2],1), round(elementsShale[i][j][3],2)])
            #
            results.append(resultsShale)
            #print("Results Shale:")
            #for j in range(0, len(resultsShale)):
            #    print(j+1, ":", resultsShale[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Shale = open("./outputs/Geochemistry_Shale.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Shale.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("C")+","+str("O")+","+str("Na")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(resultsShale)):
                Geochemistry_Shale.write(str(resultsShale[i][0])+","+str(resultsShale[i][1])+","+str(resultsShale[i][5][1]/resultsShale[i][1])+","+str(resultsShale[i][10][1]/resultsShale[i][1]/resultsShale[i][1])+","+str(resultsShale[i][3][1]/resultsShale[i][1])+","+str(resultsShale[i][11][1]/resultsShale[i][1])+","+str(resultsShale[i][7][1]/resultsShale[i][1])+","+str(resultsShale[i][4][1]/resultsShale[i][1])+","+str(resultsShale[i][2][1]/resultsShale[i][1])+","+str(resultsShale[i][8][1]/resultsShale[i][1])+","+str(resultsShale[i][9][1]/resultsShale[i][1])+","+str(resultsShale[i][6][1]/resultsShale[i][1])+"\n")
        #
        if "granite" in lithologies:
            ## GRANITE
            #
            QzGranite = []
            KfsGranite = []
            PlGranite = []
            BtGranite = []
            MsGranite = []
            ActGranite = []
            TrGranite = []

            for i in range(0, len(dataQz)):
                if dataQz[i][0] == "granite":
                    QzGranite.append(dataQz[i][1])
            for i in range(0, len(dataKfs)):
                if dataKfs[i][0] == "granite":
                    KfsGranite.append(dataKfs[i][1])
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "granite":
                    PlGranite.append(dataPl[i][1])
            for i in range(0, len(data_Bt)):
                if data_Bt[i][0] == "granite":
                    BtGranite.append(data_Bt[i][1])
            for i in range(0, len(dataMs)):
                if dataMs[i][0] == "granite":
                    MsGranite.append(dataMs[i][1])
            for i in range(0, len(dataAct)):
                if dataAct[i][0] == "granite":
                    ActGranite.append(dataAct[i][1])
            for i in range(0, len(dataTr)):
                if dataTr[i][0] == "granite":
                    TrGranite.append(dataTr[i][1])
            for i in range(0, len(granite)):
                granite[i].extend(QzGranite[i])
                granite[i].extend(KfsGranite[i])
                granite[i].extend(PlGranite[i])
                granite[i].extend(BtGranite[i])
                granite[i].extend(MsGranite[i])
                granite[i].extend(ActGranite[i])
                granite[i].extend(TrGranite[i])
            elementsGranite = []
            for i in range(0, len(granite)):
                elementsGranite.append([])
            for i in range(0, len(granite)):
                for j in range(0, len(granite[i])):
                    if [granite[i][j][0],0] in elementsGranite[i]:
                        pass
                    else:
                        elementsGranite[i].append([granite[i][j][0], 0])
            #print("Elements Granite:", elementsGranite[0])
            for i in range(0, len(granite)):
                for j in range(0, len(granite[i])):
                    for k in range(0, len(elementsGranite[i])):
                        if elementsGranite[i][k][0] == granite[i][j][0]:
                            elementsGranite[i][k][1] += granite[i][j][1]
                        else:
                            pass
            mGranite = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "granite":
                    mGranite.append(round(self.input[i][2],3))
            #print("m:", mGranite)
            dataM = []
            for i in range(0, len(elementsGranite)):
                m = 0
                for j in range(0, len(elementsGranite[i])):
                    m += round(elementsGranite[i][j][1], 3)
                dataM.append(round(m, 3))
            # print("dataM", dataM)
            for i in range(0, len(elementsGranite)):
                for j in range(0, len(elementsGranite[i])):
                    if elementsGranite[i][j][0] == "H":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemH[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "O":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemO[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "F":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemF[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "Na":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemNa[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "Mg":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemMg[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "Al":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemAl[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "Si":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemSi[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "K":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemK[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "Ca":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemCa[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
                    elif elementsGranite[i][j][0] == "Fe":
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/(0.001*chemFe[2]))
                        elementsGranite[i][j].append(elementsGranite[i][j][1]/dataM[i]*100)
            resultsGranite = []
            for i in range(0, len(elementsGranite)):
                resultsGranite.append(["granite", round(mGranite[i],3)])
            for i in range(0, len(elementsGranite)):
                for j in range(0, len(elementsGranite[i])):
                    resultsGranite[i].append([elementsGranite[i][j][0], round(elementsGranite[i][j][1],3), round(elementsGranite[i][j][2],1), round(elementsGranite[i][j][3],2)])
            #
            results.append(resultsGranite)
            #print("Results Granite:")
            #for j in range(0, len(resultsGranite)):
            #    print(j+1, ":", resultsGranite[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Granite = open("./outputs/Geochemistry_Granite.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Granite.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("O")+","+str("F")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("S")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(resultsGranite)):
                Geochemistry_Granite.write(str(resultsGranite[i][0])+","+str(resultsGranite[i][1])+","+str(resultsGranite[i][8][1]/resultsGranite[i][1])+","+str(resultsGranite[i][3][1]/resultsGranite[i][1])+","+str(resultsGranite[i][9][1]/resultsGranite[i][1])+","+str(resultsGranite[i][6][1]/resultsGranite[i][1])+","+str(resultsGranite[i][5][1]/resultsGranite[i][1])+","+str(resultsGranite[i][2][1]/resultsGranite[i][1])+","+str(resultsGranite[i][10][1]/resultsGranite[i][1])+","+str(resultsGranite[i][4][1]/resultsGranite[i][1])+","+str(resultsGranite[i][7][1]/resultsGranite[i][1])+"\n")
        #
        if "halite" in lithologies:
            ## HALITE
            #
            HlHalite = []   # NaCl
            AnhHalite = []  # CaSO4
            GpHalite = []   # CaSO4*H2O
            SylHalite = []  # KCl
            IltHalite = []  # (K,H3,O)(Al,Mg,Fe)2(Si,Al)4O10[(OH)2,(H2O)]

            for i in range(0, len(dataHl)):
                if dataHl[i][0] == "halite":
                    HlHalite.append(dataHl[i][1])
            for i in range(0, len(dataAnh)):
                if dataAnh[i][0] == "halite":
                    AnhHalite.append(dataAnh[i][1])
            for i in range(0, len(dataGp)):
                if dataGp[i][0] == "halite":
                    GpHalite.append(dataGp[i][1])
            for i in range(0, len(dataSyl)):
                if dataSyl[i][0] == "halite":
                    SylHalite.append(dataSyl[i][1])
            for i in range(0, len(dataIlt)):
                if dataIlt[i][0] == "halite":
                    IltHalite.append(dataIlt[i][1])
            for i in range(0, len(halite)):
                halite[i].extend(HlHalite[i])
                halite[i].extend(AnhHalite[i])
                halite[i].extend(GpHalite[i])
                halite[i].extend(SylHalite[i])
                halite[i].extend(IltHalite[i])
            elementsHalite = []
            for i in range(0, len(halite)):
                elementsHalite.append([])
            for i in range(0, len(halite)):
                for j in range(0, len(halite[i])):
                    if [halite[i][j][0],0] in elementsHalite[i]:
                        pass
                    else:
                        elementsHalite[i].append([halite[i][j][0], 0])
            #print("Elements Halite:", elementsHalite[0])
            for i in range(0, len(halite)):
                for j in range(0, len(halite[i])):
                    for k in range(0, len(elementsHalite[i])):
                        if elementsHalite[i][k][0] == halite[i][j][0]:
                            elementsHalite[i][k][1] += halite[i][j][1]
                        else:
                            pass
            mHalite = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "halite":
                    mHalite.append(round(self.input[i][2],3))
            #print("m:", mHalite)
            dataM = []
            for i in range(0, len(elementsHalite)):
                m = 0
                for j in range(0, len(elementsHalite[i])):
                    m += round(elementsHalite[i][j][1], 3)
                dataM.append(round(m, 3))
            # print("dataM", dataM)
            for i in range(0, len(elementsHalite)):
                for j in range(0, len(elementsHalite[i])):
                    if elementsHalite[i][j][0] == "H":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemH[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "O":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemO[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Na":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemNa[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Mg":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemMg[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Al":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemAl[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Si":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemSi[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "S":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemS[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Cl":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemCl[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "K":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemK[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Ca":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemCa[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)
                    elif elementsHalite[i][j][0] == "Fe":
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/(0.001*chemFe[2]))
                        elementsHalite[i][j].append(elementsHalite[i][j][1]/dataM[i]*100)

            resultsHalite = []
            for i in range(0, len(elementsHalite)):
                resultsHalite.append(["halite", round(dataM[i],3)])
            for i in range(0, len(elementsHalite)):
                for j in range(0, len(elementsHalite[i])):
                    resultsHalite[i].append([elementsHalite[i][j][0], round(elementsHalite[i][j][1],3), round(elementsHalite[i][j][2],1), round(elementsHalite[i][j][3],2)])
            #
            results.append(resultsHalite)
            #print("Results Halite:")
            #for j in range(0, len(resultsHalite)):
            #    print(j+1, ":", resultsHalite[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Halite = open("./outputs/Geochemistry_Halite.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Halite.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("O")+","+str("Na")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("S")+","+str("Cl")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(resultsHalite)):
                Geochemistry_Halite.write(str(resultsHalite[i][0])+","+str(resultsHalite[i][1])+","+str(resultsHalite[i][8][1]/resultsHalite[i][1])+","+str(resultsHalite[i][3][1]/resultsHalite[i][1])+","+str(resultsHalite[i][9][1]/resultsHalite[i][1])+","+str(resultsHalite[i][6][1]/resultsHalite[i][1])+","+str(resultsHalite[i][5][1]/resultsHalite[i][1])+","+str(resultsHalite[i][2][1]/resultsHalite[i][1])+","+str(resultsHalite[i][10][1]/resultsHalite[i][1])+","+str(resultsHalite[i][4][1]/resultsHalite[i][1])+","+str(resultsHalite[i][7][1]/resultsHalite[i][1])+"\n")
        #
        if "anhydrite" in lithologies:
            ## ANHYDRITE
            #
            AnhAnhydrite = [] # CaSO4
            HlAnhydrite = []  # NaCl
            CalAnhydrite = []  # CaCO3
            GnAnhydrite = []  # PbS
            CcpAnhydrite = []  # CuFeS2
            MolAnhydrite = []  # MoS2
            PyAnhydrite = []  # FeS2

            for i in range(0, len(dataAnh)):
                if dataAnh[i][0] == "anhydrite":
                    AnhAnhydrite.append(dataAnh[i][1])
            for i in range(0, len(dataHl)):
                if dataHl[i][0] == "anhydrite":
                    HlAnhydrite.append(dataHl[i][1])
            for i in range(0, len(dataCal)):
                if dataCal[i][0] == "anhydrite":
                    CalAnhydrite.append(dataCal[i][1])
            for i in range(0, len(dataGn)):
                if dataGn[i][0] == "anhydrite":
                    GnAnhydrite.append(dataGn[i][1])
            for i in range(0, len(dataCcp)):
                if dataCcp[i][0] == "anhydrite":
                    CcpAnhydrite.append(dataCcp[i][1])
            for i in range(0, len(dataMol)):
                if dataMol[i][0] == "anhydrite":
                    MolAnhydrite.append(dataMol[i][1])
            for i in range(0, len(dataPy)):
                if dataPy[i][0] == "anhydrite":
                    PyAnhydrite.append(dataPy[i][1])
            for i in range(0, len(anhydrite)):
                anhydrite[i].extend(AnhAnhydrite[i])
                anhydrite[i].extend(HlAnhydrite[i])
                anhydrite[i].extend(CalAnhydrite[i])
                anhydrite[i].extend(GnAnhydrite[i])
                anhydrite[i].extend(CcpAnhydrite[i])
                anhydrite[i].extend(MolAnhydrite[i])
                anhydrite[i].extend(PyAnhydrite[i])
            elementsAnhydrite = []
            for i in range(0, len(anhydrite)):
                elementsAnhydrite.append([])
            for i in range(0, len(anhydrite)):
                for j in range(0, len(anhydrite[i])):
                    if [anhydrite[i][j][0], 0] in elementsAnhydrite[i]:
                        pass
                    else:
                        elementsAnhydrite[i].append([anhydrite[i][j][0], 0])
            # print("Elements Anhydrite:", elementsAnhydrite[0])
            for i in range(0, len(anhydrite)):
                for j in range(0, len(anhydrite[i])):
                    for k in range(0, len(elementsAnhydrite[i])):
                        if elementsAnhydrite[i][k][0] == anhydrite[i][j][0]:
                            elementsAnhydrite[i][k][1] += anhydrite[i][j][1]
                        else:
                            pass
            mAnhydrite = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "anhydrite":
                    mAnhydrite.append(round(self.input[i][2], 3))
            # print("m:", mAnhydrite)
            dataM = []
            for i in range(0, len(elementsAnhydrite)):
                m = 0
                for j in range(0, len(elementsAnhydrite[i])):
                    m += round(elementsAnhydrite[i][j][1], 3)
                dataM.append(round(m, 3))
            # print("dataM", dataM)
            for i in range(0, len(elementsAnhydrite)):
                for j in range(0, len(elementsAnhydrite[i])):
                    if elementsAnhydrite[i][j][0] == "C":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemC[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "O":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemO[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Na":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemNa[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "S":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemS[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Cl":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemCl[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Ca":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemCa[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Fe":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemFe[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Cu":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemCu[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Mo":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemMo[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM [i] * 100)
                    elif elementsAnhydrite[i][j][0] == "Pb":
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / (0.001 * chemPb[2]))
                        elementsAnhydrite[i][j].append(elementsAnhydrite[i][j][1] / dataM[i] * 100)

            resultsAnhydrite = []
            for i in range(0, len(elementsAnhydrite)):
                resultsAnhydrite.append(["anhydrite", round(dataM[i], 3)])
            for i in range(0, len(elementsAnhydrite)):
                for j in range(0, len(elementsAnhydrite[i])):
                    resultsAnhydrite[i].append([elementsAnhydrite[i][j][0], round(elementsAnhydrite[i][j][1], 3), round(elementsAnhydrite[i][j][2], 1), round(elementsAnhydrite[i][j][3], 2)])
            #
            results.append(resultsAnhydrite)
            #print("Results Anhydrite:")
            #for j in range(0, len(resultsAnhydrite)):
            #    print(j + 1, ":", resultsAnhydrite[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Anhydrite = open("./outputs/Geochemistry_Anhydrite.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Anhydrite.write(str("rock")+","+str("sample mass")+","+str("C")+","+str("O")+","+str("Na")+","+str("S")+","+str("Cl")+","+str("Ca")+","+str("Fe")+","+str("Cu")+","+str("Mo")+","+str("Pb")+"\n")
            for i in range(0, len(resultsAnhydrite)):
                Geochemistry_Anhydrite.write(str(resultsAnhydrite[i][0])+","+str(resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][7][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][4][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][5][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][3][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][6][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][2][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][10][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][9][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][11][1]/resultsAnhydrite[i][1])+","+str(resultsAnhydrite[i][8][1]/resultsAnhydrite[i][1])+"\n")
        #
        if "ore" in lithologies:
            ## ORE (IRON)
            #
            Mag_ore_Fe = []
            Hem_ore_Fe = []
            Qz_ore_Fe = []
            Cal_ore_Fe = []
            Kfs_ore_Fe = []
            Pl_ore_Fe = []
            for i in range(0, len(data_Mag)):
                if data_Mag[i][0] == "ore":
                    Mag_ore_Fe.append(data_Mag[i][1])
            for i in range(0, len(data_Hem)):
                if data_Hem[i][0] == "ore":
                    Hem_ore_Fe.append(data_Hem[i][1])
            for i in range(0, len(dataQz)):
                if dataQz[i][0] == "ore":
                    Qz_ore_Fe.append(dataQz[i][1])
            for i in range(0, len(dataCal)):
                if dataCal[i][0] == "ore":
                    Cal_ore_Fe.append(dataCal[i][1])
            for i in range(0, len(dataKfs)):
                if dataKfs[i][0] == "ore":
                    Kfs_ore_Fe.append(dataKfs[i][1])
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "ore":
                    Pl_ore_Fe.append(dataPl[i][1])
            for i in range(0, len(ore_Fe)):
                ore_Fe[i].extend(Mag_ore_Fe[i])
                ore_Fe[i].extend(Hem_ore_Fe[i])
                ore_Fe[i].extend(Qz_ore_Fe[i])
                ore_Fe[i].extend(Cal_ore_Fe[i])
                ore_Fe[i].extend(Kfs_ore_Fe[i])
                ore_Fe[i].extend(Pl_ore_Fe[i])
            elements_ore_Fe = []
            for i in range(0, len(ore_Fe)):
                elements_ore_Fe.append([])
            for i in range(0, len(ore_Fe)):
                for j in range(0, len(ore_Fe[i])):
                    if [ore_Fe[i][j][0],0] in elements_ore_Fe[i]:
                        pass
                    else:
                        elements_ore_Fe[i].append([ore_Fe[i][j][0], 0])
            #print("Elements Ore (Fe):", elements_ore_Fe[0])
            for i in range(0, len(ore_Fe)):
                for j in range(0, len(ore_Fe[i])):
                    for k in range(0, len(elements_ore_Fe[i])):
                        if elements_ore_Fe[i][k][0] == ore_Fe[i][j][0]:
                            elements_ore_Fe[i][k][1] += ore_Fe[i][j][1]
                        else:
                            pass
            m_ore_Fe = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "ore":
                    m_ore_Fe.append(round(self.input[i][2],3))
            #print("m:", m_ore_Fe)
            dataM = []
            for i in range(0, len(elements_ore_Fe)):
                m = 0
                for j in range(0, len(elements_ore_Fe[i])):
                    m += round(elements_ore_Fe[i][j][1], 3)
                dataM.append(round(m, 3))
            #print("dataM", dataM)
            for i in range(0, len(elements_ore_Fe)):
                for j in range(0, len(elements_ore_Fe[i])):
                    if elements_ore_Fe[i][j][0] == "C":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemC[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "O":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemO[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "Na":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemNa[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "Al":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemAl[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "Si":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemSi[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "K":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemK[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "Ca":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemCa[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
                    elif elements_ore_Fe[i][j][0] == "Fe":
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/(0.001*chemFe[2]))
                        elements_ore_Fe[i][j].append(elements_ore_Fe[i][j][1]/dataM[i]*100)
            results_ore_Fe = []
            for i in range(0, len(elements_ore_Fe)):
                results_ore_Fe.append(["ore", round(dataM[i],3)])
            for i in range(0, len(elements_ore_Fe)):
                for j in range(0, len(elements_ore_Fe[i])):
                    results_ore_Fe[i].append([elements_ore_Fe[i][j][0], round(elements_ore_Fe[i][j][1],3), round(elements_ore_Fe[i][j][2],1), round(elements_ore_Fe[i][j][3],2)])
            #
            results.append(results_ore_Fe)
            #print("Results Ore (Fe):")
            #for j in range(0, len(results_ore_Fe)):
            #    print(j+1, ":", results_ore_Fe[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Ore_Fe = open("./outputs/Geochemistry_Ore_Fe.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Ore_Fe.write(str("rock")+","+str("sample mass")+","+str("C")+","+str("O")+","+str("Na")+","+str("Al")+","+str("Si")+","+str("K")+","+str("Ca")+","+str("Fe")+"\n")
            for i in range(0, len(results_ore_Fe)):
                Geochemistry_Ore_Fe.write(str(results_ore_Fe[i][0])+","+str(results_ore_Fe[i][1])+","+str(results_ore_Fe[i][6][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][3][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][7][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][9][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][4][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][8][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][5][1]/results_ore_Fe[i][1])+","+str(results_ore_Fe[i][2][1]/results_ore_Fe[i][1])+"\n")
        #
        if "basalt" in lithologies:
            ## Basalt
            #
            basalt_Pl = []
            basalt_En = []
            basalt_Fs = []
            basalt_Bt = []
            basalt_Act = []
            basalt_Tr = []
            basalt_Ol = []
            for i in range(0, len(dataPl)):
                if dataPl[i][0] == "basalt":
                    basalt_Pl.append(dataPl[i][1])
            for i in range(0, len(data_En)):
                if data_En[i][0] == "basalt":
                    basalt_En.append(data_En[i][1])
            for i in range(0, len(data_Fs)):
                if data_Fs[i][0] == "basalt":
                    basalt_Fs.append(data_Fs[i][1])
            for i in range(0, len(data_Bt)):
                if data_Bt[i][0] == "basalt":
                    basalt_Bt.append(data_Bt[i][1])
            for i in range(0, len(dataAct)):
                if dataAct[i][0] == "basalt":
                    basalt_Act.append(dataAct[i][1])
            for i in range(0, len(dataTr)):
                if dataTr[i][0] == "basalt":
                    basalt_Tr.append(dataTr[i][1])
            for i in range(0, len(data_Ol)):
                if data_Ol[i][0] == "basalt":
                    basalt_Ol.append(data_Ol[i][1])
            for i in range(0, len(basalt)):
                basalt[i].extend(basalt_Pl[i])
                basalt[i].extend(basalt_En[i])
                basalt[i].extend(basalt_Fs[i])
                basalt[i].extend(basalt_Bt[i])
                basalt[i].extend(basalt_Act[i])
                basalt[i].extend(basalt_Tr[i])
                basalt[i].extend(basalt_Ol[i])
            elements_basalt = []
            for i in range(0, len(basalt)):
                elements_basalt.append([])
            for i in range(0, len(basalt)):
                for j in range(0, len(basalt[i])):
                    if [basalt[i][j][0],0] in elements_basalt[i]:
                        pass
                    else:
                        elements_basalt[i].append([basalt[i][j][0], 0])
            #print("Elements Basalt:", elements_basalt[0])
            for i in range(0, len(basalt)):
                for j in range(0, len(basalt[i])):
                    for k in range(0, len(elements_basalt[i])):
                        if elements_basalt[i][k][0] == basalt[i][j][0]:
                            elements_basalt[i][k][1] += basalt[i][j][1]
                        else:
                            pass
            m_basalt = []
            for i in range(2, len(self.input)):
                if self.input[i][0] == "basalt":
                    m_basalt.append(round(self.input[i][2],3))
            #print("m:", m_basalt)
            dataM = []
            for i in range(0, len(elements_basalt)):
                m = 0
                for j in range(0, len(elements_basalt[i])):
                    m += round(elements_basalt[i][j][1], 3)
                dataM.append(round(m, 3))
            #print("dataM", dataM)
            for i in range(0, len(elements_basalt)):
                for j in range(0, len(elements_basalt[i])):
                    if elements_basalt[i][j][0] == "H":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemH[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "O":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemO[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "F":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemF[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Na":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemNa[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Mg":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemMg[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Al":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemAl[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Si":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemSi[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "K":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemK[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Ca":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemCa[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Mn":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemMn[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
                    elif elements_basalt[i][j][0] == "Fe":
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/(0.001*chemFe[2]))
                        elements_basalt[i][j].append(elements_basalt[i][j][1]/dataM[i]*100)
            results_basalt = []
            for i in range(0, len(elements_basalt)):
                results_basalt.append(["basalt", round(dataM[i],3)])
            for i in range(0, len(elements_basalt)):
                for j in range(0, len(elements_basalt[i])):
                    results_basalt[i].append([elements_basalt[i][j][0], round(elements_basalt[i][j][1],3), round(elements_basalt[i][j][2],1), round(elements_basalt[i][j][3],2)])
            #
            #print(results_basalt)
            results.append(results_basalt)
            #print("Results Basalt:")
            #for j in range(0, len(results_basalt)):
            #    print(j+1, ":", results_basalt[j])
            #
            # WRITING A CSV FILE
            try:
                Geochemistry_Basalt = open("./outputs/Geochemistry_Basalt.csv","w")
            except:
                print("Error")
                sys.exit(0)
            Geochemistry_Basalt.write(str("rock")+","+str("sample mass")+","+str("H")+","+str("O")+","+str("F")+","+str("Na")+","+str("Mg")+","+str("Al")+","+str("Si")+","+str("K")+","+str("Ca")+","+str("Mn")+","+str("Fe")+"\n")
            for i in range(0, len(results_basalt)):
                Geochemistry_Basalt.write(str(results_basalt[i][0])+","+str(results_basalt[i][1])+","+str(results_basalt[i][9][1]/results_basalt[i][1])+","+str(results_basalt[i][6][1]/results_basalt[i][1])+","+str(results_basalt[i][10][1]/results_basalt[i][1])+","+str(results_basalt[i][2][1]/results_basalt[i][1])+","+str(results_basalt[i][7][1]/results_basalt[i][1])+","+str(results_basalt[i][4][1]/results_basalt[i][1])+","+str(results_basalt[i][5][1]/results_basalt[i][1])+","+str(results_basalt[i][11][1]/results_basalt[i][1])+","+str(results_basalt[i][3][1]/results_basalt[i][1])+","+str(results_basalt[i][12][1]/results_basalt[i][1])+","+str(results_basalt[i][8][1]/results_basalt[i][1])+"\n")
        #
        return results
