#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		geochemistry.py
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
from modules.elements import elements

class dataanalysis:
    #
    def __init__(self, lithologies, sequences):
        self.lithologies = lithologies
        self.sequences = sequences
    #
    def sortGeodata(self):
        data = []
        for i in range(0, len(self.lithologies)):
            data.append([self.lithologies[i], [], [], [], [], []])
        #
        for i in range(0, len(self.sequences)):
            for j in range(0, len(self.lithologies)):
                if self.sequences[i][0] == self.lithologies[j]:
                    data[j][1].append(self.sequences[i][4])     # Density
                    data[j][2].append(self.sequences[i][5][0])  # P-wave velocity
                    data[j][3].append(self.sequences[i][5][1])  # S-wave velocity
                    data[j][4].append(self.sequences[i][6])     # Gamma ray
                    data[j][5].append(self.sequences[i][7])     # Neutron porosity
        #
        return data
    #

class geophysics:
    #
    def __init__(self, sequences):
        self.sequences = sequences
    #
    def calculatePressure(self):
        data = [0]
        g = 9.81
        for i in range(0, len(self.sequences)):
            data.append(data[i]+self.sequences[i][4]*1000*g*(self.sequences[i][3]-self.sequences[i][2]))
        #
        return data
    #
    def calculateTemperature(self):
        data = [randint(0, 25)]
        a = 30/1000
        for i in range(0, len(self.sequences)):
            data.append(round(data[0]+self.sequences[i][3]*a, 1))
        #
        return data

class elementanalysis:
    #
    def __init__(self, n, m):
        self.n = n
        self.m = m
    #
    def analyzeQuartz(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        #
        data = []
        m = chemSi[2]+2*chemO[2]
        wSi = chemSi[2]/m
        wO = (2*chemO[2])/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Qz":
                    mSi = wSi*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemSi[0], mSi], [chemO[0], mO]]])
        #
        return data
    #
    def analyzeCalcite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCa = elements.Ca(self)
        chemC = elements.C(self)
        chemO = elements.O(self)
        #
        data = []
        m = chemCa[2]+chemC[2]+3*chemO[2]
        wCa = chemCa[2]/m
        wC = chemC[2]/m
        wO = (3*chemO[2])/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Cal":
                    mCa = wCa*self.m[i][j][1]
                    mC = wC*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemCa[0], mCa], [chemC[0], mC], [chemO[0], mO]]])
        #
        return data
    #
    def analyzeDolomite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCa = elements.Ca(self)
        chemMg = elements.Mg(self)
        chemC = elements.C(self)
        chemO = elements.O(self)
        #
        data = []
        m = chemCa[2]+chemMg[2]+2*(chemC[2]+3*chemO[2])
        wCa = chemCa[2]/m
        wMg = chemMg[2]/m
        wC = 2*chemC[2]/m
        wO = 6*chemO[2]/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Dol":
                    mCa = wCa*self.m[i][j][1]
                    mMg = wMg*self.m[i][j][1]
                    mC = wC*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemCa[0], mCa], [chemMg[0], mMg], [chemC[0], mC], [chemO[0], mO]]])
        #
        return data
    #
    def analyzeOrthoclase(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = elements.K(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        #
        data = []
        m = chemK[2]+chemAl[2]+3*chemSi[2]+8*chemO[2]
        wK = chemK[2]/m
        wAl = chemAl[2]/m
        wSi = 3*chemSi[2]/m
        wO = 8*chemO[2]/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Or":
                    mK = wK*self.m[i][j][1]
                    mAl = wAl*self.m[i][j][1]
                    mSi = wSi*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemK[0], mK], [chemAl[0], mAl], [chemSi[0], mSi], [chemO[0], mO]]])
        #
        return data
    #
    def analyzeBiotite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = elements.K(self)
        chemMg = elements.Mg(self)
        chemFe = elements.Fe(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        chemH = elements.H(self)
        chemF = elements.F(self)
        #
        data = []
        m = chemK[2]+2.5*chemMg[2]+0.5*chemFe[2]+chemAl[2]+3*chemSi[2]+10*chemO[2]+1.75*(chemO[2]+chemH[2])+0.25*chemF[2]
        wK = chemK[2]/m
        wMg = 2.5*chemMg[2]/m
        wFe = 0.5*chemFe[2]/m
        wAl = chemAl[2]/m
        wSi = 3*chemSi[2]/m
        wO = (10*chemO[2]+1.75*chemO[2])/m
        wH = 1.75*chemH[2]/m
        wF = 0.25*chemF[2]/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Bt":
                    mK = wK*self.m[i][j][1]
                    mMg = wMg*self.m[i][j][1]
                    mFe = wFe*self.m[i][j][1]
                    mAl = wAl*self.m[i][j][1]
                    mSi = wSi*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    mH = wH*self.m[i][j][1]
                    mF = wF*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemK[0], mK], [chemMg[0], mMg], [chemFe[0], mFe], [chemAl[0], mAl], [chemSi[0], mSi], [chemO[0], mO], [chemH[0], mH], [chemF[0], mF]]])
        #
        return data
    #
    def analyzeIllite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = elements.K(self)
        chemH = elements.H(self)
        chemO = elements.O(self)
        chemAl = elements.Al(self)
        chemMg = elements.Mg(self)
        chemFe = elements.Fe(self)
        chemSi = elements.Si(self)
        #
        data = []
        m = 0.6*chemK[2]+0.4*(3*chemH[2]+chemO[2])+1.3*chemAl[2]+0.3*chemMg[2]+0.1*chemFe[2]+3.5*chemSi[2]+10*chemO[2]+2*(chemO[2]+chemH[2])+2*chemH[2]+chemO[2]
        wK = (0.6*chemK[2])/m
        wMg = (0.3*chemMg[2])/m
        wAl = (1.3*chemAl[2])/m
        wFe = (0.1*chemFe[2])/m
        wSi = (3.5*chemSi[2])/m
        wH = (0.4*3*chemH[2]+2*2*chemH[2])/m
        wO = (0.4*chemO[2]+10*chemO[2]+2*chemO[2]+chemO[2])/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Ilt":
                    mK = wK*self.m[i][j][1]
                    mMg = wMg*self.m[i][j][1]
                    mAl = wAl*self.m[i][j][1]
                    mFe = wFe*self.m[i][j][1]
                    mSi = wSi*self.m[i][j][1]
                    mH = wH*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemK[0], mK], [chemMg[0], mMg], [chemAl[0], mAl], [chemFe[0], mFe], [chemSi[0], mSi], [chemH[0], mH], [chemO[0], mO]]])
        #
        return data
    #
    def analyzeChlorite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemFe = elements.Fe(self)
        chemMg = elements.Mg(self)
        chemO = elements.O(self)
        chemAl = elements.Al(self)
        chemH = elements.H(self)
        chemSi = elements.Si(self)
        #
        data = []
        m = 3*chemFe[2]+1.5*chemMg[2]+chemAl[2]+0.5*chemFe[2]+3*chemSi[2]+chemAl[2]+12*chemO[2]+6*(chemO[2]+chemH[2])
        wFe = ((3+0.5)*chemFe[2])/m
        wMg = (1.5*chemMg[2])/m
        wAl = (2*chemAl[2])/m
        wSi = (3*chemSi[2])/m
        wO = ((12+6)*chemO[2])/m
        wH = (6*chemH[2])/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Chl":
                    mFe = wFe*self.m[i][j][1]
                    mMg = wMg*self.m[i][j][1]
                    mAl = wAl*self.m[i][j][1]
                    mSi = wSi*self.m[i][j][1]
                    mO = wO*self.m[i][j][1]
                    mH = wH*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemFe[0], mFe], [chemMg[0], mMg], [chemAl[0], mAl], [chemSi[0], mSi], [chemO[0], mO], [chemH[0], mH]]])
        #
        return data
    #
    def analyzeGlauconite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = elements.K(self)
        chemNa = elements.Na(self)
        chemMg = elements.Mg(self)
        chemFe = elements.Fe(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        chemH = elements.H(self)
        #
        data = []
        m = 0.6*chemK[2]+0.05*chemNa[2]+1.5*chemFe[2]+0.4*chemMg[2]+0.3*chemAl[2]+3.8*chemSi[2]+10*chemO[2]+2*(chemO[2]+chemH[2])
        wK = (0.6*chemK[2])/m
        wNa = (0.05*chemNa[2])/m
        wMg = (0.4*chemMg[2])/m
        wAl = (0.3*chemAl[2])/m
        wFe = ((1.3+0.2)*chemFe[2])/m
        wSi = (3.8*chemSi[2])/m
        wO = ((10+2)*chemO[2])/m
        wH = (2*chemH[2])/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Glt":
                    dataK = wK*self.m[i][j][1]
                    dataNa = wNa*self.m[i][j][1]
                    dataMg = wMg*self.m[i][j][1]
                    dataAl = wAl*self.m[i][j][1]
                    dataFe = wFe*self.m[i][j][1]
                    dataSi = wSi*self.m[i][j][1]
                    dataO = wO*self.m[i][j][1]
                    dataH = wH*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemK[0], dataK], [chemNa[0], dataNa], [chemMg[0], dataMg], [chemAl[0], dataAl], [chemFe[0], dataFe], [chemSi[0], dataSi], [chemO[0], dataO], [chemH[0], dataH]]])
        #
        return data
    #
    def analyzeAlbite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemNa = elements.K(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        #
        data = []
        m = chemNa[2] + chemAl[2] + 3*chemSi[2] + 8*chemO[2]
        wNa = (chemNa[2])/m
        wAl = (chemAl[2])/m
        wSi = (3*chemSi[2])/m
        wO = (8*chemO[2])/m
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Ab":
                    dataNa = wNa*self.m[i][j][1]
                    dataAl = wAl*self.m[i][j][1]
                    dataSi = wSi*self.m[i][j][1]
                    dataO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemNa[0], dataNa], [chemAl[0], dataAl], [chemSi[0], dataSi], [chemO[0], dataO]]])
        #
        return data
    #
    def analyzeAnorthite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCa = elements.K(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        #
        data = []
        M = chemCa[2] + 2*chemAl[2] + 2*chemSi[2] + 8*chemO[2]
        wCa = (chemCa[2])/M
        wAl = (2*chemAl[2])/M
        wSi = (2*chemSi[2])/M
        wO = (8*chemO[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "An":
                    dataCa = wCa*self.m[i][j][1]
                    dataAl = wAl*self.m[i][j][1]
                    dataSi = wSi*self.m[i][j][1]
                    dataO = wO*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemCa[0], dataCa], [chemAl[0], dataAl], [chemSi[0], dataSi], [chemO[0], dataO]]])
        #
        return data
    #
    def analyzeMuscovite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = elements.K(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        chemH = elements.H(self)
        chemF = elements.F(self)
        #
        data = []
        M = chemK[2]+3*chemAl[2]+3*chemSi[2]+10*chemO[2]+1.8*(chemO[2]+chemH[2])+0.2*chemF[2]
        wK = (chemK[2])/M
        wAl = (3*chemAl[2])/M
        wSi = (3*chemSi[2])/M
        wO = ((10+1.8)*chemO[2])/M
        wH = (1.8*chemH[2])/M
        wF = (0.2*chemF[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Ms":
                    dataK = wK*self.m[i][j][1]
                    dataAl = wAl*self.m[i][j][1]
                    dataSi = wSi*self.m[i][j][1]
                    dataO = wO*self.m[i][j][1]
                    dataH = wH*self.m[i][j][1]
                    dataF = wF*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemK[0], dataK], [chemAl[0], dataAl], [chemSi[0], dataSi], [chemO[0], dataO], [chemH[0], dataH], [chemF[0], dataF]]])
        #
        return data
    #
    def analyzePyrite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemFe = elements.Fe(self)
        chemS = elements.S(self)
        #
        data = []
        M = chemFe[2] + 2*chemS[2]
        wFe = (chemFe[2])/M
        wS = (2*chemS[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Py":
                    dataFe = wFe*self.m[i][j][1]
                    dataS = wS*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemFe[0], dataFe], [chemS[0], dataS]]])
        #
        return data
    #
    def analyzeHalite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemNa = elements.Na(self)
        chemCl = elements.Cl(self)
        #
        data = []
        M = chemNa[2] + chemCl[2]
        wNa = (chemNa[2])/M
        wCl = (chemCl[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Hl":
                    dataNa = wNa*self.m[i][j][1]
                    dataCl = wCl*self.m[i][j][1]
                    data.append([self.m[i][0], [[chemNa[0], dataNa], [chemCl[0], dataCl]]])
        #
        return data
        #
    def analyzeAnhydrite(self):
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCa = elements.Ca(self)
        chemS = elements.S(self)
        chemO = elements.O(self)
        #
        data = []
        M = chemCa[2] + chemS[2] + 4*chemO[2]
        wCa = (chemCa[2]) / M
        wS = (chemS[2]) / M
        wO = (chemO[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Anh":
                    dataCa = wCa * self.m[i][j][1]
                    dataS = wS * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemCa[0], dataCa], [chemS[0], dataS], [chemO[0], dataO]]])
        #
        return data
    #
    def analyzeGypsum(self): # CaSO4*H2O
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCa = elements.Ca(self)
        chemS = elements.S(self)
        chemO = elements.O(self)
        chemH = elements.H(self)
        #
        data = []
        M = chemCa[2] + chemS[2] + 4*chemO[2] + 2*chemH[2] + chemO[2]
        wCa = (chemCa[2]) / M
        wS = (chemS[2]) / M
        wO = (chemO[2]) / M
        wH = (chemH[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Gp":
                    dataCa = wCa * self.m[i][j][1]
                    dataS = wS * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    dataH = wH * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemCa[0], dataCa], [chemS[0], dataS], [chemO[0], dataO], [chemH[0], dataH]]])
        #
        return data
    #
    def analyzeSylvite(self):  # KCl
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = elements.K(self)
        chemCl = elements.Cl(self)
        #
        data = []
        M = chemK[2] + chemCl[2]
        wK = (chemK[2]) / M
        wCl = (chemCl[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Syl":
                    dataK = wK * self.m[i][j][1]
                    dataCl = wCl * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemK[0], dataK], [chemCl[0], dataCl]]])
        return data
        #
    def analyzeGalena(self):  # PbS
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemPb = elements.Pb(self)
        chemS = elements.S(self)
        #
        data = []
        M = chemPb[2] + chemS[2]
        wPb = (chemPb[2]) / M
        wS = (chemS[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Gn":
                    dataPb = wPb * self.m[i][j][1]
                    dataS = wS * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemPb[0], dataPb], [chemS[0], dataS]]])
        return data
    #
    def analyzeChalcopyrite(self):  # CuFeS2
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCu = elements.Cu(self)
        chemFe = elements.Fe(self)
        chemS = elements.S(self)
        #
        data = []
        M = chemCu[2] + chemFe[2] + 2*chemS[2]
        wCu = (chemCu[2]) / M
        wFe = (chemFe[2]) / M
        wS = (chemS[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Ccp":
                    dataCu = wCu * self.m[i][j][1]
                    dataFe = wFe * self.m[i][j][1]
                    dataS = wS * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemCu[0], dataCu], [chemFe[0], dataFe], [chemS[0], dataS]]])
        return data
    #
    def analyzeMolybdenite(self):  # MoS2
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMo = elements.Mo(self)
        chemS = elements.S(self)
        #
        data = []
        M = chemMo[2] + 2*chemS[2]
        wMo = (chemMo[2]) / M
        wS = (chemS[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Mol":
                    dataMo = wMo * self.m[i][j][1]
                    dataS = wS * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemMo[0], dataMo], [chemS[0], dataS]]])
        return data
    #
    def analyzeKaolinite(self):  # Al2(OH)4Si2O5
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemAl = elements.Al(self)
        chemO = elements.O(self)
        chemH = elements.H(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        #
        data = []
        M = round(2*chemAl[2] + 4*(chemO[2]+chemH[2]) + 2*chemSi[2] + 5*chemO[2], 3)
        wAl = (chemAl[2]) / M
        wO = (chemO[2]) / M
        wH = (chemH[2]) / M
        wSi = (chemSi[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Kln":
                    dataAl = wAl * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    dataH = wH * self.m[i][j][1]
                    dataSi = wSi * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemAl[0], dataAl], [chemO[0], dataO], [chemH[0], dataH], [chemSi[0], dataSi]]])
        return data
    #
    def analyzeAlkalifeldspar(self, sequences):  # Na(x)K(1-x)AlSi3O8
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemNa = elements.Na(self)
        chemK = elements.K(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        self.sequences = sequences
        #
        data = []
        for i in range(0, len(self.sequences)):
            for j in range(0, len(self.sequences[i][9])):
                if self.sequences[i][9][j][0] == "Kfs" or self.sequences[i][9][j][0] == "Afs":
                    M = self.sequences[i][9][j][2]
        wNa = (chemNa[2]) / M
        wK = (chemK[2]) / M
        wAl = (chemAl[2]) / M
        wSi = (chemSi[2]) / M
        wO = (chemO[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Kfs":
                    dataNa = wNa * self.m[i][j][1]
                    dataK = wK * self.m[i][j][1]
                    dataAl = wAl * self.m[i][j][1]
                    dataSi = wSi * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemNa[0], dataNa], [chemK[0], dataK], [chemAl[0], dataAl], [chemSi[0], dataSi], [chemO[0], dataO]]])
        return data
    #
    def analyzePlagioclase(self, sequences):  # Na(x)Ca(1-x)Al(2-x)Si(2+x)O8
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemNa = elements.Na(self)
        chemCa = elements.Ca(self)
        chemAl = elements.Al(self)
        chemSi = elements.Si(self)
        chemO = elements.O(self)
        self.sequences = sequences
        #
        data = []
        for i in range(0, len(self.sequences)):
            for j in range(0, len(self.sequences[i][9])):
                if self.sequences[i][9][j][0] == "Pl":
                    M = self.sequences[i][9][j][2]
        wNa = (chemNa[2]) / M
        wCa = (chemCa[2]) / M
        wAl = (chemAl[2]) / M
        wSi = (chemSi[2]) / M
        wO = (chemO[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Pl":
                    dataNa = wNa * self.m[i][j][1]
                    dataCa = wCa * self.m[i][j][1]
                    dataAl = wAl * self.m[i][j][1]
                    dataSi = wSi * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemNa[0], dataNa], [chemCa[0], dataCa], [chemAl[0], dataAl], [chemSi[0], dataSi], [chemO[0], dataO]]])
        return data
    #
    def analyzeSiderite(self):  # FeCO3
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemFe = elements.Fe(self)
        chemC = elements.C(self)
        chemO = elements.O(self)
        #
        data = []
        M = round(chemFe[2] + chemC[2] + 3*chemO[2], 3)
        wFe = (chemFe[2]) / M
        wC = (chemC[2]) / M
        wO = (chemO[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Sd":
                    dataFe = wFe * self.m[i][j][1]
                    dataC = wC * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemFe[0], dataFe], [chemC[0], dataC], [chemO[0], dataO]]])
        return data
    #
    def analyzeMagnesite(self):  # MgCO3
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMg = elements.Mg(self)
        chemC = elements.C(self)
        chemO = elements.O(self)
        #
        data = []
        M = round(chemMg[2] + chemC[2] + 3*chemO[2], 3)
        wMg = (chemMg[2]) / M
        wC = (chemC[2]) / M
        wO = (chemO[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Mgs":
                    dataMg = wMg * self.m[i][j][1]
                    dataC = wC * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    data.append([self.m[i][0], [[chemMg[0], dataMg], [chemC[0], dataC], [chemO[0], dataO]]])
        return data
    #
    def analyzeActinolite(self, sequences):  # Ca2(Mg,Fe)5Si8O22(OH)2
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Ca = elements.Ca(self)
        chem_Mg = elements.Mg(self)
        chem_Fe = elements.Fe(self)
        chem_Si = elements.Si(self)
        chem_O = elements.O(self)
        chem_H = elements.H(self)
        self.sequences = sequences
        #
        data = []
        for i in range(0, len(self.sequences)):
            for j in range(0, len(self.sequences[i][9])):
                if self.sequences[i][9][j][0] == "Act":
                    M = round(2*chem_Ca[2] + 5*(self.sequences[i][9][j][2]*chem_Mg[2]+(1-self.sequences[i][9][j][2])*chem_Fe[2]) + 8*chem_Si[2] + 22*chem_O[2] + 2*(chem_O[2]+chem_H[2]), 3)
        wCa = (chem_Ca[2]) / M
        wMg = (chem_Mg[2]) / M
        wFe = (chem_Fe[2]) / M
        wSi = (chem_Si[2]) / M
        wO = (chem_O[2]) / M
        wH = (chem_H[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Act":
                    dataCa = wCa * self.m[i][j][1]
                    dataMg = wMg * self.m[i][j][1]
                    dataFe = wFe * self.m[i][j][1]
                    dataSi = wSi * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    dataH = wH * self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_Ca[0], dataCa], [chem_Mg[0], dataMg], [chem_Fe[0], dataFe], [chem_Si[0], dataSi], [chem_O[0], dataO], [chem_H[0], dataH]]])
        return data
    #
    def analyzeTremolite(self):  # Ca2Mg5Si8O22(OH)2
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Ca = elements.Ca(self)
        chem_Mg = elements.Mg(self)
        chem_Si = elements.Si(self)
        chem_O = elements.O(self)
        chem_H = elements.H(self)
        #
        data = []
        M = round(2*chem_Ca[2] + 5*chem_Mg[2] + 2*(4*chem_Si[2]+11*chem_O[2]+(chem_O[2]+chem_H[2])), 3)
        wCa = (chem_Ca[2]) / M
        wMg = (chem_Mg[2]) / M
        wSi = (chem_Si[2]) / M
        wO = (chem_O[2]) / M
        wH = (chem_H[2]) / M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Tr":
                    dataCa = wCa * self.m[i][j][1]
                    dataMg = wMg * self.m[i][j][1]
                    dataSi = wSi * self.m[i][j][1]
                    dataO = wO * self.m[i][j][1]
                    dataH = wH * self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_Ca[0], dataCa], [chem_Mg[0], dataMg], [chem_Si[0], dataSi], [chem_O[0], dataO], [chem_H[0], dataH]]])
        return data
    #
    def analyzeMagnetite(self):  # Fe3O4
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Fe = elements.Fe(self)
        chem_O = elements.O(self)
        #
        data = []
        M = round(3*chem_Fe[2] + 4*chem_O[2], 3)
        w_Fe = (chem_Fe[2])/M
        w_O = (chem_O[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Mag":
                    data_Fe = w_Fe*self.m[i][j][1]
                    data_O = w_O*self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_Fe[0], data_Fe], [chem_O[0], data_O]]])
        return data
    #
    def analyzeHematite(self):  # Fe2O3
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Fe = elements.Fe(self)
        chem_O = elements.O(self)
        #
        data = []
        M = round(2*chem_Fe[2] + 3*chem_O[2], 3)
        w_Fe = (chem_Fe[2])/M
        w_O = (chem_O[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Hem":
                    data_Fe = w_Fe*self.m[i][j][1]
                    data_O = w_O*self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_Fe[0], data_Fe], [chem_O[0], data_O]]])
        return data
    #
    def analyze_biotite(self, sequences):  # K [Mg(a)Fe(1-a)]3 Al(b)Fe(1-b) [Al(c)Si(1-c)]3 O10 [OH(d)F(1-d)]2
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_H = elements.H(self)
        chem_O = elements.O(self)
        chem_F = elements.F(self)
        chem_Mg = elements.Mg(self)
        chem_Al = elements.Al(self)
        chem_Si = elements.Si(self)
        chem_K = elements.K(self)
        chem_Fe = elements.Fe(self)
        self.sequences = sequences
        #
        data = []
        for i in range(0, len(self.sequences)):
            for j in range(0, len(self.sequences[i][9])):
                if self.sequences[i][9][j][0] == "Bt":
                    M = self.sequences[i][9][j][2]
        w_H = (chem_H[2])/M
        w_O = (chem_O[2])/M
        w_F = (chem_F[2])/M
        w_Mg = (chem_Mg[2])/M
        w_Al = (chem_Al[2])/M
        w_Si = (chem_Si[2])/M
        w_K = (chem_K[2])/M
        w_Fe = (chem_Fe[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Bt":
                    data_H = w_H*self.m[i][j][1]
                    data_O = w_O*self.m[i][j][1]
                    data_F = w_F*self.m[i][j][1]
                    data_Mg = w_Mg*self.m[i][j][1]
                    data_Al = w_Al*self.m[i][j][1]
                    data_Si = w_Si*self.m[i][j][1]
                    data_K = w_K*self.m[i][j][1]
                    data_Fe = w_Fe*self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_H[0], data_H], [chem_O[0], data_O], [chem_F[0], data_F], [chem_Mg[0], data_Mg], [chem_Al[0], data_Al], [chem_Si[0], data_Si], [chem_K[0], data_K], [chem_Fe[0], data_Fe]]])
        return data
    #
    def analyze_enstatite(self):  # Mg2Si2O6
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Mg = elements.Mg(self)
        chem_Si = elements.Si(self)
        chem_O = elements.O(self)
        #
        data = []
        M = round(2*chem_Mg[2] + 2*chem_Si[2] + 6*chem_O[2],3)
        w_Mg = (chem_Mg[2])/M
        w_Si = (chem_Si[2])/M
        w_O = (chem_O[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "En":
                    data_Mg = w_Mg*self.m[i][j][1]
                    data_Si = w_Si*self.m[i][j][1]
                    data_O = w_O*self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_Mg[0], data_Mg], [chem_Si[0], data_Si], [chem_O[0], data_O]]])
        return data
    #
    def analyze_ferrosilite(self):  # FeSiO3
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_Fe = elements.Fe(self)
        chem_Si = elements.Si(self)
        chem_O = elements.O(self)
        #
        data = []
        M = round(chem_Fe[2] + chem_Si[2] + 3*chem_O[2], 3)
        w_Fe = (chem_Fe[2])/M
        w_Si = (chem_Si[2])/M
        w_O = (chem_O[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Fs":
                    data_Fe = w_Fe*self.m[i][j][1]
                    data_Si = w_Si*self.m[i][j][1]
                    data_O = w_O*self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_Fe[0], data_Fe], [chem_Si[0], data_Si], [chem_O[0], data_O]]])
        return data
    #
    def analyze_olivine(self, sequences):  # (Mg,Fe,Mn)2 Si O4
        # input = n
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chem_O = elements.O(self)
        chem_Mg = elements.Mg(self)
        chem_Si = elements.Si(self)
        chem_Mn = elements.Mn(self)
        chem_Fe = elements.Fe(self)
        self.sequences = sequences
        #
        data = []
        for i in range(0, len(self.sequences)):
            for j in range(0, len(self.sequences[i][9])):
                if self.sequences[i][9][j][0] == "Ol":
                    M = self.sequences[i][9][j][2]
        w_O = (chem_O[2])/M
        w_Mg = (chem_Mg[2])/M
        w_Si = (chem_Si[2])/M
        w_Mn = (chem_Mn[2])/M
        w_Fe = (chem_Fe[2])/M
        #
        for i in range(0, len(self.m)):
            for j in range(2, len(self.m[i])):
                if self.m[i][j][0] == "Ol":
                    data_O = w_O*self.m[i][j][1]
                    data_Mg = w_Mg*self.m[i][j][1]
                    data_Si = w_Si*self.m[i][j][1]
                    data_Mn = w_Mn*self.m[i][j][1]
                    data_Fe = w_Fe*self.m[i][j][1]
                    data.append([self.m[i][0], [[chem_O[0], data_O], [chem_Mg[0], data_Mg], [chem_Si[0], data_Si], [chem_Mn[0], data_Mn], [chem_Fe[0], data_Fe]]])
        return data
#
class Fractions:
    #
    def __init__(self,):
        pass
    #
    def calculate_volume_fraction(self, mineralogy, w):
        self.mineralogy = mineralogy
        self.w = w
        #
        V0 = np.sum([self.w[i]/self.mineralogy[i][2] for i in range(len(self.mineralogy))])
        #
        phi = []
        for i in range(len(self.mineralogy)):
            phi.append(self.w[i]/(self.mineralogy[i][2]*V0))
        #
        return phi