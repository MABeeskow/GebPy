#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		geochemistry.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.06.2024

#-----------------------------------------------

## MODULES
import re
import numpy as np
from numpy import round
import scipy.constants as const
import matplotlib.pyplot as plt
from random import *
import random as rd
from modules import minerals
from modules.elements import elements
from modules.chemistry import PeriodicSystem

class MassSpectrometry():
    #
    def __init__(self, data):
        """
        :param data: List/Array of a mineralogical dataset, e.g. quartz
        """
        self.data = data
        MassSpectrometry.show_plot(self)
    #
    def show_plot(self):
        u = const.physical_constants["atomic mass constant"][0]
        na = const.physical_constants["Avogadro constant"][0]
        print(u)
        atoms = self.data[-1]
        print(elements)
        results = []
        for i in range(len(atoms)):
            results.append([atoms[i][0], (PeriodicSystem(name=atoms[i][0]).get_data()[2]*u)/(atoms[i][1]), atoms[i][2]*na])
        print(results)
        results = np.array(results, dtype=object)
        print(results[:, 1])
        fig, ax = plt.subplots(dpi=100)
        ax.scatter(x=results[:, 1], y=results[:, 2])
        plt.show()

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
#
class MineralChemistry():
    #
    def __init__(self, mineral=None, traces=None, w_traces=None, molar_mass_pure=None, majors=None):
        self.mineral = mineral
        self.traces = traces
        self.w_traces = w_traces    # [ [Atom name, amount] ]
        self.molar_mass_pure = molar_mass_pure
        self.majors = majors
        self.compounds = [["H", ["H2O", 2, 1, "O"]], ["Li", ["Li2O", 2, 1, "O"]], ["B", ["B2O3", 2, 3, "O"]],
                          ["Na", ["Na2O", 2, 1, "O"]], ["Mg", ["MgO", 1, 1, "O"]], ["Al", ["Al2O3", 2, 3, "O"]],
                          ["Si", ["SiO2", 1, 2, "O"]], ["K", ["K2O", 2, 1, "O"]], ["Ti", ["TiO2", 1, 2, "O"]],
                          ["Fe", ["Fe2O3", 2, 3, "O"]], ["Cu", ["CuO", 1, 1, "O"]], ["Ga", ["Ga2O3", 2, 3, "O"]],
                          ["Ge", ["GeO2", 1, 2, "O"]], ["As", ["As2O3", 2, 3, "O"]], ["Ag", ["Ag2O", 2, 1, "O"]],
                          ["Sn", ["SnO2", 1, 2, "O"]]]
    #
    def calculate_molar_mass_compounds(self):
        molar_mass_compounds = []
        for item in self.compounds:
            for element in self.traces:
                if element in item:
                    molar_mass_trace_compound = item[1][1]*PeriodicSystem(name=item[0]).get_data()[2] + item[1][2]*\
                                                PeriodicSystem(name=item[1][-1]).get_data()[2]
                    molar_mass_compounds.append([item[0], item[1], molar_mass_trace_compound])
        return molar_mass_compounds
    #
    def calculate_molar_volume(self, volume_cell, z):
        N_A = const.physical_constants["Avogadro constant"][0]
        V_m = (N_A * volume_cell)/(z)
        #
        return V_m
    #
    def calculate_molar_mass(self):
        molar_mass = 0
        molar_masses = []
        if len(self.w_traces) > 0:
            molar_mass_tracer = MineralChemistry(traces=self.w_traces[:, 0]).calculate_molar_mass_compounds()
            for i in range(len(self.w_traces)):
                molar_mass += self.w_traces[i][2]*molar_mass_tracer[i][2]
                molar_masses.append(round(self.w_traces[i][2]*molar_mass_tracer[i][2], 6))
            molar_masses.append(round((1-np.sum(self.w_traces[:, 2]))*self.molar_mass_pure, 6))
            molar_mass += (1-np.sum(self.w_traces[:, 2]))*self.molar_mass_pure
        else:
            molar_mass += self.molar_mass_pure
            for i in range(len(self.majors)):
                molar_masses.append(round(self.majors[i][2]*self.majors[i][3], 6))
        #
        amounts_helper = []
        helper_sum = 0
        amounts = []
        if len(self.w_traces) > 0:
            for i in range(len(self.w_traces)):
                helper_sum += self.w_traces[i][2]*molar_mass
                amounts_helper.append([self.w_traces[i][0], helper_sum])
                amounts.append([self.w_traces[i][0], self.w_traces[i][1], round(self.w_traces[i][2], 6)])
            for i in range(len(self.majors)):
                if self.majors[i][0] != "O":
                    helper_sum += (1-np.sum(self.w_traces[:, 2]))*self.majors[i][3]
                    amounts_helper.append([self.majors[i][0], helper_sum])
                    amounts.append([self.majors[i][0], self.majors[i][1], round((1-np.sum(self.w_traces[:, 2]))*
                                                             self.majors[i][3]/molar_mass, 6)])
                else:
                    pass
            amounts_helper.append(["O", molar_mass-helper_sum])
            amounts.append(["O", 8, round(amounts_helper[-1][1]/molar_mass, 6)])
        else:
            for i in range(len(self.majors)):
                amounts.append([self.majors[i][0], self.majors[i][1], round((self.majors[i][2]*
                                                          self.majors[i][3])/self.molar_mass_pure, 6)])
        amounts = np.array(amounts, dtype=object)
        amounts = amounts[amounts[:, 1].argsort()]
        amounts = amounts.tolist()
        #
        return molar_mass, amounts
#
class Compounds:
    #
    def __init__(self, formula):
        self.formula = formula
        #
    def split_formula(self):
        compound = {}
        chemistry = {}
        for item in self.formula:
            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item)
            if key == None:
                key = re.search(r"([A-Z][a-z]?)(\d*)", item)
            if len(key.groups()) == 4:
                compound[key.group()] = {}
                compound[key.group()]["Total"] = 0
                chemistry[key.group(1)] = PeriodicSystem(name=key.group(1)).get_data()[2]
                if key.group(3) not in chemistry:
                    chemistry[key.group(3)] = PeriodicSystem(name=key.group(3)).get_data()[2]
                #
                if key.group(2) == "":
                    compound[key.group()][key.group(1)] = [1, PeriodicSystem(name=key.group(1)).get_data()[2]]
                    compound[key.group()]["Total"] += round(1*chemistry[key.group(1)], 3)
                else:
                    compound[key.group()][key.group(1)] = [int(key.group(2)), PeriodicSystem(name=key.group(1)).get_data()[2]]
                    compound[key.group()]["Total"] += round(int(key.group(2))*PeriodicSystem(name=key.group(1)).get_data()[2], 3)
                if key.group(4) == "":
                    compound[key.group()][key.group(3)] = [1, PeriodicSystem(name=key.group(3)).get_data()[2]]
                    compound[key.group()]["Total"] += round(1*chemistry[key.group(3)], 3)
                else:
                    compound[key.group()][key.group(3)] = [int(key.group(4)), PeriodicSystem(name=key.group(3)).get_data()[2]]
                    compound[key.group()]["Total"] += round(int(key.group(4))*PeriodicSystem(name=key.group(3)).get_data()[2], 3)
                #
                compound[key.group()]["Total"] = round(compound[key.group()]["Total"], 3)
                compound[key.group()][key.group(1)].append(
                    compound[key.group()][key.group(1)][0]*chemistry[key.group(1)]/compound[key.group()]["Total"])
                compound[key.group()][key.group(3)].append(compound[key.group()][key.group(3)][0]*
                                                           chemistry[key.group(3)]/compound[key.group()]["Total"])
            elif len(key.groups()) == 2:
                compound[key.group()] = {}
                compound[key.group()]["Total"] = 0
                chemistry[key.group(1)] = PeriodicSystem(name=key.group(1)).get_data()[2]
                #
                if key.group(2) == "":
                    compound[key.group()][key.group(1)] = [1, PeriodicSystem(name=key.group(1)).get_data()[2]]
                    compound[key.group()]["Total"] += round(1*chemistry[key.group(1)], 3)
                else:
                    compound[key.group()][key.group(1)] = [int(key.group(2)), PeriodicSystem(name=key.group(1)).get_data()[2]]
                    compound[key.group()]["Total"] += round(int(key.group(2))*PeriodicSystem(name=key.group(1)).get_data()[2], 3)
                #
                compound[key.group()]["Total"] = round(compound[key.group()]["Total"], 3)
                compound[key.group()][key.group(1)].append(round(compound[key.group()][key.group(1)][0]*chemistry[key.group(1)]/compound[key.group()]["Total"], 4))
        #
        return compound
#
class TraceElements:
    #
    def __init__(self, tracer):
        self.traces_data = tracer
    #
    def calculate_composition_oxides(self, var_oxides, var_composition, var_mineral, var_elements):
        ## Split Formulas
        compounds_splitted = Compounds(formula=var_oxides).split_formula()
        #
        cond_w = False
        cond_x = False
        while cond_w == False and cond_x == False:
            w_total = 0
            x_total = 0
            M = 0
            #
            for oxide in var_composition:
                M += var_composition[oxide]*10**(-6)*compounds_splitted[oxide]["Total"]
            element_list = np.sort(var_elements)
            final_comp = {}
            for element in element_list:
                final_comp[element] = {}
                final_comp[element]["w"] = 0
                final_comp[element]["x"] = 0
            for oxide in var_composition:
                for element in element_list:
                    if element in compounds_splitted[oxide]:
                        final_comp[element]["w"] += round(var_composition[oxide]*10**(-6)*compounds_splitted[oxide][element][2], 6)
                        final_comp[element]["x"] += round(var_composition[oxide]*10**(-6)*compounds_splitted[oxide][element][0], 6)
            #
            if var_mineral == "Quartz":
                final_comp["O"]["w"] = 1
                final_comp["O"]["x"] = 2
                final_comp["Si"]["x"] = 1
                for element in final_comp:
                    if element != "O":
                        final_comp["O"]["w"] -= final_comp[element]["w"]
                        if element != "Si":
                            final_comp["Si"]["x"] -= final_comp[element]["x"]
            #
            elif var_mineral == "Magnetite":
                final_comp["O"]["w"] = 1
                final_comp["O"]["x"] = 4
                final_comp["Fe"]["x"] = 3
                for element in final_comp:
                    if element != "O":
                        final_comp["O"]["w"] -= final_comp[element]["w"]
                        if element != "Fe":
                            final_comp["Fe"]["x"] -= final_comp[element]["x"]
            #
            elif var_mineral == "Boehmite":
                final_comp["O"]["w"] = 1
                final_comp["H"]["x"] = 1
                final_comp["O"]["x"] = 2
                final_comp["Al"]["x"] = 1
                for element in final_comp:
                    if element != "O":
                        final_comp["O"]["w"] -= final_comp[element]["w"]
                        if element not in ["Al", "H"]:
                            final_comp["Al"]["x"] -= final_comp[element]["x"]
            #
            for element in final_comp:
                final_comp[element]["w"] = round(final_comp[element]["w"], 6)
                final_comp[element]["x"] = round(final_comp[element]["x"], 6)
                w_total += final_comp[element]["w"]
                x_total += final_comp[element]["x"]
            #
            if np.isclose(w_total, 1.0000) == True:
                cond_w = True
            #
            if var_mineral == "Quartz":
                if np.isclose(x_total, 3.0000) == True:
                    cond_x = True
            elif var_mineral == "Magnetite":
                if np.isclose(x_total, 7.0000) == True:
                    cond_x = True
            elif var_mineral == "Boehmite":
                if np.isclose(x_total, 4.0000) == True:
                    cond_x = True
        #
        return final_comp
    #
    def calculate_composition_sulfides(self, var_elements, var_composition, var_mineral, var_x=None):
        cond_w = False
        cond_x = False
        while cond_w == False and cond_x == False:
            w_total = 0
            x_total = 0
            M = 0
            #
            for element in var_composition:
                M += var_composition[element]*10**(-6)
            element_list = np.sort(var_elements)
            final_comp = {}
            for element in element_list:
                final_comp[element] = {}
                final_comp[element]["w"] = 0
                final_comp[element]["x"] = 0
            for element in element_list:
                final_comp[element]["w"] += round(
                    var_composition[element]*10**(-6), 6)
                final_comp[element]["x"] += round(
                    var_composition[element]*10**(-6), 6)
            #
            if var_mineral == "Cobaltite":
                final_comp["S"]["w"] = 1
                final_comp["S"]["x"] = 1
                final_comp["As"]["x"] = 1
                final_comp["Co"]["x"] = 1
                for element in final_comp:
                    if element != "S":
                        final_comp["S"]["w"] -= final_comp[element]["w"]
                        if element in ["Sb", "Fe"]:
                            final_comp["As"]["x"] -= final_comp[element]["x"]
                        elif element in ["Cu", "Pb", "Ni"]:
                            final_comp["Co"]["x"] -= final_comp[element]["x"]
            #
            if var_mineral == "Marmatite":
                final_comp["S"]["w"] = 1
                final_comp["S"]["x"] = 1
                final_comp["Fe"]["x"] = var_x
                final_comp["Zn"]["x"] = 1 - var_x
                for element in final_comp:
                    if element != "S":
                        final_comp["S"]["w"] -= final_comp[element]["w"]
                        if element in ["Mn", "Cd", "Hg", "In", "Tl", "Ga", "Ge", "Sb", "Sn", "Pb", "Ag", "Co"]:
                            final_comp["Zn"]["x"] -= final_comp[element]["x"]
            #
            for element in final_comp:
                final_comp[element]["w"] = round(final_comp[element]["w"], 6)
                final_comp[element]["x"] = round(final_comp[element]["x"], 6)
                w_total += final_comp[element]["w"]
                x_total += final_comp[element]["x"]
            #
            if np.isclose(w_total, 1.0000) == True:
                cond_w = True
            #
            if var_mineral == "Cobaltite":
                if np.isclose(x_total, 3.0000) == True:
                    cond_x = True
            elif var_mineral == "Marmatite":
                if np.isclose(x_total, 2.0000) == True:
                    cond_x = True
        #
        return final_comp

    def calculate_composition_quartz(self):
        parts_qz = 1000000
        oxides = ["SiO2"]

        cond_w = 0
        cond_x = 0
        while cond_w == False and cond_x == False:
            w_total = 0
            x_total = 0
            #
            element_list = ["O", "Si"]
            composition = {}
            composition[oxides[0]] = parts_qz
            #
            trace_groups = {}
            trace_groups["W"] = []
            trace_groups["X"] = []
            trace_groups["Y"] = []
            trace_groups["Z"] = []
            trace_combinations = {}
            trace_combinations["singles"] = []
            trace_combinations["couples"] = []
            trace_combinations["individual"] = []
            if type(self.traces_data) == dict:
                for tracer in list(self.traces_data.keys()):
                    if tracer in ["Ti", "Ge", "Sn"]:
                        element_list.append(tracer)
                        compound = tracer+str("O2")
                        trace_groups["X"].append(compound)
                        trace_combinations["singles"].append([compound])
                    elif tracer in ["Al", "Fe", "B", "As", "Ga"]:
                        element_list.append(tracer)
                        compound = tracer+str("2O3")
                        trace_groups["Y"].append(compound)
                        trace_combinations["couples"].append([compound])
                    elif tracer in ["H", "Li", "Na", "Ag", "K"]:
                        element_list.append(tracer)
                        compound = tracer+str("2O")
                        trace_groups["Z"].append(compound)
                    elif tracer in ["Mg", "Be", "Ca", "Mn", "Cu", "Zn"]:
                        element_list.append(tracer)
                        compound = tracer+str("O")
                        trace_groups["W"].append(compound)
                        trace_combinations["couples"].append([compound])
                    #
                    trace_combinations["individual"].append([tracer, compound])
            #
            for index_y, couple in enumerate(trace_combinations["couples"], start=0):
                for index_z, item_z in enumerate(trace_groups["Z"], start=0):
                    if item_z not in couple and index_y == index_z:
                        if len(couple) == 2:
                            del couple[-1]
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)
                    elif item_z not in couple and len(couple) == 1 and index_y-index_z != 1:
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)
            #
            for key in trace_combinations:
                if key != "individual":
                    pass
                else:
                    for item in trace_combinations[key]:
                        trace_element = item[0]
                        compound = item[1]
                        oxides.append(compound)
                        val_min = self.traces_data[trace_element]["Min"]
                        val_max = self.traces_data[trace_element]["Max"]
                        mean = (val_min + val_max)/2
                        sigma = (mean - val_min)/3
                        #
                        condition = False
                        while condition == False:
                            amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0])
                            if amount_ppm >= 0:
                                condition = True

                        composition[compound] = amount_ppm
                        composition[oxides[0]] -= amount_ppm
            #
            # print("Groups:", trace_groups)
            # print("Families:", trace_combinations)
            #
            results = Compounds(formula=oxides).split_formula()
            # print("Composition:", composition)
            # print("Results:", results)
            M = 0
            for oxide in composition:
                M += composition[oxide]*10**(-6) * results[oxide]["Total"]
            element_list = np.sort(element_list)
            final_comp = {}
            for element in element_list:
                final_comp[element] = {}
                final_comp[element]["w"] = 0
                final_comp[element]["x"] = 0
            for oxide in composition:
                for element in element_list:
                    if element in results[oxide]:
                        final_comp[element]["w"] += round(composition[oxide]*10**(-6) * results[oxide][element][2], 6)
                        final_comp[element]["x"] += round(composition[oxide]*10**(-6)*results[oxide][element][0], 6)
            final_comp["O"]["w"] = 1
            final_comp["O"]["x"] = 2
            final_comp["Si"]["x"] = 1
            for element in final_comp:
                if element != "O":
                    final_comp["O"]["w"] -= final_comp[element]["w"]
                    if element != "Si":
                        final_comp["Si"]["x"] -= final_comp[element]["x"]
            for element in final_comp:
                final_comp[element]["w"] = round(final_comp[element]["w"], 6)
                final_comp[element]["x"] = round(final_comp[element]["x"], 6)
                w_total += final_comp[element]["w"]
                x_total += final_comp[element]["x"]
            #
            if w_total == 1:
                cond_w = True
            if x_total == 3:
                cond_x = True

        #print("Final:", final_comp)
        #print("Final:", round(w_total, 6), round(x_total, 6))
        #
        return final_comp

    def calculate_composition_orthoclase(self):
        parts_mineral = 1000000
        oxides = ["Al2O3", "SiO2", "K2O"]

        cond_w = 0
        cond_x = 0
        while cond_w == False and cond_x == False:
            w_total = 0
            x_total = 0

            element_list = ["O", "Al", "Si", "K"]
            composition = {}
            composition["Al2O3"] = int(0.1832*parts_mineral)
            composition["SiO2"] = int(0.6476*parts_mineral)
            composition["K2O"] = int(0.1692*parts_mineral)

            trace_groups = {}
            trace_groups["W"] = []
            trace_groups["Y"] = []
            trace_groups["Z"] = []
            trace_combinations = {}
            trace_combinations["singles"] = []
            trace_combinations["couples"] = []
            trace_combinations["individual"] = []

            if type(self.traces_data) == dict:
                for tracer in list(self.traces_data.keys()):
                    if tracer in ["Ba", "Ca", "Fe"]:
                        element_list.append(tracer)
                        compound = tracer+str("O")
                        trace_groups["W"].append(compound)
                        trace_combinations["couples"].append([compound])
                        charge = 2
                    elif tracer in ["Na", "Rb"]:
                        element_list.append(tracer)
                        compound = tracer+str("2O")
                        trace_groups["Z"].append(compound)
                        charge = 1

                    trace_combinations["individual"].append([tracer, compound, charge])

            for index_y, couple in enumerate(trace_combinations["couples"], start=0):
                for index_z, item_z in enumerate(trace_groups["Z"], start=0):
                    if item_z not in couple and index_y == index_z:
                        if len(couple) == 2:
                            del couple[-1]
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)
                    elif item_z not in couple and len(couple) == 1 and index_y-index_z != 1:
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)

            for key in trace_combinations:
                if key != "individual":
                    pass
                else:
                    for item in trace_combinations[key]:
                        trace_element = item[0]
                        compound = item[1]
                        charge = item[2]
                        oxides.append(compound)
                        val_min = self.traces_data[trace_element]["Min"]
                        val_max = self.traces_data[trace_element]["Max"]
                        mean = (val_min + val_max)/2
                        sigma = (mean - val_min)/3

                        condition = False
                        while condition == False:
                            amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0])
                            if amount_ppm >= 0:
                                condition = True

                        composition[compound] = amount_ppm
                        composition["K2O"] -= amount_ppm
            #
            # print("Groups:", trace_groups)
            # print("Families:", trace_combinations)
            #
            results = Compounds(formula=oxides).split_formula()
            #print("Composition:", composition)
            #print("Results:", results)
            M = 0
            for oxide in composition:
                M += composition[oxide]*10**(-6) * results[oxide]["Total"]
            element_list = np.sort(element_list)
            final_comp = {}
            for element in element_list:
                final_comp[element] = {}
                final_comp[element]["w"] = 0
                final_comp[element]["x"] = 0
            for oxide in composition:
                for element in element_list:
                    if element in results[oxide]:
                        final_comp[element]["w"] += round(composition[oxide]*10**(-6) * results[oxide][element][2], 6)
                        final_comp[element]["x"] += round(composition[oxide]*10**(-6)*results[oxide][element][0], 6)
            final_comp["O"]["w"] = 1
            final_comp["O"]["x"] = 8
            final_comp["Al"]["x"] = 1
            final_comp["Si"]["x"] = 3
            final_comp["K"]["x"] = 1
            for element in final_comp:
                if element != "O":
                    final_comp["O"]["w"] -= final_comp[element]["w"]
                    if element not in ["K", "Al", "Si"]:
                        final_comp["K"]["x"] -= final_comp[element]["x"]
            for element in final_comp:
                final_comp[element]["w"] = round(final_comp[element]["w"], 6)
                final_comp[element]["x"] = round(final_comp[element]["x"], 6)
                w_total += final_comp[element]["w"]
                x_total += final_comp[element]["x"]
            #
            if w_total == 1:
                cond_w = True
            if x_total == 13:
                cond_x = True

        #print("Final:", final_comp)
        #print("Final:", round(w_total, 6), round(x_total, 6))
        #
        return final_comp

    def calculate_composition_illite(self):
        parts_mineral = 1000000
        oxides = ["H2O", "Al2O3", "SiO2", "K2O"]

        cond_w = 0
        cond_x = 0
        while cond_w == False and cond_x == False:
            w_total = 0
            x_total = 0

            element_list = ["H", "O", "Al", "Si", "K"]
            composition = {}
            composition["H2O"] = int(0.0468*parts_mineral)
            composition["Al2O3"] = int(0.3509*parts_mineral)
            composition["SiO2"] = int(0.5228*parts_mineral)
            composition["K2O"] = int(0.0795*parts_mineral)

            trace_groups = {}
            trace_groups["W"] = []
            trace_groups["Y"] = []
            trace_groups["Z"] = []
            trace_combinations = {}
            trace_combinations["singles"] = []
            trace_combinations["couples"] = []
            trace_combinations["individual"] = []

            if type(self.traces_data) == dict:
                for tracer in list(self.traces_data.keys()):
                    if tracer in ["Mg", "Fe"]:
                        element_list.append(tracer)
                        compound = tracer+str("O")
                        trace_groups["W"].append(compound)
                        trace_combinations["couples"].append([compound])
                        charge = 2

                    trace_combinations["individual"].append([tracer, compound, charge])

            for index_y, couple in enumerate(trace_combinations["couples"], start=0):
                for index_z, item_z in enumerate(trace_groups["Z"], start=0):
                    if item_z not in couple and index_y == index_z:
                        if len(couple) == 2:
                            del couple[-1]
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)
                    elif item_z not in couple and len(couple) == 1 and index_y-index_z != 1:
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)

            for key in trace_combinations:
                if key != "individual":
                    pass
                else:
                    for item in trace_combinations[key]:
                        trace_element = item[0]
                        compound = item[1]
                        charge = item[2]
                        oxides.append(compound)
                        val_min = self.traces_data[trace_element]["Min"]
                        val_max = self.traces_data[trace_element]["Max"]
                        mean = (val_min + val_max)/2
                        sigma = (mean - val_min)/3

                        condition = False
                        while condition == False:
                            amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0])
                            if amount_ppm >= 0:
                                condition = True

                        composition[compound] = amount_ppm
                        composition["K2O"] -= amount_ppm
            #
            # print("Groups:", trace_groups)
            # print("Families:", trace_combinations)
            #
            results = Compounds(formula=oxides).split_formula()
            #print("Composition:", composition)
            #print("Results:", results)
            M = 0
            for oxide in composition:
                M += composition[oxide]*10**(-6) * results[oxide]["Total"]
            element_list = np.sort(element_list)
            final_comp = {}
            for element in element_list:
                final_comp[element] = {}
                final_comp[element]["w"] = 0
                final_comp[element]["x"] = 0
            for oxide in composition:
                for element in element_list:
                    if element in results[oxide]:
                        final_comp[element]["w"] += round(composition[oxide]*10**(-6) * results[oxide][element][2], 6)
                        final_comp[element]["x"] += round(composition[oxide]*10**(-6)*results[oxide][element][0], 6)
            final_comp["O"]["w"] = 1
            final_comp["O"]["x"] = 12
            final_comp["Al"]["x"] = 2.65
            final_comp["Si"]["x"] = 3.35
            final_comp["K"]["x"] = 0.65
            final_comp["H"]["x"] = 2
            for element in final_comp:
                if element != "O":
                    final_comp["O"]["w"] -= final_comp[element]["w"]
                    if element not in ["K", "Al", "Si", "H"]:
                        final_comp["K"]["x"] -= final_comp[element]["x"]
            for element in final_comp:
                final_comp[element]["w"] = round(final_comp[element]["w"], 6)
                final_comp[element]["x"] = round(final_comp[element]["x"], 6)
                w_total += final_comp[element]["w"]
                x_total += final_comp[element]["x"]
            #
            if w_total == 1:
                cond_w = True
            if x_total == 20.65:
                cond_x = True

        #print("Final:", final_comp)
        #print("Final:", round(w_total, 6), round(x_total, 6))
        #
        return final_comp
    #
    def calculate_composition_apatite_f(self):
        parts_ap = 1000000
        oxides = ["CaO", "PO4"]
        element_list = ["O", "Ca", "P", "F"]
        #
        cond_w = 0
        cond_x = 0
        while cond_w == False and cond_x == False:
            w_total = 0
            x_total = 0
            #
            element_list = ["O", "Ca", "P", "F"]
            composition = {}
            composition[oxides[0]] = parts_ap
            #
            trace_groups = {}
            trace_groups["W"] = []
            trace_groups["X"] = []
            trace_groups["Y"] = []
            trace_groups["Z"] = []
            trace_combinations = {}
            trace_combinations["singles"] = []
            trace_combinations["couples"] = []
            for tracer in self.tracer:
                if tracer in ["Ti", "Zr", "Hf", "Th"]:  # 4+
                    element_list.append(tracer)
                    compound = tracer+str("O2")
                    trace_groups["X"].append(compound)
                    trace_combinations["singles"].append([compound])
                elif tracer in ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Cr", "As"]:   # 3+
                    element_list.append(tracer)
                    compound = tracer+str("2O3")
                    trace_groups["Y"].append(compound)
                    trace_combinations["couples"].append([compound])
                elif tracer in ["Cl", "H", "Rb"]:    # 1+
                    element_list.append(tracer)
                    compound = tracer+str("2O")
                    trace_groups["Z"].append(compound)
                elif tracer in ["Mn", "Co", "Sr", "Ba", "Pb"]:    # 2+
                    element_list.append(tracer)
                    compound = tracer+str("O")
                    trace_groups["W"].append(compound)
                    trace_combinations["couples"].append([compound])
            for index_y, couple in enumerate(trace_combinations["couples"], start=0):
                for index_z, item_z in enumerate(trace_groups["Z"], start=0):
                    if item_z not in couple and index_y == index_z:
                        if len(couple) == 2:
                            del couple[-1]
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)
                    elif item_z not in couple and len(couple) == 1 and index_y-index_z != 1:
                        if couple[0] in trace_groups["W"]:
                            key = re.search(r"([A-Z][a-z]?)(\d*)([A-Z][a-z]?)(\d*)?", item_z)
                            couple.append(key.group(1)+key.group(2))
                        else:
                            couple.append(item_z)
            #
            for key in trace_combinations:
                for item in trace_combinations[key]:
                    amount = rd.uniform(1*10**(-6), 0.1)
                    amount_ppm = int(amount*10000)
                    for compound in item:
                        oxides.append(compound)
                        composition[compound] = amount_ppm
                        composition[oxides[0]] -= amount_ppm
                    item.append(amount_ppm)
            #
            #print("Groups:", trace_groups)
            #print("Families:", trace_combinations)
            #
            results = Compounds(formula=oxides).split_formula()
            #print("Composition:", composition)
            #print("Results:", results)
            M = 0
            for oxide in composition:
                M += composition[oxide]*10**(-6) * results[oxide]["Total"]
            element_list = np.sort(element_list)
            final_comp = {}
            for element in element_list:
                final_comp[element] = {}
                final_comp[element]["w"] = 0
                final_comp[element]["x"] = 0
            for oxide in composition:
                for element in element_list:
                    if element in results[oxide]:
                        #final_comp[element]["w"] += composition[oxide]*10**(-6) * results[oxide][element][0]*results[oxide][element][1]/M
                        final_comp[element]["w"] += round(composition[oxide]*10**(-6) * results[oxide][element][2], 6)
                        final_comp[element]["x"] += round(composition[oxide]*10**(-6)*results[oxide][element][0], 6)
            final_comp["O"]["w"] = 1
            final_comp["O"]["x"] = 2
            final_comp["Ca"]["x"] = 1
            for element in final_comp:
                if element != "O":
                    final_comp["O"]["w"] -= final_comp[element]["w"]
                    if element != "Ca":
                        final_comp["Ca"]["x"] -= final_comp[element]["x"]
            for element in final_comp:
                final_comp[element]["w"] = round(final_comp[element]["w"], 6)
                final_comp[element]["x"] = round(final_comp[element]["x"], 6)
                w_total += final_comp[element]["w"]
                x_total += final_comp[element]["x"]
            #
            if w_total == 1:
                cond_w = True
            if x_total == 3:
                cond_x = True

        #print("Final:", final_comp)
        #print("Final:", round(w_total, 6), round(x_total, 6))
        #
        return final_comp