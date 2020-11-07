#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		Seismology_Single-Wiggle.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		05.11.2020

# -----------------------------------------------

## MODULES
import sys
sys.path.append('../modules')
import numpy as np
from numpy import sqrt, exp, sin, cos, arcsin, arccos, pi
from random import *
import random as rd
from scipy import constants
import  matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from scipy import interpolate
import pywt
#import wiggle
from modules import carbonates
from modules import siliciclastics
from modules.sequences import surface, sand, subsurface
from modules import elements
from modules import minerals
from modules.geochemistry import dataanalysis, elementanalysis
from modules.geophysics import geophysics, seismics
from modules import core
from modules.core import sample
from modules import plotting

## CONSTANTS / PROPERTIES
pi = constants.pi
avogadro = constants.Avogadro
u = constants.atomic_mass
# Example = [name, [density], [p-wave velocity], [s-wave velocity], [gamma ray], [neutron]]
shale = ["shale", [2.0, 2.4, np.random.normal((1.8 + 2.8) / 2, 0.1)], [1100, 2500, np.mean([1100, 2500])],
         [200, 800, np.mean([200, 800])], [90, 450, np.mean([90, 450])], [25, 75, np.mean([25, 75])],
         [0.0, 0.1, np.random.normal((0.0 + 0.1) / 2, 0.1)]]
#
mineralEnstatite = minerals.inosilicates.enstatite("")


#
def createRicker(f, t, t0, sigma):
    #
    # a = ((t-t0)**2)/(sigma**2)
    a = (np.pi ** 2) * (f ** 2) * ((t - t0) ** 2)
    x0 = (1 - a) * exp(-a / 2)
    #
    return x0


# sampleFreq = 500
# t = np.arange(0, 0.1, 1/sampleFreq)
# x  = createRicker(100, t, 1/100, pi)
# lines = plt.plot(t, x, "-x")
# plt.show()

def drawWiggle(ax, x, t, xOffset):
    #
    wiggle = ax.plot(x + xOffset, t, color="black")
    if len(t) > 1:
        tracefill = np.array(x)
        tracefill[0] = 0.0
        tracefill[-1] = 0.0
        tracefill[np.nonzero(x > x[0])] = 0
        fill, = ax.fill(tracefill + xOffset, t, color="white")
        tracefill = np.array(x)
        tracefill[0] = 0.0
        tracefill[-1] = 0.0
        tracefill[np.nonzero(x > x[0])] = 0
        fill, = ax.fill(tracefill + xOffset, t, color="black")


#
# fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
# drawWiggle(ax, x, t, 0)
# drawWiggle(ax, x, t, 1)
# drawWiggle(ax, x, t, 2)
# ax.invert_yaxis()
# plt.show()

## CREATE SEQUENCE STRATIGRAPHY
#
seqCounter = 0
actualThickness = 0
maxThickness = 1500
parts = 20
sequences = []
initialsequence = randint(0, 2)
#
# sigmaShale = [0.2, 350, 60, 1]
#
if initialsequence == 0:
    data = surface("water", actualThickness)
    water = data.createSurface()
    sequences.append(water)
    actualThickness += water[1]
    #
    data = sand(actualThickness)
    wetSand = data.createWetSand()
    sequences.append(wetSand)
    actualThickness += wetSand[1]
    seqCounter += 1
elif initialsequence == 1:
    data = surface("dry sand", actualThickness)
    drySand = data.createSurface()
    sequences.append(drySand)
    actualThickness += drySand[1]
    #
    data = sand(actualThickness)
    wetSand = data.createWetSand()
    sequences.append(wetSand)
    actualThickness += wetSand[1]
    seqCounter += 1
elif initialsequence == 2:
    data = surface("soil", actualThickness)
    Soil = data.createSurface()
    sequences.append(Soil)
    actualThickness += Soil[1]
    #
    data = sand(actualThickness)
    wetSand = data.createWetSand()
    sequences.append(wetSand)
    actualThickness += wetSand[1]
    seqCounter += 1
while actualThickness < maxThickness:
    #
    ##################
    # AFTER WET SAND #
    ##################
    if sequences[seqCounter][0] == "wet sand":
        magicnumber = rd.randint(0, 100)
        if magicnumber in range(0,25):
            #############
            # SANDSTONE #
            #############
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataSandstone = data.createSiliciclastics()
            for i in range(0, len(dataSandstone)):
                sequences.append(dataSandstone[i])
                actualThickness += dataSandstone[i][1]
                seqCounter += 1
        elif magicnumber in range(25,50):
            #############
            # CARBONATE #
            #############
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataCarbonate = data.createCarbonates()
            for i in range(0, len(dataCarbonate)):
                sequences.append(dataCarbonate[i])
                actualThickness += dataCarbonate[i][1]
                seqCounter += 1
        elif magicnumber in range(50,75):
            #########
            # SHALE #
            #########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataShale = data.createShales()
            for i in range(0, len(dataShale)):
                sequences.append(dataShale[i])
                actualThickness += dataShale[i][1]
                seqCounter += 1
        elif magicnumber in range(75,90):
            #############
            # ANHYDRITE #
            #############
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataAnhydrite = data.createEvaporites("Anhydrite")
            for i in range(0, len(dataAnhydrite)):
                sequences.append(dataAnhydrite[i])
                actualThickness += dataAnhydrite[i][1]
                seqCounter += 1
        elif magicnumber in range(90,100):
            ##############
            # ORE (IRON) #
            ##############
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            data_ore_Fe = data.create_ores()
            for i in range(0, len(data_ore_Fe)):
                sequences.append(data_ore_Fe[i])
                actualThickness += data_ore_Fe[i][1]
                seqCounter += 1
    #
    ###################
    # AFTER SANDSTONE #
    ###################
    elif sequences[seqCounter][0] == "sandstone":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(0,50):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                for i in range(0, len(dataSandstone)):
                    sequences.append(dataSandstone[i])
                    actualThickness += dataSandstone[i][1]
                    seqCounter += 1
            elif magicnumber in range(50,65):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates()
                for i in range(0, len(dataCarbonate)):
                    sequences.append(dataCarbonate[i])
                    actualThickness += dataCarbonate[i][1]
                    seqCounter += 1
            elif magicnumber in range(65,80):
                #########
                # SHALE #
                #########
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataShale = data.createShales()
                for i in range(0, len(dataShale)):
                    sequences.append(dataShale[i])
                    actualThickness += dataShale[i][1]
                    seqCounter += 1
            elif magicnumber in range(80,95):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite", "Anhydrite"]
                magicnumber = rd.randint(0, 1)
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(95,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ###################
    # AFTER LIMESTONE #
    ###################
    elif sequences[seqCounter][0] == "limestone":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(50,65):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                for i in range(0, len(dataSandstone)):
                    sequences.append(dataSandstone[i])
                    actualThickness += dataSandstone[i][1]
                    seqCounter += 1
            elif magicnumber in range(0,33):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates(keyword="limestone")
                for i in range(0, len(dataCarbonate)):
                    sequences.append(dataCarbonate[i])
                    actualThickness += dataCarbonate[i][1]
                    seqCounter += 1
            elif magicnumber in range(33,50):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates(keyword="dolomite")
                for i in range(0, len(dataCarbonate)):
                    sequences.append(dataCarbonate[i])
                    actualThickness += dataCarbonate[i][1]
                    seqCounter += 1
            elif magicnumber in range(65,80):
                #########
                # SHALE #
                #########
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataShale = data.createShales()
                for i in range(0, len(dataShale)):
                    sequences.append(dataShale[i])
                    actualThickness += dataShale[i][1]
                    seqCounter += 1
            elif magicnumber in range(90,95):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite", "Anhydrite"]
                magicnumber = rd.randint(0, 1)
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(95,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ##################
    # AFTER DOLOMITE #
    ##################
    elif sequences[seqCounter][0] == "dolomite":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(50,65):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                for i in range(0, len(dataSandstone)):
                    sequences.append(dataSandstone[i])
                    actualThickness += dataSandstone[i][1]
                    seqCounter += 1
            elif magicnumber in range(33,50):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates(keyword="limestone")
                for i in range(0, len(dataCarbonate)):
                    sequences.append(dataCarbonate[i])
                    actualThickness += dataCarbonate[i][1]
                    seqCounter += 1
            elif magicnumber in range(0,33):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates(keyword="dolomite")
                for i in range(0, len(dataCarbonate)):
                    sequences.append(dataCarbonate[i])
                    actualThickness += dataCarbonate[i][1]
                    seqCounter += 1
            elif magicnumber in range(65,80):
                #########
                # SHALE #
                #########
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataShale = data.createShales()
                for i in range(0, len(dataShale)):
                    sequences.append(dataShale[i])
                    actualThickness += dataShale[i][1]
                    seqCounter += 1
            elif magicnumber in range(90,95):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite", "Anhydrite"]
                magicnumber = rd.randint(0, 1)
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(95,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ###############
    # AFTER SHALE #
    ###############
    elif sequences[seqCounter][0] == "shale":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(0,35):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                nData = len(dataSandstone)
                if dataSandstone[0][8] == "water" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "oil" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "gas" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
            elif magicnumber in range(35,70):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates()
                nData = len(dataCarbonate)
                if dataCarbonate[0][8] == "water" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "oil" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "gas" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
            elif magicnumber in range(70,90):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite", "Anhydrite"]
                magicnumber = rd.randint(0, 1)
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(90,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ################
    # AFTER HALITE #
    ################
    elif sequences[seqCounter][0] == "halite":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(0,30):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                nData = len(dataSandstone)
                if dataSandstone[0][8] == "water" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "oil" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "gas" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
            elif magicnumber in range(30,60):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates()
                nData = len(dataCarbonate)
                if dataCarbonate[0][8] == "water" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "oil" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "gas" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
            elif magicnumber in range(60,90):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite"]
                magicnumber = 0
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(90,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ###################
    # AFTER ANHYDRITE #
    ###################
    elif sequences[seqCounter][0] == "anhydrite":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(0,30):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                nData = len(dataSandstone)
                if dataSandstone[0][8] == "water" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "oil" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "gas" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
            elif magicnumber in range(30,60):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates()
                nData = len(dataCarbonate)
                if dataCarbonate[0][8] == "water" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "oil" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "gas" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
            elif magicnumber in range(60,90):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite"]
                magicnumber = 0
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(90,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ####################
    # AFTER ORE (IRON) #
    ####################
    elif sequences[seqCounter][0] == "ore":
        if actualThickness < 0.90 * maxThickness:
            magicnumber = rd.randint(0, 100)
            if magicnumber in range(0,30):
                #############
                # SANDSTONE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataSandstone = data.createSiliciclastics()
                nData = len(dataSandstone)
                if dataSandstone[0][8] == "water" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "oil" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
                elif dataSandstone[0][8] == "gas" and len(dataSandstone) > 0:
                    for i in range(0, len(dataSandstone)):
                        sequences.append(dataSandstone[i])
                        actualThickness += dataSandstone[i][1]
                        seqCounter += 1
            elif magicnumber in range(30,60):
                #############
                # CARBONATE #
                #############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                dataCarbonate = data.createCarbonates()
                nData = len(dataCarbonate)
                if dataCarbonate[0][8] == "water" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "oil" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
                elif dataCarbonate[0][8] == "gas" and len(dataCarbonate) > 0:
                    for i in range(0, len(dataCarbonate)):
                        sequences.append(dataCarbonate[i])
                        actualThickness += dataCarbonate[i][1]
                        seqCounter += 1
            elif magicnumber in range(60,75):
                ##############
                # EVAPORITES #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                rocks = ["Halite", "Anhydrite"]
                magicnumber = rd.randint(0, 1)
                dataEvaporites = data.createEvaporites(rocks[magicnumber])
                for i in range(0, len(dataEvaporites)):
                    sequences.append(dataEvaporites[i])
                    actualThickness += dataEvaporites[i][1]
                    seqCounter += 1
            elif magicnumber in range(75,100):
                ##############
                # ORE (IRON) #
                ##############
                rocktype = sequences[seqCounter][0]
                data = subsurface(actualThickness, rocktype, parts)
                data_ore_Fe = data.create_ores()
                for i in range(0, len(data_ore_Fe)):
                    sequences.append(data_ore_Fe[i])
                    actualThickness += data_ore_Fe[i][1]
                    seqCounter += 1
        else:
            ###########
            # IGNEOUS #
            ###########
            rocktype = sequences[seqCounter][0]
            data = subsurface(actualThickness, rocktype, parts)
            dataIgneous = data.createIgneousRocks()
            for i in range(0, len(dataIgneous)):
                sequences.append(dataIgneous[i])
                actualThickness += dataIgneous[i][1]
                seqCounter += 1
    #
    ################
    # AFTER BASALT #
    ################
    elif sequences[seqCounter][0] == "basalt":
        ###########
        # IGNEOUS #
        ###########
        rocktype = sequences[seqCounter][0]
        data = subsurface(actualThickness, rocktype, parts)
        dataIgneous = data.createIgneousRocks(rock="basalt")
        for i in range(0, len(dataIgneous)):
            sequences.append(dataIgneous[i])
            actualThickness += dataIgneous[i][1]
            seqCounter += 1
    #
    #################
    # AFTER GRANITE #
    #################
    elif sequences[seqCounter][0] == "granite":
        ###########
        # IGNEOUS #
        ###########
        rocktype = sequences[seqCounter][0]
        data = subsurface(actualThickness, rocktype, parts)
        dataIgneous = data.createIgneousRocks(rock="granite")
        for i in range(0, len(dataIgneous)):
            sequences.append(dataIgneous[i])
            actualThickness += dataIgneous[i][1]
            seqCounter += 1
#
## EXPORT SEQUENCES
#
print("")
print("FINISHED SEQUENCE STRATIGRAPHY:")
sequenceLog = []
lithologies = []
nSequences = len(sequences)
thicknessData = []
topData = []
bottomData = [0]
depthData = [0]
dataDepth = [[0], [0], [0]]
densitiesData = []
densitiesBottom = []
velocities = []
vP = []
vS = []
vFluid = []
logGR = []
phiN = []
phi = []
mineralogy = []
nu = []
PE = []
PEF = []
AI = []
for i in range(0, nSequences):
    sequenceLog.append(sequences[i][0])
    thicknessData.append(sequences[i][1])
    depthData.append(depthData[i] - thicknessData[i])
    topData.append(-sequences[i][2])
    bottomData.append(-sequences[i][3])
    dataDepth[0].append(depthData[i] - thicknessData[i])
    dataDepth[1].append(sequences[i][1])
    dataDepth[2].append(dataDepth[0][i] - dataDepth[1][i + 1] / 2)
    densitiesData.append(sequences[i][4])
    densitiesBottom.append(sequences[i][4])
    velocities.append(sequences[i][5])
    vP.append(velocities[i][0])
    vS.append(velocities[i][1])
    vFluid.append(velocities[i][3])
    logGR.append(sequences[i][6])
    phiN.append(sequences[i][7])
    nu.append(sequences[i][10])
    PE.append(sequences[i][11])
    AI.append(sequences[i][4]*sequences[i][5][0])
    print(i + 1, "unit:", sequences[i])
    if sequences[i][0] not in lithologies:
        lithologies.append(sequences[i][0])
AI.append(AI[-1])
for i in range(2, nSequences):
    mineralogy.append(sequences[i][9])
print("Final thickness:", round(sum(thicknessData),1))
topData.append(-sum(thicknessData))
sortedUnits = []
for i in range(0, len(lithologies)):
    sortedUnits.append([lithologies[i]])
for i in range(0, len(lithologies)):
    for j in range(0, len(sequences)):
        if sequences[j][0] in sortedUnits[i][0]:
            sortedUnits[i].append([-sequences[j][2], -sequences[j][3]])
        else:
            continue
#
## DATA
print("")
print("Lithologies:", len(lithologies))
print(lithologies)
# print("")
# print(len(dataDepth[0]), dataDepth[0])
# print(len(dataDepth[1]), dataDepth[1])
# print(len(dataDepth[2]), dataDepth[2])
# print("Thickness of the units:", len(thicknessData))
# print(thicknessData)
# print("")
# print("Depth of the units:", len(depthData)-1)
del depthData[-1]
# print(depthData)
# print("")
# print("Top of the units:", len(topData))
# print(topData)
# print("")
# print("Bottom of the units:", len(bottomData))
# print(bottomData)
# print("")
# print("Bottom densities of the units:", len(densitiesBottom)+1)
densitiesBottom.append(sequences[nSequences - 1][4])
# print(densitiesBottom)
# print("")
# print("Densities of the units:", len(densitiesData))
# print(densitiesData)
# print("")
# print("Velocities (p-wave) of the units:", len(vP))
vP.append(sequences[nSequences - 1][5][0])
# print(vP)
# print("")
# print("GR log:", len(logGR))
logGR.append(sequences[nSequences - 1][6])
# print(logGR)
# print("")
# print("Neutron log:", len(phiN))
phiN.append(sequences[nSequences - 1][7])
# print(phiN)
# print("")
# print("Mineralogy:", len(mineralogy))
mineralogy.append(sequences[nSequences - 1][7])
# print(mineralogy)
nu.append(sequences[nSequences-1][10])
PE.append(sequences[nSequences-1][11])

## GEOPHYISCS
# print("")
data = geophysics(sequences)
dataPressure = data.calculatePressure()
# print("Pressure:", len(dataPressure))
# print(dataPressure)
# print("")
data = geophysics(sequences)
dataTemperature = data.calculateTemperature()
# print("Temperature:", len(dataTemperature))
# print(dataTemperature)
# print("")
data = geophysics(sequences)
dataTTI = data.calculateTTI()
# print("Transit Time Interval:", len(dataTTI))
# print(dataTTI)

## SEISMICS
data = seismics(sequences)
dataImpedance = data.calculateImpedance()
data = seismics(dataImpedance)
dataReflection = data.calculateReflection()
RC = dataReflection
RC.append(dataReflection[-1])
data = seismics(velocities)
dtc, dts = data.calculate_travel_times()
# data = seismics(dataReflection)
# dataWavelet = data.calculateRickerWavelet()
data = seismics(sequences)
dataT0 = data.calculatet0()
data = seismics(dataReflection)
dataTrace = data.calculateTrace()
# data = seismics(dataWavelet)
# dataWiggle = data.seismic_wiggle()
# print("")
#print("Impedance:", len(dataImpedance))
#data = [round(i, 2) for i in dataImpedance]
#print(data)
#print(AI)
# print("")
#print("Reflection:", len(dataReflection))
#data = [round(i, 2) for i in dataReflection]
#print(data)
# print("")
# print("Wavelet (Ricker):", len(dataWavelet))
# data = [round(i,2) for i in dataWavelet]
# print(data)
# print("")
# print("Minimum travel time:", len(dataT0))
# data = [round(i,2) for i in dataT0]
# print(dataT0)
# print("")
# print("Trace:", len(dataTrace[0]), len(dataTrace[0][0]))
# data = [round(i,2) for i in dataTrace[0]]
# print(dataTrace[0][0])

## ANALYSIS
maxDepth = -sum(thicknessData)

## STRATIGRAPHY
pattern = ["-", "+", "x", "\\", "*", "o", "O", ".", "/"]
nUnits = []
usedLithologies = []
for i in range(0, len(sortedUnits)):
    if sortedUnits[i][0] == "water":
        sortedUnits[i].append("//")
        sortedUnits[i].append("royalblue")
        usedLithologies.append(["water", "", "royalblue"])
    elif sortedUnits[i][0] == "dry sand":
        sortedUnits[i].append("..")
        sortedUnits[i].append("gold")
        usedLithologies.append(["dry sand", "", "gold"])
    elif sortedUnits[i][0] == "wet sand":
        sortedUnits[i].append("..")
        sortedUnits[i].append("palegoldenrod")
        usedLithologies.append(["wet sand", "", "palegoldenrod"])
    elif sortedUnits[i][0] == "scree":
        sortedUnits[i].append("xx")
        sortedUnits[i].append("darkkhaki")
        usedLithologies.append(["scree", "", "darkkhaki"])
    elif sortedUnits[i][0] == "soil":
        sortedUnits[i].append("xx")
        sortedUnits[i].append("saddlebrown")
        usedLithologies.append(["soil", "", "saddlebrown"])
    elif sortedUnits[i][0] == "sandstone":
        sortedUnits[i].append("oo")
        sortedUnits[i].append("burlywood")
        usedLithologies.append(["sandstone", "", "burlywood"])
    elif sortedUnits[i][0] == "limestone":
        sortedUnits[i].append("++")
        sortedUnits[i].append("lightblue")
        usedLithologies.append(["limestone", "", "lightblue"])
    elif sortedUnits[i][0] == "dolomite":
        sortedUnits[i].append("++")
        sortedUnits[i].append("violet")
        usedLithologies.append(["dolomite", "", "violet"])
    elif sortedUnits[i][0] == "shale":
        sortedUnits[i].append("---")
        sortedUnits[i].append("palegreen")
        usedLithologies.append(["shale", "", "palegreen"])
    elif sortedUnits[i][0] == "halite":
        sortedUnits[i].append("**")
        sortedUnits[i].append("thistle")
        usedLithologies.append(["halite", "", "thistle"])
    elif sortedUnits[i][0] == "anhydrite":
        sortedUnits[i].append("**")
        sortedUnits[i].append("pink")
        usedLithologies.append(["anhydrite", "", "pink"])
    elif sortedUnits[i][0] == "granite":
        sortedUnits[i].append("**")
        sortedUnits[i].append("coral")
        usedLithologies.append(["granite", "", "coral"])
    elif sortedUnits[i][0] == "basalt":
        sortedUnits[i].append("++")
        sortedUnits[i].append("silver")
        usedLithologies.append(["basalt", "", "silver"])
    elif sortedUnits[i][0] == "ore":
        sortedUnits[i].append("++")
        sortedUnits[i].append("brown")
        usedLithologies.append(["ore", "", "brown"])
    nUnits.append(len(sortedUnits[i]) - 3)
usedPattern = []
usedColors = []
legendLith = []
for i in range(0, len(lithologies)):
    usedPattern.append(usedLithologies[i][1])
    usedColors.append(usedLithologies[i][2])
    legendLith.append(mpatches.Patch(facecolor=usedLithologies[i][2], hatch=usedLithologies[i][1], label=usedLithologies[i][0]))
# print("")
# print("Stratigraphy:", sortedUnits)
# print("nUnits:", nUnits)
# print("Used Lithologies:", usedLithologies)
# fig, ax = plt.subplots(1)
# for i in range(0, len(nUnits)):
#     for j in range(1, nUnits[i] + 1):
#         # ax.hist(np.linspace(sortedUnits[i][j][0], sortedUnits[i][j][1]), bins=10, color=sortedUnits[i][-1], hatch=sortedUnits[i][-2], orientation="horizontal")
#         ax.hist(np.linspace(sortedUnits[i][j][0], sortedUnits[i][j][1]), bins=10, color=sortedUnits[i][-1],
#                 orientation="horizontal")
# plt.title("Sequence Stratigraphy")
# plt.ylabel("Depth [m]")
# ax.set_xticklabels([])
# plt.xlim(0, 1)
# plt.ylim(min(depthData), 1.1 * max(depthData))
# plt.yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
# plt.grid(True)
# plt.grid(which="major", linewidth=1.0)
# plt.grid(which="minor", linewidth=0.5)
# ax.set_axisbelow(True)
# ax.legend(handles=legendLith, loc="lower center", bbox_to_anchor=(0.5, -0.2), shadow=True, ncol=3, prop={'size': 8})
# plt.show()

## SEQUENCE LOG
# seqX = []
# seqY = []
# units = []
# for i in range(0, len(thicknessData)):
#     seqX.append(thicknessData[i])
#     seqY.append(depthData[i])
#     if sequences[i][0] == "water":
#         units.append([thicknessData[i], depthData[i], " "])
#     elif sequences[i][0] == "scree":
#         units.append([thicknessData[i], depthData[i], "O"])
#     elif sequences[i][0] == "soil":
#         units.append([thicknessData[i], depthData[i], "O"])
#     elif sequences[i][0] == "dry sand":
#         units.append([thicknessData[i], depthData[i], "+"])
#     elif sequences[i][0] == "wet sand":
#         units.append([thicknessData[i], depthData[i], "x"])
#     elif sequences[i][0] == "sandstone":
#         units.append([thicknessData[i], depthData[i], "."])
#     elif sequences[i][0] == "limestone":
#         units.append([thicknessData[i], depthData[i], "-"])
#     elif sequences[i][0] == "dolomite":
#         units.append([thicknessData[i], depthData[i], "-"])
#     elif sequences[i][0] == "shale":
#         units.append([thicknessData[i], depthData[i], "/"])
#     elif sequences[i][0] == "halite":
#         units.append([thicknessData[i], depthData[i], "/"])
#     elif sequences[i][0] == "anhydrite":
#         units.append([thicknessData[i], depthData[i], "/"])
#     elif sequences[i][0] == "granite":
#         units.append([thicknessData[i], depthData[i], "*"])
#     elif sequences[i][0] == "basalt":
#         units.append([thicknessData[i], depthData[i], "+"])
#     elif sequences[i][0] == "ore":
#         units.append([thicknessData[i], depthData[i], "*"])
# seqX.append(0.0)
# seqY.append(maxDepth)

## Pressure-Temperature LOG
# x = []
# for i in range(0, len(dataPressure)):
#     x.append(dataPressure[i] / 1000000)
# y = dataDepth[2]
# fig, ax1 = plt.subplots()
# ax1.plot(x, y, color="green", linewidth=2)
# ax1.set_xlabel("Pressure p [MPa]")
# ax1.set_ylabel("Depth [m]")
# ax1.set_xlim(0.9 * min(x), 1.1 * max(x))
# ax1.set_ylim(min(depthData), 1.1 * max(depthData))
# ax1.set_xticks(np.arange(round(0.9 * min(x), 1), round(1.1 * max(x), 1), 5))
# ax1.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
# ax1.grid(True)
# ax1.grid(which="major", linewidth=1.0)
# ax1.grid(which="minor", linewidth=0.5)
# plt.rc('axes', axisbelow=True)
# ax2 = ax1.twiny()
# ax2.plot(dataTemperature, y, color="red", linewidth=2)
# ax2.set_xlabel("Temperature T [°C]")
# ax2.set_xticks(np.arange(round(0.9 * min(dataTemperature), 1), round(1.1 * max(dataTemperature), 1), 5))
# plots = [mpatches.Patch(color="green", label="Pressure"), mpatches.Patch(color="red", label="Temperature")]
# ax2.legend(handles=plots, loc="lower left", ncol=3, prop={'size': 8})
# plt.show()

## DENSITY VS. VELOCITY
# x = np.asarray(vP)
# y = np.asarray(densitiesBottom)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     elif sequenceLog[i] == "ore":
#         coordOre[0].append(x[i])
#         coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# if len(coordOre[0]) > 0:
#     plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(scatterpoints=1, prop={'size': 8}, loc="upper left", ncol=3)
# plt.title("Relationship between density and p-wave velocity")
# plt.xlabel("$v_P$ velocity [m/s]")
# plt.ylabel("Density $\\varrho$ [g/cm^3]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# plt.ylim(1.8, 3.4)
# plt.xticks(np.arange(round(0.9 * min(x) / 500, 0) * 500, round(1.1 * max(x) / 500, 0) * 500, 1000))
# plt.yticks(np.arange(1.8, 3.5, 0.1))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

## DENSITY VS. PHOTOELECTRIC
# x = np.asarray(densitiesBottom)
# y = np.asarray(PE)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     elif sequenceLog[i] == "ore":
#         coordOre[0].append(x[i])
#         coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# if len(coordOre[0]) > 0:
#     plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(scatterpoints=1, prop={'size': 8}, loc="upper left", ncol=2)
# plt.title("Relationship between density and photoelectricity")
# plt.xlabel("Density $\\varrho$ [g/cm^3]")
# plt.ylabel("PE (barns/electron))")
# plt.yscale('log')
# plt.xlim(1.8, 3.4)
# #plt.ylim(0, 50)
# plt.xticks(np.arange(1.8, 3.5, 0.1))
# #plt.yticks(np.arange(0, 50, 5))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

# ## POISSON'S RATIO VS. IMPEDANCE
# X = []
# Y = []
# # xData = vP
# # del xData[-1]
# for i in range(1, len(dataImpedance)):
#     Y.append(dataImpedance[i] / 1000000)
# # for j in range(1, len(xData)):
# # Y.append(xData[j]/vS[j])
# for j in range(1, len(sequences)):
#     X.append(sequences[j][10])
# x = np.asarray(X)
# y = np.asarray(Y)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     #elif sequenceLog[i] == "ore":
#     #    coordOre[0].append(x[i])
#     #    coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# #if len(coordOre[0]) > 0:
# #    plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(loc="upper left", scatterpoints=1, prop={'size': 8}, ncol=2)
# plt.xlabel("Poisson's ratio $\\nu$ [1]")
# plt.ylabel("Impedance $Z$ [$MNs/m^3$]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# plt.ylim(min(y), 1.1 * max(y))
# plt.xticks(np.arange(0.8 * min(x), 1.2 * max(x), 0.2))
# plt.yticks(np.arange(round(0.9 * min(y) / 0.2, 0) * 0.2, round(1.1 * max(y) / 0.2, 0) * 0.2, 1))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

# ## POISSON'S RATIO VS. vP/vS
# X = []
# Y = []
# xData = vP
# del xData[-1]
# for i in range(1, len(sequences)):
#     Y.append(sequences[i][5][0] / sequences[i][5][1])
# # for j in range(1, len(xData)):
# # Y.append(xData[j]/vS[j])
# for j in range(1, len(sequences)):
#     X.append(sequences[j][10])
# x = np.asarray(X)
# y = np.asarray(Y)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     elif sequenceLog[i] == "ore":
#         coordOre[0].append(x[i])
#         coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# if len(coordOre[0]) > 0:
#     plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(loc="upper left", scatterpoints=1, prop={'size': 8}, ncol=3)
# plt.xlabel("Poisson's ratio $\\nu$ [1]")
# plt.ylabel("$v_p/v_s$ [1]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# # plt.ylim(min(y), 1.1*max(y))
# plt.xticks(np.arange(0.8 * min(x), 1.2 * max(x), 0.2))
# # plt.yticks(np.arange(round(0.9*min(y)/0.2,0)*0.2, round(1.1*max(y)/0.2,0)*0.2, 0.1))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()
#
# ## IMPEDANCE VS. vP/vS
# X = []
# Y = []
# for i in range(1, len(dataImpedance)):
#     X.append(dataImpedance[i] / 1000000)
# for i in range(1, len(sequences)):
#     Y.append(sequences[i][5][0] / sequences[i][5][1])
# x = np.asarray(X)
# y = np.asarray(Y)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     elif sequenceLog[i] == "ore":
#         coordOre[0].append(x[i])
#         coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# if len(coordOre[0]) > 0:
#     plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(scatterpoints=1, prop={'size': 8}, ncol=3)
# plt.xlabel("Impedance $Z$ [$MNs/m^3$]")
# plt.ylabel("$v_p/v_s$ [1]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# plt.ylim(min(y), 1.1 * max(y))
# plt.xticks(np.arange(0.8 * min(x), 1.2 * max(x), 2))
# plt.yticks(np.arange(round(0.9 * min(y) / 0.2, 0) * 0.2, round(1.1 * max(y) / 0.2, 0) * 0.2, 0.2))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()
#
# ## MINERALOGY
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordGranite = [[], []]
# for i in range(2, len(sequenceLog)):
#     if sequenceLog[i] == "sandstone":
#         n = 6
#         for j in range(0, n):
#             if mineralogy[i - 2][j][0] == "Qz":
#                 coordSandstone[0].append(mineralogy[i - 2][j][1])
#             elif mineralogy[i - 2][j][0] == "Cal":
#                 coordSandstone[1].append(mineralogy[i - 2][j][1])
#     elif sequenceLog[i] == "limestone":
#         n = 4
#         for j in range(0, n):
#             if mineralogy[i - 2][j][0] == "Qz":
#                 coordLimestone[0].append(mineralogy[i - 2][j][1])
#             elif mineralogy[i - 2][j][0] == "Cal":
#                 coordLimestone[1].append(mineralogy[i - 2][j][1])
#     elif sequenceLog[i] == "dolomite":
#         n = 4
#         for j in range(0, n):
#             if mineralogy[i - 2][j][0] == "Qz":
#                 coordDolomite[0].append(mineralogy[i - 2][j][1])
#             elif mineralogy[i - 2][j][0] == "Cal":
#                 coordDolomite[1].append(mineralogy[i - 2][j][1])
#     elif sequenceLog[i] == "shale":
#         n = 7
#         for j in range(0, n):
#             if mineralogy[i - 2][j][0] == "Qz":
#                 coordShale[0].append(mineralogy[i - 2][j][1])
#             elif mineralogy[i - 2][j][0] == "Cal":
#                 coordShale[1].append(mineralogy[i - 2][j][1])
# # elif sequenceLog[i] == "granite":
# #	n = 6
# #	for j in range(0, n):
# #		if mineralogy[i-2][j][0] == "Qz":
# #			coordGranite[0].append(mineralogy[i-2][j][1])
# #		elif mineralogy[i-2][j][0] == "Cal":
# #			coordGranite[1].append(mineralogy[i-2][j][1])
# plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# # plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# plt.legend(scatterpoints=1, prop={'size': 8}, ncol=3)
# plt.xlabel("$x_{Qz}$")
# plt.ylabel("$x_{Cal}$")
# plt.xlim(0.0, 1.0)
# plt.ylim(0.0, 1.0)
# plt.xticks(np.arange(0.0, 1.1, 0.1))
# plt.yticks(np.arange(0.0, 1.1, 0.1))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

## POROSITY VS. VELOCITY
phiNedited = phiN
del phiNedited[-1]
# x = np.asarray(vP)
# y = np.asarray(phiNedited)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     elif sequenceLog[i] == "ore":
#         coordOre[0].append(x[i])
#         coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# if len(coordOre[0]) > 0:
#     plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(scatterpoints=1, prop={'size': 8}, loc="upper left", ncol=3)
# plt.title("Relationship between porosity and p wave velocity")
# plt.xlabel("$v_P$ velocity [m/s]")
# plt.ylabel("Neutron porosity [%]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# plt.ylim(0, 50)
# plt.xticks(np.arange(round(0.9 * min(x) / 500, 0) * 500, round(1.1 * max(x) / 500, 0) * 500, 1000))
# plt.yticks(np.arange(0, 60, 10))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

## vP vs. vS
# x = np.asarray(vP)
# y = np.asarray(vS)
# coordDrySand = [[], []]
# coordWetSand = [[], []]
# coordScree = [[], []]
# coordSoil = [[], []]
# coordSandstone = [[], []]
# coordLimestone = [[], []]
# coordDolomite = [[], []]
# coordShale = [[], []]
# coordHalite = [[], []]
# coordAnhydrite = [[], []]
# coordGranite = [[], []]
# coordBasalt = [[], []]
# coordOre = [[], []]
# for i in range(0, len(sequenceLog)):
#     if sequenceLog[i] == "water":
#         continue
#     elif sequenceLog[i] == "dry sand":
#         coordDrySand[0].append(x[i])
#         coordDrySand[1].append(y[i])
#     elif sequenceLog[i] == "wet sand":
#         coordWetSand[0].append(x[i])
#         coordWetSand[1].append(y[i])
#     elif sequenceLog[i] == "scree":
#         coordScree[0].append(x[i])
#         coordScree[1].append(y[i])
#     elif sequenceLog[i] == "soil":
#         coordSoil[0].append(x[i])
#         coordSoil[1].append(y[i])
#     elif sequenceLog[i] == "sandstone":
#         coordSandstone[0].append(x[i])
#         coordSandstone[1].append(y[i])
#     elif sequenceLog[i] == "limestone":
#         coordLimestone[0].append(x[i])
#         coordLimestone[1].append(y[i])
#     elif sequenceLog[i] == "dolomite":
#         coordDolomite[0].append(x[i])
#         coordDolomite[1].append(y[i])
#     elif sequenceLog[i] == "shale":
#         coordShale[0].append(x[i])
#         coordShale[1].append(y[i])
#     elif sequenceLog[i] == "halite":
#         coordHalite[0].append(x[i])
#         coordHalite[1].append(y[i])
#     elif sequenceLog[i] == "anhydrite":
#         coordAnhydrite[0].append(x[i])
#         coordAnhydrite[1].append(y[i])
#     elif sequenceLog[i] == "granite":
#         coordGranite[0].append(x[i])
#         coordGranite[1].append(y[i])
#     elif sequenceLog[i] == "basalt":
#         coordBasalt[0].append(x[i])
#         coordBasalt[1].append(y[i])
#     elif sequenceLog[i] == "ore":
#         coordOre[0].append(x[i])
#         coordOre[1].append(y[i])
# if len(coordDrySand[0]) > 0:
#     plt.scatter(coordDrySand[0], coordDrySand[1], c="gold", label="dry sand", alpha=0.8)
# else:
#     pass
# if len(coordSoil[0]) > 0:
#     plt.scatter(coordSoil[0], coordSoil[1], c="saddlebrown", label="soil", alpha=0.8)
# else:
#     pass
# plt.scatter(coordWetSand[0], coordWetSand[1], c="palegoldenrod", label="wet sand", alpha=0.8)
# if len(coordSandstone[0]) > 0:
#     plt.scatter(coordSandstone[0], coordSandstone[1], c="peru", label="sandstone", alpha=0.8)
# if len(coordLimestone[0]) > 0:
#     plt.scatter(coordLimestone[0], coordLimestone[1], c="steelblue", label="limestone", alpha=0.8)
# if len(coordDolomite[0]) > 0:
#     plt.scatter(coordDolomite[0], coordDolomite[1], c="darkviolet", label="dolomite", alpha=0.8)
# if len(coordShale[0]) > 0:
#     plt.scatter(coordShale[0], coordShale[1], c="limegreen", label="shale", alpha=0.8)
# if len(coordHalite[0]) > 0:
#     plt.scatter(coordHalite[0], coordHalite[1], c="thistle", label="halite", alpha=0.8)
# if len(coordAnhydrite[0]) > 0:
#     plt.scatter(coordAnhydrite[0], coordAnhydrite[1], c="pink", label="anhydrite", alpha=0.8)
# if len(coordGranite[0]) > 0:
#     plt.scatter(coordGranite[0], coordGranite[1], c="coral", label="granite", alpha=0.8)
# if len(coordBasalt[0]) > 0:
#     plt.scatter(coordBasalt[0], coordBasalt[1], c="silver", label="basalt", alpha=0.8)
# if len(coordOre[0]) > 0:
#     plt.scatter(coordOre[0], coordOre[1], c="brown", label="ore", alpha=0.8)
# plt.legend(scatterpoints=1, prop={'size': 8}, loc="lower right", ncol=3)
# plt.title("Relationship between the P- and S-wave velocities")
# plt.xlabel("$v_P$ [m/s]")
# plt.ylabel("$v_S$ [m/s]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# plt.ylim(0.9 * min(y), 1.1 * max(y))
# plt.xticks(np.arange(round(0.9 * min(x) / 500, 0) * 500, round(1.1 * max(x) / 500, 0) * 500 + 500, 1000))
# plt.yticks(np.arange(round(0.9 * min(y) / 500, 0) * 500, round(1.1 * max(y) / 500, 0) * 500, 500))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

# ## vP vs. vS vs. vF
# x = np.asarray(vP)
# y = np.asarray(vS)
# z = np.asarray(vFluid)
# plt.scatter(x, y, c=z, cmap="viridis", norm=None, s=30, alpha=0.8)
# cb = plt.colorbar()
# cb.set_label("$v_F$ [m/s]")
# plt.title("Influence of the fluid velocity on the seismic velocities")
# plt.xlabel("$v_P$ [m/s]")
# plt.ylabel("$v_S$ [m/s]")
# plt.xlim(0.9 * min(x), 1.1 * max(x))
# plt.ylim(0.9 * min(y), 1.1 * max(y))
# plt.xticks(np.arange(round(0.9 * min(x) / 500, 0) * 500, round(1.1 * max(x) / 500, 0) * 500 + 500, 1000))
# plt.yticks(np.arange(round(0.9 * min(y) / 500, 0) * 500, round(1.1 * max(y) / 500, 0) * 500, 500))
# plt.grid(True)
# plt.rc('axes', axisbelow=True)
# plt.show()

################
# Scatter Plot #
#  Lithology   #
################
dataPlot = plotting.scatter()
dataPlotHist = plotting.histogram()
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Qz", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Cal", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Kfs", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Pl", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Bt", lith="sandstone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Glt", lith="sandstone")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="limestone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="limestone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="limestone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Cal", lith="limestone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Dol", lith="limestone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Qz", lith="limestone")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Or", lith="limestone")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="dolomite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="dolomite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="dolomite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Dol", lith="dolomite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Cal", lith="dolomite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Qz", lith="dolomite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Or", lith="dolomite")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Qz", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Kln", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Chl", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Ilt", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Cal", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Dol", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Kfs", lith="shale")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Pl", lith="shale")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Hl", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Anh", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Gp", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Syl", lith="halite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Ilt", lith="halite")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Anh", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Hl", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Cal", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Gn", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Ccp", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Mol", lith="anhydrite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Py", lith="anhydrite")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Qz", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Or", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Ab", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="An", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Bt", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Ms", lith="granite")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Py", lith="granite")
#
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="GR", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="NPHI", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="POISSON", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="An", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Wo", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="En", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Fs", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Tr", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Fo", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Fa", lith="basalt")
# dataPlot.scatterPlotLith(sequences, x="VP", y="RHOB", z="Tep", lith="basalt")

# ROCK DETERMINATION
#found = False
#while found == False:
#    i = rd.randint(0, len(sequences))
#    if sequences[i][0] == "sandstone":
#        measurements = np.array([sequences[i][4], sequences[i][6], sequences[i][11]])
#        feature_list = ["RHOB", "GR", "PEF"]
#        X = RD()
#        x = X.create_sandstone_by_properties(measurements, feature_list)
#        print(x)
#        found = True
#    else:
#        found = False

## MULTIPLOT
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey=True, figsize=(9, 12))
fig.subplots_adjust(wspace=0.25)
# 1
for i in range(0, len(nUnits)):
    for j in range(1, nUnits[i] + 1):
        ax1.hist(np.linspace(sortedUnits[i][j][0], sortedUnits[i][j][1]), bins=10, color=sortedUnits[i][-1], orientation="horizontal")
ax1.set_ylabel("Depth [m]")
ax1.set_xticklabels([])
ax1.set_ylim(min(depthData), 1.1 * max(depthData))
ax1.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax1.grid(True)
ax1.grid(which="major", linewidth=1.0)
ax1.grid(which="minor", linewidth=0.5)
ax1.set_axisbelow(True)
ax1.legend(handles=legendLith, loc="lower center", bbox_to_anchor=(0.25, -0.125), shadow=True, ncol=2, prop={'size': 8})
# 2
ax2.plot(logGR, dataDepth[2], color="#00549F", linewidth=2)
ax2.set_xlabel("Gamma ray (API)")
ax2.set_xlim(0, max(logGR))
ax2.set_xticks(np.arange(0, max(logGR)+50, 100))
ax2.set_ylim(min(depthData), 1.1 * max(depthData))
ax2.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax2.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax2.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax2.minorticks_on()
plt.rc('axes', axisbelow=True)
# 3
vPedit = vP
del vPedit[-1]
vPedit.insert(0, vP[0])
vP_km = [i*0.001 for i in vP]
ax3.plot(vP_km, dataDepth[0], color="#00549F", linewidth=2)
ax3.set_xlabel("vP (km/s)")
ax3.set_xlim(0, 8.5)
ax3.set_xticks(np.arange(0, 8.5, 2.0))
ax3.set_ylim(min(depthData), 1.1 * max(depthData))
ax3.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax3.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax3.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax3.minorticks_on()
#
vSedit = vS
vSedit.insert(0, vS[0])
vS_km = [i*0.001 for i in vS]
ax3_2 = ax3.twiny()
ax3_2.plot(vS_km, dataDepth[0], color="#CC071E", linewidth=2)
ax3_2.set_xlabel("vS (km/s)")
ax3_2.set_xlim(0, 8.5)
ax3_2.set_xticks(np.arange(0, 8.5, 2.0))
ax3_2.minorticks_on()
plt.rc('axes', axisbelow=True)
plots = [mpatches.Patch(color="#00549F", label="VP"), mpatches.Patch(color="#CC071E", label="VS")]
ax3_2.legend(handles=plots, loc="lower right", bbox_to_anchor=(1, -0.125), prop={'size': 8})
# 4
plotRHOB = ax4.plot(densitiesBottom, dataDepth[0], color="#57AB27", linewidth=2, label="Densitiy log")
ax4.set_xlabel("Density [g/cm^3]")
ax4.set_xlim(1.70, 3.30)
ax4.set_xticks(np.around(np.linspace(1.70, 3.30, 4, endpoint=True), decimals=2))
ax4.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax4.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax4.minorticks_on()
#
ax4_2 = ax4.twiny()
phiN.append(sequences[nSequences - 1][7])
plotPHIN = ax4_2.plot(phiN, dataDepth[0], color="#00549F", linewidth=2, linestyle="-", label="Neutron log")
ax4_2.set_xlabel("Neutron porosity [%]")
ax4_2.set_xlim(57, -27)
ax4_2.set_xticks(np.around(np.linspace(57, -27, 6, endpoint=True), decimals=0))
ax4_2.set_ylim(min(depthData), 1.1 * max(depthData))
ax4_2.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax4_2.minorticks_on()
plt.rc('axes', axisbelow=True)
plots = [mpatches.Patch(color="#57AB27", label="RHOB"), mpatches.Patch(color="#00549F", label="PHIN")]
ax4_2.legend(handles=plots, loc="lower right", bbox_to_anchor=(1, -0.125), prop={'size': 8})
# 5
ax5.plot(PE, dataDepth[2], color="#006165", linewidth=2)
ax5.set_xlabel("PE (barns/electron)")
ax5.set_xscale('log')
ax5.set_xlim(0, 40)
ax5.set_ylim(min(depthData), 1.1 * max(depthData))
ax5.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax5.get_xaxis().set_minor_formatter(mpl.ticker.ScalarFormatter())
#ax5.xaxis.set_major_formatter(NullFormatter())
ax5.xaxis.set_minor_formatter(NullFormatter())
ax5.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax5.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax5.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax5.minorticks_on()
plt.rc('axes', axisbelow=True)
#
plt.savefig("./outputs/GebPy_Log_01.png", dpi=300)
plt.show()

## MULTIPLOT (SEISMIC)
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey=True, figsize=(9, 12))
fig.subplots_adjust(wspace=0.25)
# 1
for i in range(0, len(nUnits)):
    for j in range(1, nUnits[i] + 1):
        ax1.hist(np.linspace(sortedUnits[i][j][0], sortedUnits[i][j][1]), bins=10, color=sortedUnits[i][-1], orientation="horizontal")
ax1.set_ylabel("Depth [m]")
ax1.set_xticklabels([])
ax1.set_ylim(min(depthData), 1.1 * max(depthData))
ax1.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax1.grid(True)
ax1.grid(which="major", linewidth=1.0)
ax1.grid(which="minor", linewidth=0.5)
ax1.set_axisbelow(True)
ax1.legend(handles=legendLith, loc="lower center", bbox_to_anchor=(0.25, -0.125), shadow=True, ncol=2, prop={'size': 8})
# 2
plotRHOB = ax2.plot(densitiesBottom, dataDepth[0], color="#57AB27", linewidth=2, label="Densitiy log")
ax2.set_xlabel("Density [g/cm^3]")
ax2.set_xlim(1.70, 3.30)
ax2.set_xticks(np.around(np.linspace(1.70, 3.30, 4, endpoint=True), decimals=2))
ax2.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax2.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax2.minorticks_on()
#
ax2_2 = ax2.twiny()
plotPHIN = ax2_2.plot(phiN, dataDepth[0], color="#00549F", linewidth=2, linestyle="-", label="Neutron log")
ax2_2.set_xlabel("Neutron porosity [%]")
ax2_2.set_xlim(57, -27)
ax2_2.set_xticks(np.around(np.linspace(57, -27, 6, endpoint=True), decimals=0))
ax2_2.set_ylim(min(depthData), 1.1 * max(depthData))
ax2_2.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax2_2.minorticks_on()
plt.rc('axes', axisbelow=True)
plots = [mpatches.Patch(color="#57AB27", label="RHOB"), mpatches.Patch(color="#00549F", label="PHIN")]
ax2_2.legend(handles=plots, loc="lower right", bbox_to_anchor=(1, -0.125), prop={'size': 8})
# 3
vPedit = vP
vP_km = [i*0.001 for i in vP]
ax3.plot(vP_km, dataDepth[0], color="#00549F", linewidth=2)
ax3.set_xlabel("vP (km/s)")
ax3.set_xlim(0, 8)
ax3.set_xticks(np.arange(0, 8.5, 2.0))
ax3.set_ylim(min(depthData), 1.1 * max(depthData))
ax3.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax3.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax3.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax3.minorticks_on()
#
vSedit = vS
vS_km = [i*0.001 for i in vS]
ax3_2 = ax3.twiny()
ax3_2.plot(vS_km, dataDepth[0], color="#CC071E", linewidth=2)
ax3_2.set_xlabel("vS (km/s)")
ax3_2.set_xlim(0, 8)
ax3_2.set_xticks(np.arange(0, 8.5, 2.0))
ax3_2.minorticks_on()
plt.rc('axes', axisbelow=True)
plots = [mpatches.Patch(color="#00549F", label="VP"), mpatches.Patch(color="#CC071E", label="VS")]
ax3_2.legend(handles=plots, loc="lower right", bbox_to_anchor=(1, -0.125), prop={'size': 8})
# 4
AI = np.asarray(AI)/1000
ax4.plot(AI, dataDepth[2], color="#006165", linewidth=2)
ax4.set_xlabel("AI (kNs/m$3$)")
#ax4.set_xscale('log')
ax4.set_xlim(0, 25)
ax4.set_xticks(np.around(np.linspace(0, 30, 6, endpoint=True), decimals=0))
ax4.set_ylim(min(depthData), 1.1 * max(depthData))
ax4.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax4.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax4.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax4.minorticks_on()
plt.rc('axes', axisbelow=True)
# 5
ax5.plot(RC, dataDepth[2], color="#006165", linewidth=2)
ax5.set_xlabel("R (1)")
#ax5.set_xscale('symlog')
#ax5.set_xlim(min(RC), max(RC))
ax5.set_xlim(-0.5, 0.5)
ax5.set_ylim(min(depthData), 1.1 * max(depthData))
ax5.set_yticks(np.arange(0, -1.05 * sum(thicknessData), -200))
ax5.grid(b=True, which='major', color="black", linestyle="solid", linewidth=0.4, alpha=1)
ax5.grid(b=True, which='minor', color="black", linestyle="solid", linewidth=0.2, alpha=1)
ax5.minorticks_on()
plt.rc('axes', axisbelow=True)
#
plt.savefig("./outputs/GebPy_Log_02.png", dpi=300)
plt.show()

## GEOCHEMISTRY
data = sample(sequences, 0.05)
dataCore = data.createCore()
data = core.geochemistry(dataCore, sequences)
results = data.calculateConcentration()
data = dataanalysis(lithologies, sequences)
sortedGeodata = data.sortGeodata()
nLiths = len(legendLith)
print("")
print("DATA ANALYSIS")
for i in range(0, nLiths):
    print("{0:9s}".format(sortedGeodata[i][0]) + ":", round(np.mean(sortedGeodata[i][1]), 3), round(np.mean(sortedGeodata[i][2]), 1), round(np.mean(sortedGeodata[i][3]), 1), round(np.mean(sortedGeodata[i][4]), 1), round(np.mean(sortedGeodata[i][5]), 1))