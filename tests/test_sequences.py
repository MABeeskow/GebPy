#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_sequences.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		14.01.2021

# -----------------------------------------------

## MODULES
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from modules import sequences

## TESTING
# Test soil generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_soil = data.create_soil()
for i in range(len(data_soil)):
    print(data_soil[i])

print("")
# Test sand generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sand = data.create_sand()
for i in range(len(data_sand)):
    print(data_sand[i])

print("")
# Test sandstone generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone()
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(fluid="gas")
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
data = sequences.SedimentaryBasin()
data_sandstone = data.create_sandstone(fluid="oil")
for i in range(len(data_sandstone)):
    print(data_sandstone[i])

print("")
# Test shale generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_shale = data.create_shale()
for i in range(len(data_shale)):
    print(data_shale[i])

print("")
# Test rock salt generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_rocksalt = data.create_rocksalt()
for i in range(len(data_rocksalt)):
    print(data_rocksalt[i])

print("")
# Test sedimentary basin generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_sedbasin = data.create_sedimentary_basin(maximum_thickness=500)
rock = []
top = []
bottom = []
rho = []
vP = []
vS = []
vPvS = []
phi = []
gr = []
pe = []
rock_sorted = []
rho_sorted = []
vP_sorted = []
vS_sorted = []
vPvS_sorted = []
phi_sorted = []
gr_sorted = []
pe_sorted = []
for i in range(len(data_sedbasin)):
    if data_sedbasin[i][0][0] not in rock_sorted:
        rock_sorted.append(data_sedbasin[i][0][0])
        rho_sorted.append([data_sedbasin[i][0][0], []])
        vP_sorted.append([data_sedbasin[i][0][0], []])
        vS_sorted.append([data_sedbasin[i][0][0], []])
        vPvS_sorted.append([data_sedbasin[i][0][0], []])
        phi_sorted.append([data_sedbasin[i][0][0], []])
        gr_sorted.append([data_sedbasin[i][0][0], []])
        pe_sorted.append([data_sedbasin[i][0][0], []])
for i in range(len(data_sedbasin)):
    for j in range(len(data_sedbasin[i])):
        for k in range(len(rock_sorted)):
            if data_sedbasin[i][0][0] == rho_sorted[k][0]:
                rho_sorted[k][1].append(data_sedbasin[i][j][4][1][0])
                vP_sorted[k][1].append(data_sedbasin[i][j][4][3][0])
                vS_sorted[k][1].append(data_sedbasin[i][j][4][3][1])
                vPvS_sorted[k][1].append(data_sedbasin[i][j][4][3][0]/data_sedbasin[i][j][4][3][1])
                phi_sorted[k][1].append(data_sedbasin[i][j][4][4][0])
                gr_sorted[k][1].append(data_sedbasin[i][j][4][6][0])
                pe_sorted[k][1].append(data_sedbasin[i][j][4][6][1])
for i in range(len(data_sedbasin)):
    for j in range(len(data_sedbasin[i])):
        rock.append(data_sedbasin[i][j][0])
        top.append(data_sedbasin[i][j][2])
        bottom.append(data_sedbasin[i][j][3])
        rho.append(data_sedbasin[i][j][4][1][0])
        vP.append(data_sedbasin[i][j][4][3][0])
        vS.append(data_sedbasin[i][j][4][3][1])
        vPvS.append(data_sedbasin[i][j][4][3][0]/data_sedbasin[i][j][4][3][1])
        phi.append(data_sedbasin[i][j][4][4][0])
        gr.append(data_sedbasin[i][j][4][6][0])
        pe.append(data_sedbasin[i][j][4][6][1])
    print(data_sedbasin[i])

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vP_sorted[i][1], label=rock_sorted[i], alpha=0.8)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$v_P$")
plt.legend()
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vS_sorted[i][1], label=rock_sorted[i], alpha=0.8)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$v_S$")
plt.legend()
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vPvS_sorted[i][1], label=rock_sorted[i], alpha=0.8)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$v_P/v_S$")
plt.legend()
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], phi_sorted[i][1], label=rock_sorted[i], alpha=0.8)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$\\varphi$")
plt.legend()
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], gr_sorted[i][1], label=rock_sorted[i], alpha=0.8)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("GR")
plt.legend()
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], pe_sorted[i][1], label=rock_sorted[i], alpha=0.8)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("PE")
plt.legend()
plt.show()