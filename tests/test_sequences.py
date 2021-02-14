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
import  matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
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
# Test limestone generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_limestone = data.create_limestone()
for i in range(len(data_limestone)):
    print(data_limestone[i])

print("")
data = sequences.SedimentaryBasin()
data_limestone = data.create_limestone(fluid="gas")
for i in range(len(data_limestone)):
    print(data_limestone[i])

print("")
data = sequences.SedimentaryBasin()
data_limestone = data.create_limestone(fluid="oil")
for i in range(len(data_limestone)):
    print(data_limestone[i])

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
# Test granite generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_granite = data.create_granite()
for i in range(len(data_granite)):
    print(data_granite[i])

print("")
# Test basalt generation within SedimentaryBasin class
data = sequences.SedimentaryBasin()
data_basalt = data.create_basalt()
for i in range(len(data_basalt)):
    print(data_basalt[i])

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

colors = [["soil", "peru"], ["sand", "sandybrown"], ["sandstone", "tan"], ["limestone", "lightblue"], ["shale", "olivedrab"], ["rock salt", "hotpink"], ["granite", "darkorange"], ["basalt", "grey"]]
units_sorted = []
for i in range(len(rock_sorted)):
    units_sorted.append([rock_sorted[i]])
    for j in range(len(data_sedbasin)):
        for k in range(len(data_sedbasin[j])):
            if data_sedbasin[j][k][0] == units_sorted[i][0]:
                units_sorted[i].append([data_sedbasin[j][k][2], data_sedbasin[j][k][3]])
for i in range(len(colors)):
    for j in range(len(units_sorted)):
        if units_sorted[j][0] == colors[i][0] and colors[i][1] not in units_sorted[j]:
            units_sorted[j].append(colors[i][1])
print("units_sorted:", units_sorted)
legend_lithology = []
n_units = []
for i in range(len(units_sorted)):
    n_units.append(len(units_sorted[i])-2)
    legend_lithology.append(mpatches.Patch(facecolor=units_sorted[i][-1], hatch="", label=units_sorted[i][0]))
print("")
for i in range(len(n_units)):
    print(units_sorted[i][0], units_sorted[i][-1], n_units[i])
print("sum:", sum(n_units))

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vP_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], alpha=0.9)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$v_P$")
plt.legend(fontsize="x-small", framealpha=1.0)
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vS_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], alpha=0.9)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$v_S$")
plt.legend(fontsize="x-small", framealpha=1.0)
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vPvS_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], alpha=0.9)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$v_P/v_S$")
plt.legend(fontsize="x-small", framealpha=1.0)
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], phi_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], alpha=0.9)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("$\\varphi$")
plt.legend(fontsize="x-small", framealpha=1.0)
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], gr_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], alpha=0.9)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("GR")
plt.legend(fontsize="x-small", framealpha=1.0)
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], pe_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], alpha=0.9)
plt.grid()
plt.xlabel("$\\varrho$")
plt.ylabel("PE")
plt.legend(fontsize="x-small", framealpha=1.0)
plt.show()

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey='row', gridspec_kw={'wspace': 0.15}, figsize=(9, 12))
fig.subplots_adjust(wspace=0.25)
# 1
ax1.plot(gr, top, color="#00549F", linewidth=2)
ax1.set_xlabel("GR (API)")
ax1.set_ylabel("Depth (m)")
#ax1.set_xscale('log')
ax1.set_xlim(0, 300)
ax1.set_xticks(np.arange(0, 400, 100))
ax1.set_ylim(0, 500)
ax1.set_yticks(np.arange(0, 500+50, 50))
ax1.grid(True)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 2
vP_edit = [vP[i]/1000 for i in range(len(vP))]
vS_edit = [vS[i]/1000 for i in range(len(vS))]
ax2.plot(vP_edit, top, color="#00549F", linewidth=2)
ax2.set_xlabel("vP (km/s)")
ax2.set_xlim(0, 8.5)
ax2.set_xticks(np.arange(0, 8.5, 2.0))
ax2.set_ylim(0, 500)
ax2.set_yticks(np.arange(0, 500+50, 50))
ax2_2 = ax2.twiny()
ax2_2.plot(vS_edit, top, color="#CC071E", linewidth=2)
ax2_2.set_xlabel("vS (km/s)")
ax2_2.set_xlim(0, 8.5)
ax2_2.set_xticks(np.arange(0, 8.5, 2.0))
ax2_2.minorticks_on()
ax2_2.grid(True)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 3
phi_edit = [phi[i]*100 for i in range(len(phi))]
ax3.plot(rho, top, color="#57AB27", linewidth=2)
ax3.set_xlabel("$\\varrho$ [g/cm^3]")
#ax3.set_xlim(1.7, 3.3)
ax3.set_xlim(1.6, 3.2)
#ax3.set_xticks(np.around(np.linspace(1.7, 3.3, 4, endpoint=True), decimals=1))
ax3.set_xticks(np.around(np.linspace(1.6, 3.2, 4, endpoint=True), decimals=1))
ax3.set_ylim(0, 500)
ax3.set_yticks(np.arange(0, 500+50, 50))
ax3_2 = ax3.twiny()
ax3_2.plot(phi_edit, top, color="#00549F", linewidth=2)
ax3_2.set_xlabel("$\\varphi$ [1]")
#ax3_2.set_xlim(57, -27)
ax3_2.set_xlim(60, 0)
#ax3_2.set_xticks(np.around(np.linspace(57, -27, 6, endpoint=True), decimals=0))
ax3_2.set_xticks(np.around(np.linspace(60, 0, 6, endpoint=True), decimals=0))
ax3_2.minorticks_on()
ax3_2.grid(True)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 4
ax4.plot(pe, top, color="#00549F", linewidth=2)
ax4.set_xlabel("PE (barns/electron)")
ax4.set_xscale('log')
#ax4.set_xlim(0, 40)
#ax4.set_xticks(np.arange(0, max(gr)+100, 100))
ax4.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax4.get_xaxis().set_minor_formatter(mpl.ticker.ScalarFormatter())
ax4.xaxis.set_minor_formatter(NullFormatter())
ax4.set_ylim(0, 500)
ax4.set_yticks(np.arange(0, 500+50, 50))
ax4.grid(True)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 5
for i in range(len(n_units)):
    for j in range(1, n_units[i]+1):
        ax5.hist(np.linspace(units_sorted[i][j][0], units_sorted[i][j][1]), bins=len(n_units), color=units_sorted[i][-1], orientation="horizontal")
ax5.set_xlabel("Lithology")
ax5.set_xlim(0, 5)
ax5.set_xticks([])
ax5.set_ylim(0, 500)
ax5.set_yticks(np.arange(0, 500+50, 50))
ax5.margins(0.3, 0.0)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
ax5.legend(handles=legend_lithology, loc="upper right", bbox_to_anchor=(2.0, 1.0), shadow=True, ncol=1, prop={'size': 8}, frameon=False)
plt.tight_layout()
plt.savefig("Test_Stratigraphy_01.png", bbox_inches="tight")
plt.show()