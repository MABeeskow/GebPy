#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_sequences.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		03.04.2021

# -----------------------------------------------

## MODULES
import sys
import numpy as np
import random as rd
import  matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from modules import sequences, geophysics

## TESTING
# Test sedimentary basin generation within SedimentaryBasin class
max_thickness = 500
n_units = int(max_thickness/25)
data = sequences.SedimentaryBasin()
data_sedbasin = data.create_sedimentary_basin(maximum_thickness=max_thickness, n_units=n_units, csv_stratigraphy=True,
                                              csv_lithology=True, excludeRocksalt=True, excludeLimestone=True)
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
z = []
kmod = []
gmod = []
rock_sorted = []
rho_sorted = []
vP_sorted = []
vS_sorted = []
vPvS_sorted = []
phi_sorted = []
gr_sorted = []
pe_sorted = []
z_sorted = []
kmod_sorted = []
gmod_sorted = []
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
        z_sorted.append([data_sedbasin[i][0][0], []])
        kmod_sorted.append([data_sedbasin[i][0][0], []])
        gmod_sorted.append([data_sedbasin[i][0][0], []])
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
                z_sorted[k][1].append(data_sedbasin[i][j][4][1][0]*data_sedbasin[i][j][4][3][0])
                kmod_sorted[k][1].append(data_sedbasin[i][j][4][2][0])
                gmod_sorted[k][1].append(data_sedbasin[i][j][4][2][1])
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
        z.append(data_sedbasin[i][j][4][1][0]*data_sedbasin[i][j][4][3][0])
        kmod.append(data_sedbasin[i][j][4][2][0])
        gmod.append(data_sedbasin[i][j][4][2][1])
    print(data_sedbasin[i])

seismic_trace = geophysics.Seismology().create_seismic_trace(data_all=data_sedbasin)

colors = [["soil", "peru"], ["sand", "moccasin"], ["sandstone", "tan"], ["limestone", "lightblue"], ["shale", "olivedrab"], ["rock salt", "violet"], ["granite", "darkorange"], ["basalt", "grey"]]
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
    print(units_sorted[i][0], n_units[i])
print("sum:", sum(n_units), n_units)
print(legend_lithology)

## PLOTTING
fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vP_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], edgecolors="black", alpha=0.9)
plt.grid(color="grey", linestyle="dashed")
plt.xlabel("$\\varrho$ [g/cm$^3$]")
plt.ylabel("$v_P$ [m/s]")
plt.legend(loc="best", fontsize="x-small", framealpha=1.0)
ax.set_axisbelow(True)
plt.savefig("Test_01.png", bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vS_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], edgecolors="black", alpha=0.9)
plt.grid(color="grey", linestyle="dashed")
plt.xlabel("$\\varrho$ [g/cm$^3$]")
plt.ylabel("$v_S$ [m/s]")
plt.legend(loc="best", fontsize="x-small", framealpha=1.0)
ax.set_axisbelow(True)
plt.savefig("Test_02.png", bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], vPvS_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], edgecolors="black", alpha=0.9)
plt.grid(color="grey", linestyle="dashed")
plt.xlabel("$\\varrho$ [g/cm$^3$]")
plt.ylabel("$v_P/v_S$ [1]")
plt.legend(loc="best", fontsize="x-small", framealpha=1.0)
ax.set_axisbelow(True)
plt.savefig("Test_03.png", bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], phi_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], edgecolors="black", alpha=0.9)
plt.grid(color="grey", linestyle="dashed")
plt.xlabel("$\\varrho$ [g/cm$^3$]")
plt.ylabel("$\\varphi$ [1]")
plt.legend(loc="best", fontsize="x-small", framealpha=1.0)
ax.set_axisbelow(True)
plt.savefig("Test_04.png", bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], gr_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], edgecolors="black", alpha=0.9)
plt.grid(color="grey", linestyle="dashed")
plt.xlabel("$\\varrho$ [g/cm$^3$]")
plt.ylabel("GR [API]")
plt.legend(loc="best", fontsize="x-small", framealpha=1.0)
ax.set_axisbelow(True)
plt.savefig("Test_05.png", bbox_inches="tight")
plt.show()

fig, ax = plt.subplots(1, 1, dpi=100)
for i in range(len(rock_sorted)):
    plt.scatter(rho_sorted[i][1], pe_sorted[i][1], label=units_sorted[i][0], color=units_sorted[i][-1], edgecolors="black", alpha=0.9)
plt.grid(color="grey", linestyle="dashed")
plt.xlabel("$\\varrho$ [g/cm$^3$]")
plt.ylabel("PE [barns/electron]")
plt.legend(loc="best", fontsize="x-small", framealpha=1.0)
ax.set_axisbelow(True)
plt.savefig("Test_06.png", bbox_inches="tight")
plt.show()

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey='row', gridspec_kw={'wspace': 0.15}, figsize=(9, 12))
fig.subplots_adjust(wspace=0.25)
# 1
ax1.plot(gr, top, color="#00549F", linewidth=2)
ax1.set_xlabel("GR [API]")
ax1.set_ylabel("Depth [m]")
ax1.set_xlim(0, 200)
ax1.set_xticks(np.arange(0, 250, 50))
ax1.set_ylim(0, max_thickness)
ax1.set_yticks(np.arange(0, max_thickness+50, 50))
ax1.grid(color="grey", linestyle="dashed")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 2
vP_edit = [vP[i]/1000 for i in range(len(vP))]
vS_edit = [vS[i]/1000 for i in range(len(vS))]
ax2.plot(vP_edit, top, color="#00549F", linewidth=2)
ax2.set_xlabel("$v_P$ [km/s]")
ax2.set_xlim(0, 8.5)
ax2.set_xticks(np.arange(0, 8.5, 2.0))
ax2.xaxis.label.set_color("#00549F")
ax2.set_ylim(0, max_thickness)
ax2.set_yticks(np.arange(0, max_thickness+50, 50))
ax2.grid(color="grey", linestyle="dashed")
ax2_2 = ax2.twiny()
ax2_2.plot(vS_edit, top, color="#CC071E", linewidth=2)
ax2_2.set_xlabel("$v_S$ [km/s]")
ax2_2.set_xlim(0, 8.5)
ax2_2.set_xticks(np.arange(0, 8.5, 2.0))
ax2_2.minorticks_on()
ax2_2.xaxis.label.set_color("#CC071E")
ax2_2.grid(color="grey", linestyle="dashed")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 3
phi_edit = [phi[i]*100 for i in range(len(phi))]
ax3.plot(rho, top, color="#57AB27", linewidth=2)
ax3.set_xlabel("$\\varrho$ [g/cm$^3$]")
ax3.set_xlim(1.6, 3.2)
ax3.set_xticks(np.around(np.linspace(1.6, 3.2, 4, endpoint=True), decimals=1))
ax3.xaxis.label.set_color("#57AB27")
ax3.set_ylim(0, max_thickness)
ax3.set_yticks(np.arange(0, max_thickness+50, 50))
ax3.grid(color="grey", linestyle="dashed")
ax3_2 = ax3.twiny()
ax3_2.plot(phi_edit, top, color="#00549F", linewidth=2)
ax3_2.set_xlabel("$\\varphi$ [1]")
ax3_2.set_xlim(60, 0)
ax3_2.set_xticks(np.around(np.linspace(60, 0, 6, endpoint=True), decimals=0))
ax3_2.minorticks_on()
ax3_2.xaxis.label.set_color("#00549F")
ax3_2.grid(color="grey", linestyle="dashed")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 4
ax4.plot(pe, top, color="#00549F", linewidth=2)
ax4.set_xlabel("PE [barns/electron]")
ax4.set_xscale("log")
ax4.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax4.get_xaxis().set_minor_formatter(mpl.ticker.ScalarFormatter())
ax4.xaxis.set_minor_formatter(NullFormatter())
ax4.set_ylim(0, max_thickness)
ax4.set_yticks(np.arange(0, max_thickness+50, 50))
ax4.grid(color="grey", linestyle="dashed", which="both")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 5
for i in range(len(n_units)):
    for j in range(1, n_units[i]+1):
        ax5.hist(np.linspace(units_sorted[i][j][0], units_sorted[i][j][1]), bins=len(n_units), color=units_sorted[i][-1], orientation="horizontal")
ax5.set_xlabel("Lithology")
ax5.set_xlim(0, 5)
ax5.set_xticks([])
ax5.set_ylim(0, max_thickness)
ax5.set_yticks(np.arange(0, max_thickness+50, 50))
ax5.margins(0.3, 0.0)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
ax5.legend(handles=legend_lithology, loc="upper right", bbox_to_anchor=(2.0, 1.0), shadow=True, ncol=1, prop={'size': 8}, frameon=False)
plt.tight_layout()
plt.savefig("Test_Stratigraphy_01.png", bbox_inches="tight")
plt.show()

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey='row', gridspec_kw={'wspace': 0.15}, figsize=(9, 12))
fig.subplots_adjust(wspace=0.25)
# 1
ax1.plot(np.array(z)/1000, top, color="#00549F", linewidth=2)
ax1.set_xlabel("Z [kNs/m$^3$]")
ax1.set_ylabel("Depth [m]")
#ax1.set_xlim(0, 300)
#ax1.set_xticks(np.arange(0, 400, 100))
ax1.set_ylim(0, max_thickness)
ax1.set_yticks(np.arange(0, max_thickness+50, 50))
ax1.grid(color="grey", linestyle="dashed")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 2
ax2.plot(vP_edit, top, color="#00549F", linewidth=2)
ax2.set_xlabel("$v_P$ [km/s]")
ax2.set_xlim(0, 8.5)
ax2.set_xticks(np.arange(0, 8.5, 2.0))
ax2.xaxis.label.set_color("#00549F")
ax2.set_ylim(0, max_thickness)
ax2.set_yticks(np.arange(0, max_thickness+50, 50))
ax2.grid(color="grey", linestyle="dashed")
ax2_2 = ax2.twiny()
ax2_2.plot(vS_edit, top, color="#CC071E", linewidth=2)
ax2_2.set_xlabel("$v_S$ [km/s]")
ax2_2.set_xlim(0, 8.5)
ax2_2.set_xticks(np.arange(0, 8.5, 2.0))
ax2_2.minorticks_on()
ax2_2.xaxis.label.set_color("#CC071E")
ax2_2.grid(color="grey", linestyle="dashed")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 3
ax3.plot(kmod, top, color="#57AB27", linewidth=2)
ax3.set_xlabel("K [GPa]")
ax3.set_xlim(0, 60)
ax3.set_xticks(np.arange(0, 60, 20))
ax3.xaxis.label.set_color("#57AB27")
ax3.set_ylim(0, max_thickness)
ax3.set_yticks(np.arange(0, max_thickness+50, 50))
ax3.grid(color="grey", linestyle="dashed")
ax3_2 = ax3.twiny()
ax3_2.plot(gmod, top, color="#00549F", linewidth=2)
ax3_2.set_xlabel("G [GPa]")
ax3_2.set_xlim(0, 60)
ax3_2.set_xticks(np.arange(0, 60, 20))
ax3_2.minorticks_on()
ax3_2.xaxis.label.set_color("#00549F")
ax3_2.grid(color="grey", linestyle="dashed")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 4
ax4.plot(seismic_trace, top, color="#000000", linewidth=2)
ax4.fill_betweenx(top, 0.0, seismic_trace, where=(seismic_trace>0.0), color="#00549F")
ax4.fill_betweenx(top, 0.0, seismic_trace, where=(seismic_trace<0.0), color="#CC071E")
ax4.set_xlabel("Seismic trace")
ax4.set_ylim(0, max_thickness)
ax4.set_yticks(np.arange(0, max_thickness+50, 50))
ax4.grid(color="grey", linestyle="dashed", which="both")
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
# 5
for i in range(len(n_units)):
    for j in range(1, n_units[i]+1):
        ax5.hist(np.linspace(units_sorted[i][j][0], units_sorted[i][j][1]), bins=len(n_units), color=units_sorted[i][-1], orientation="horizontal")
ax5.set_xlabel("Lithology")
ax5.set_xlim(0, 5)
ax5.set_xticks([])
ax5.set_ylim(0, max_thickness)
ax5.set_yticks(np.arange(0, max_thickness+50, 50))
ax5.margins(0.3, 0.0)
plt.gca().invert_yaxis()
plt.rc('axes', axisbelow=True)
ax5.legend(handles=legend_lithology, loc="upper right", bbox_to_anchor=(2.0, 1.0), shadow=True, ncol=1, prop={'size': 8}, frameon=False)
plt.tight_layout()
plt.savefig("Test_Stratigraphy_02.png", bbox_inches="tight")
plt.show()