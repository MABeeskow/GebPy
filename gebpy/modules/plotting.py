#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		plotting.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		22.05.2020

#-----------------------------------------------

## MODULES
import numpy as np
import matplotlib.pyplot as plt

class scatter:
    #
    def __init__(self):
        pass
    #
    def scatterPlotLith(self, sequences, x="VP", y="RHOB", z="GR", lith="sandstone"):
        self.sequences = sequences
        self.paramX = x
        self.paramY = y
        self.paramZ = z
        self.lith = lith
        #
        x = []
        y = []
        z = []
        for i in range(0, len(self.sequences)):
            if self.sequences[i][0] == self.lith:
                if self.paramX == "RHOB":
                    x.append(self.sequences[i][4])
                elif self.paramX == "GR":
                    x.append(self.sequences[i][6])
                elif self.paramX == "NPHI":
                    x.append(self.sequences[i][7])
                elif self.paramX == "POISSON":
                    x.append(self.sequences[i][10])
                elif self.paramX == "VP":
                    x.append(self.sequences[i][5][0])
                elif self.paramX == "VS":
                    x.append(self.sequences[i][5][1])
                elif self.paramX == "VF":
                    x.append(self.sequences[i][5][3])
                elif self.paramX == "Qz":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Qz":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Cal":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Cal":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Dol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Dol":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Or":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Or":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Bt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Bt":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Glt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Glt":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Chl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Chl":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Ilt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ilt":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Anh":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Anh":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Hl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Hl":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Gp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gp":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Syl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Syl":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Gn":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gn":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Ccp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ccp":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Mol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Mol":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Py":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Py":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Ab":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ab":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "An":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "An":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Ms":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ms":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Wo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Wo":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "En":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "En":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Fs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fs":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Tr":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tr":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Fo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fo":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Fa":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fa":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Tep":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tep":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Kln":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kln":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Kfs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kfs":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Afs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Afs":
                            x.append(self.sequences[i][9][j][1])
                elif self.paramX == "Pl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Pl":
                            x.append(self.sequences[i][9][j][1])
        for i in range(0, len(self.sequences)):
            if self.sequences[i][0] == self.lith:
                if self.paramY == "RHOB":
                    y.append(self.sequences[i][4])
                elif self.paramY == "GR":
                    y.append(self.sequences[i][6])
                elif self.paramY == "NPHI":
                    y.append(self.sequences[i][7])
                elif self.paramY == "POISSON":
                    y.append(self.sequences[i][10])
                elif self.paramY == "VP":
                    y.append(self.sequences[i][5][0])
                elif self.paramY == "VS":
                    y.append(self.sequences[i][5][1])
                elif self.paramY == "VF":
                    y.append(self.sequences[i][5][3])
                elif self.paramY == "Qz":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Qz":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Cal":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Cal":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Dol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Dol":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Or":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Or":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Bt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Bt":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Glt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Glt":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Chl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Chl":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Ilt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ilt":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Anh":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Anh":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Hl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Hl":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Gp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gp":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Syl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Syl":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Gn":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gn":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Ccp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ccp":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Mol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Mol":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Py":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Py":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Ab":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ab":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "An":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "An":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Ms":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ms":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Wo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Wo":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "En":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "En":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Fs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fs":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Tr":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tr":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Fo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fo":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Fa":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fa":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Tep":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tep":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Kln":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kln":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Kfs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kfs":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Afs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Afs":
                            y.append(self.sequences[i][9][j][1])
                elif self.paramY == "Pl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Pl":
                            y.append(self.sequences[i][9][j][1])
        for i in range(0, len(self.sequences)):
            if self.sequences[i][0] == self.lith:
                if self.paramZ == "RHOB":
                    z.append(self.sequences[i][4])
                elif self.paramZ == "GR":
                    z.append(self.sequences[i][6])
                elif self.paramZ == "NPHI":
                    z.append(self.sequences[i][7])
                elif self.paramZ == "POISSON":
                    z.append(self.sequences[i][10])
                elif self.paramZ == "VP":
                    z.append(self.sequences[i][5][0])
                elif self.paramZ == "VS":
                    z.append(self.sequences[i][5][1])
                elif self.paramZ == "VF":
                    z.append(self.sequences[i][5][3])
                elif self.paramZ == "Qz":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Qz":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Cal":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Cal":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Dol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Dol":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Or":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Or":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Bt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Bt":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Glt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Glt":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Chl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Chl":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Ilt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ilt":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Anh":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Anh":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Hl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Hl":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Gp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gp":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Syl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Syl":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Gn":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gn":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Ccp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ccp":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Mol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Mol":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Py":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Py":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Ab":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ab":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "An":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "An":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Ms":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ms":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Wo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Wo":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "En":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "En":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Fs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fs":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Tr":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tr":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Fo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fo":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Fa":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fa":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Tep":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tep":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Kln":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kln":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Kfs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kfs":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Afs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Afs":
                            z.append(self.sequences[i][9][j][1])
                elif self.paramZ == "Pl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Pl":
                            z.append(self.sequences[i][9][j][1])
            else:
                pass
        #
        if len(x) > 0:
            plt.figure()
            plot = plt.scatter(x, y, c=z, cmap="viridis", norm=None, s=30, alpha=0.9)
            cbar = plt.colorbar(plot)
            cbar.set_label(self.paramZ)
            plt.title(self.lith)
            plt.xlabel(self.paramX)
            plt.ylabel(self.paramY)
            plt.grid(color="#000000", linestyle='dashed', alpha=0.75)
            plt.rc('axes', axisbelow=True)
            plt.show()
        else:
            pass
    #
class histogram:
    #
    def __init__(self):
        pass
    #
    def histPlotLith(self, sequences, x="VP", lith="sandstone"):
        self.sequences = sequences
        self.x = x
        self.lith = lith
        #
        data = []
        #
        for i in range(0, len(self.sequences)):
            if self.sequences[i][0] == self.lith:
                if self.x == "RHOB":
                    data.append(self.sequences[i][4])
                elif self.x == "GR":
                    data.append(self.sequences[i][6])
                elif self.x == "NPHI":
                    data.append(self.sequences[i][7])
                elif self.x == "POISSON":
                    data.append(self.sequences[i][10])
                elif self.x == "VP":
                    data.append(self.sequences[i][5][0])
                elif self.x == "VS":
                    data.append(self.sequences[i][5][1])
                elif self.x == "VF":
                    data.append(self.sequences[i][5][3])
                elif self.x == "Qz":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Qz":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Cal":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Cal":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Dol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Dol":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Or":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Or":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Bt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Bt":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Glt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Glt":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Chl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Chl":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Ilt":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ilt":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Anh":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Anh":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Hl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Hl":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Gp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gp":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Syl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Syl":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Gn":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Gn":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Ccp":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ccp":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Mol":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Mol":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Py":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Py":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Ab":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ab":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "An":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "An":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Ms":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Ms":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Wo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Wo":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "En":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "En":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Fs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fs":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Tr":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tr":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Fo":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fo":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Fa":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Fa":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Tep":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Tep":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Kln":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kln":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Kfs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Kfs":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Afs":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Afs":
                            data.append(self.sequences[i][9][j][1])
                elif self.x == "Pl":
                    for j in range(0, len(self.sequences[i][9])):
                        if self.sequences[i][9][j][0] == "Pl":
                            data.append(self.sequences[i][9][j][1])
            else:
                pass
        #
        if len(data) > 0:
            plt.figure()
            plot = plt.hist(data, bins=20, color="#00549F", edgecolor="black", align="left", density=1)
            plt.title(self.lith)
            plt.xlabel(self.x)
            plt.ylabel("Probability")
            plt.grid(color="#000000", linestyle='dashed', alpha=0.75)
            plt.rc('axes', axisbelow=True)
            plt.show()
        else:
            pass
    #