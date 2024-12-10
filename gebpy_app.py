#!/usr/bin/env python
# -*-coding: utf-8 -*-

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# Name:		gebpy_app.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		05.07.2024
# License:  GPL v3.0

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

## MODULES
# external
import os, sys
from datetime import datetime
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import collections
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from modules.geophysics import Elasticity as elast
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.patches as mpatches
from matplotlib.figure import Figure
# internal
from modules.gui_elements import SimpleElements
from modules.oxides import Oxides
from modules.carbonates import Carbonates, CarbonateRocks
from modules.igneous import Plutonic, Volcanic, UltraMafic, Pyroclastic
from modules.sulfides import Sulfides
from modules.sulfates import Sulfates
from modules.halides import Halides
from modules.phospides import Phospides
from modules.phosphates import Phosphates
from modules.silicates import Phyllosilicates, Tectosilicates, Inosilicates, Nesosilicates, Sorosilicates, \
    Cyclosilicates
from modules.organics import Organics
from modules.fluids import Water
from modules.siliciclastics import SiliciclasticRocks, Geophysics
from modules.ore import OreRocks
from modules.metamorphics import GranuliteFacies, GreenschistFacies, AmphiboliteFacies
from modules.exploration import ExplorationInterface
# Sequence Stratigraphy
from modules.series import Muschelkalk, Zechstein, Buntsandstein
from modules.petrophysics import SeismicVelocities

## GUI
class GebPyGUI(tk.Frame):
    #
    def __init__(self, parent, var_screen_width, var_screen_height):
        tk.Frame.__init__(self, parent)

        self.path_gebpy = os.path.dirname(os.path.realpath(sys.argv[0]))

        var_screen_width = var_screen_width
        var_screen_height = var_screen_height

        self.str_version_number = "0.9.0"
        self.val_version = "GebPy: " + self.str_version_number + " - 10.12.2024"

        ### Container
        self.gui_elements = {}
        self.gui_elements_sub = {}
        gui_elements = ["Frame", "Label", "Button", "Radiobutton", "Checkbox", "Entry", "Option Menu", "Canvas",
                        "Figure", "Axis"]
        gui_priority = ["Static", "Temporary", "Rockbuilder Static", "Rockbuilder Temporary", "Element Concentration",
                        "Compound Concentration"]
        gui_subwindows = ["Trace Elements", "Mineralogy", "Mineral Physics", "Mineral Chemistry", "Synthetic LA-ICP-MS"]
        #
        for subwindow in gui_subwindows:
            self.gui_elements_sub[subwindow] = {}
            for priority in gui_priority:
                self.gui_elements_sub[subwindow][priority] = {}
                for gui_element in gui_elements:
                    if gui_element not in ["Canvas", "Figure", "Axis"]:
                        self.gui_elements_sub[subwindow][priority][gui_element] = []
                    else:
                        self.gui_elements_sub[subwindow][priority][gui_element] = {}
        #
        for priority in gui_priority:
            self.gui_elements[priority] = {}
            for gui_element in gui_elements:
                if gui_element not in ["Canvas", "Figure", "Axis"]:
                    self.gui_elements[priority][gui_element] = []
                else:
                    self.gui_elements[priority][gui_element] = {}
        #
        ### Colors
        self.colors_gebpy = {
            "Background": "#EFEFEF", "Navigation": "#252422", "Accent": "#EB5E28", "Option": "#CCC5B9",
            "White": "#FFFFFF", "Black": "#000000", "Accent Blue": "#118AB2"}
        #
        ### Variables
        self.gui_variables = {}
        for gui_element in gui_elements:
            self.gui_variables[gui_element] = {}
        #
        self.last_rb_analysis_mineral = tk.IntVar()
        self.last_rb_analysis_mineral.set(42)
        self.last_rb_diagram_mineral = tk.IntVar()
        self.last_rb_diagram_mineral.set(42)
        #
        self.last_rb_analysis_rock = tk.IntVar()
        self.last_rb_analysis_rock.set(42)
        self.last_rb_diagram_rock = tk.IntVar()
        self.last_rb_diagram_rock.set(42)
        #
        self.last_rb_setting = {"Mineralogy": {}, "Petrology": {}, "Stratigraphy": {}}
        # Mineralogy
        self.last_rb_setting["Mineralogy"]["General Mode"] = tk.IntVar()
        self.last_rb_setting["Mineralogy"]["General Mode"].set(42)
        self.last_rb_setting["Mineralogy"]["Mineral Physics"] = tk.IntVar()
        self.last_rb_setting["Mineralogy"]["Mineral Physics"].set(42)
        self.last_rb_setting["Mineralogy"]["Mineral Chemistry"] = tk.IntVar()
        self.last_rb_setting["Mineralogy"]["Mineral Chemistry"].set(42)
        self.last_rb_setting["Mineralogy"]["Composition Elements"] = tk.IntVar()
        self.last_rb_setting["Mineralogy"]["Composition Elements"].set(42)
        # Petrology
        self.last_rb_setting["Petrology"]["General Mode"] = tk.IntVar()
        self.last_rb_setting["Petrology"]["General Mode"].set(42)
        self.last_rb_setting["Petrology"]["Rock Physics"] = tk.IntVar()
        self.last_rb_setting["Petrology"]["Rock Physics"].set(42)
        self.last_rb_setting["Petrology"]["Rock Chemistry"] = tk.IntVar()
        self.last_rb_setting["Petrology"]["Rock Chemistry"].set(42)
        self.last_rb_setting["Petrology"]["Composition Minerals"] = tk.IntVar()
        self.last_rb_setting["Petrology"]["Composition Minerals"].set(42)
        self.last_rb_setting["Petrology"]["Composition Elements"] = tk.IntVar()
        self.last_rb_setting["Petrology"]["Composition Elements"].set(42)

        self.container_variables = {"Radiobuttons": {"Category Diagram": tk.IntVar(), "Category Data": tk.IntVar()}}
        self.container_variables["Radiobuttons"]["Category Diagram"].set(0)
        self.container_variables["Radiobuttons"]["Category Data"].set(0)

        self.var_rb = {
            "Rock analysis": {"Rock definition": tk.IntVar(), "External data": tk.IntVar(), "x-axis": tk.StringVar(),
                              "y-axis": tk.StringVar()}}

        self.chemistry_data = {
            "O": 15.999, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "K": 39.098, "Ca": 40.078,
            "Ti": 47.867, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Ga": 69.723, "Ge": 72.630, "Zr": 91.224,
            "Ba": 137.33, "B": 10.81, "Ag": 107.87, "As": 74.922, "Li": 6.94, "Rb": 85.468, "Cs": 132.91, "Sr": 87.62,
            "Sc": 44.956, "Y": 88.906, "Hf": 178.49, "V": 50.942, "Nb": 92.906, "Ta": 180.95, "Mo": 95.962, "W": 183.84,
            "Tc": 98.906, "Re": 186.21, "Ru": 101.07, "Os": 190.23, "Co": 58.933, "Rh": 102.91, "Ir": 192.22,
            "Ni": 58.693, "Pd": 106.42, "Pt": 195.08, "Cu": 63.546, "Au": 196.97, "Zn": 65.38, "Cd": 112.41,
            "Hg": 200.59, "In": 114.82, "Tl": 204.38, "C": 12.011, "Sn": 118.71, "Pb": 207.2, "N": 14.007, "Sb": 121.76,
            "Bi": 208.98, "S": 32.06, "Se": 78.96, "Te": 127.60, "Po": 209.98, "Cl": 35.45, "Br": 79.904, "I": 126.90,
            "At": 210.99, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 146.92, "Sm": 150.36,
            "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93,
            "Yb": 173.05, "Lu": 174.97, "Ac": 227.03, "Th": 232.04, "Pa": 231.04, "U": 238.05, "Be": 9.0122,
            "F": 18.998, "H": 1.008}
        self.chemistry_data_sills = {
            "O": 16.000, "Na": 22.990, "Mg": 24.300, "Al": 26.980, "Si": 28.090, "P": 30.970, "K": 39.100, "Ca": 40.080,
            "Ti": 47.870, "Cr": 52.000, "Mn": 54.940, "Fe": 55.850, "Ga": 69.720, "Ge": 72.610, "Zr": 91.220,
            "Ba": 137.300}
        self.chemistry_data_oxides = {
            "SiO2": 60.083, "Al2O3": 101.961, "Fe2O3": 159.687, "FeO": 71.844, "Na2O": 61.979, "TiO2": 79.865,
            "MnO": 70.937, "Mn2O3": 157.873, "SnO": 134.709, "Li2O": 29.879, "Ga2O3": 187.443, "B2O3": 69.617,
            "BeO": 25.0112, "GeO2": 104.628, "CaO": 56.077, "Rb2O": 186.935, "AgO": 123.869, "As2O3": 197.841,
            "Au2O": 409.939, "BaO": 153.32, "Br2O": 175.807, "Cl2O": 86.899, "Cs2O": 281.819, "CuO": 79.545,
            "PbO": 223.199, "SO3": 80.057, "Sb2O3": 291.517, "SrO": 103.619, "WO3": 231.837, "ZnO": 81.379,
            "MgO": 40.304, "K2O": 55.097, "SnO2": 150.708, "Ag2O": 231.739, "Bi2O5": 497.955, "CO2": 44.009,
            "CdO": 128.409, "Ce2O3": 328.237, "CeO2": 172.118, "CoO": 74.932, "Cr2O3": 151.989, "Dy2O3": 372.997,
            "Er2O3": 382.517, "Eu2O3": 351.917, "Gd2O3": 362.497, "HfO2": 404.977, "HgO": 216.589, "Ho2O3": 377.857,
            "In2O3": 277.637, "IrO": 208.219, "La2O3": 325.817, "Lu2O3": 397.937, "MnO2": 86.936, "MoO3": 143.959,
            "N2O5": 108.009, "Nb2O5": 265.807, "Nd2O3": 336.477, "NiO": 74.692, "OsO": 206.229, "P2O5": 141.943,
            "PbO2": 239.198, "PdO": 122.419, "Pr2O3": 329.817, "Pr6O11": 1021.449, "PtO": 211.079, "ReO": 202.209,
            "RhO": 118.909, "RuO": 117.069, "SO4": 96.056, "Sb2O5": 323.515, "Sc2O3": 137.909, "SeO3": 126.957,
            "Sm2O3": 348.717, "Ta2O5": 441.895, "Tb2O3": 365.857, "Tb4O7": 747.713, "TeO3": 175.597, "ThO2": 264.038,
            "Tl2O3": 456.757, "Tm2O3": 385.857, "UO2": 270.048, "UO3": 286.047, "U3O8": 842.142, "V2O5": 181.879,
            "Y2O3": 225.809, "Yb2O3": 394.097, "ZrO2": 123.222, "I2O4": 317.796, "I2O5": 333.795, "I4O9": 651.591,
            "I2O": 269.799, "Ni2O3": 165.383, "Co2O3": 165.863, "CrO": 67.995, "H2O": 18.015, "SO2": 64.058}

        self.conversion_factors = {
            "SiO2": (self.chemistry_data["Si"]/self.chemistry_data_oxides["SiO2"])**(-1),
            "Al2O3": (2*self.chemistry_data["Al"]/self.chemistry_data_oxides["Al2O3"])**(-1),
            "Fe2O3": (2*self.chemistry_data["Fe"]/self.chemistry_data_oxides["Fe2O3"])**(-1),
            "FeO": (self.chemistry_data["Fe"]/self.chemistry_data_oxides["FeO"])**(-1),
            "Na2O": (2*self.chemistry_data["Na"]/self.chemistry_data_oxides["Na2O"])**(-1),
            "TiO2": (self.chemistry_data["Ti"]/self.chemistry_data_oxides["TiO2"])**(-1),
            "MnO": (self.chemistry_data["Mn"]/self.chemistry_data_oxides["MnO"])**(-1),
            "Mn2O3": (2*self.chemistry_data["Mn"]/self.chemistry_data_oxides["Mn2O3"])**(-1),
            "SnO": (self.chemistry_data["Sn"]/self.chemistry_data_oxides["SnO"])**(-1),
            "Li2O": (2*self.chemistry_data["Li"]/self.chemistry_data_oxides["Li2O"])**(-1),
            "Ga2O3": (2*self.chemistry_data["Ga"]/self.chemistry_data_oxides["Ga2O3"])**(-1),
            "B2O3": (2*self.chemistry_data["B"]/self.chemistry_data_oxides["B2O3"])**(-1),
            "BeO": (self.chemistry_data["Be"]/self.chemistry_data_oxides["BeO"])**(-1),
            "GeO2": (self.chemistry_data["Ge"]/self.chemistry_data_oxides["GeO2"])**(-1),
            "CaO": (self.chemistry_data["Ca"]/self.chemistry_data_oxides["CaO"])**(-1),
            "Rb2O": (2*self.chemistry_data["Rb"]/self.chemistry_data_oxides["Rb2O"])**(-1),
            "AgO": (self.chemistry_data["Ag"]/self.chemistry_data_oxides["AgO"])**(-1),
            "As2O3": (2*self.chemistry_data["As"]/self.chemistry_data_oxides["As2O3"])**(-1),
            "Au2O": (2*self.chemistry_data["Au"]/self.chemistry_data_oxides["Au2O"])**(-1),
            "BaO": (self.chemistry_data["Ba"]/self.chemistry_data_oxides["BaO"])**(-1),
            "Br2O": (2*self.chemistry_data["Br"]/self.chemistry_data_oxides["Br2O"])**(-1),
            "Cl2O": (2*self.chemistry_data["Cl"]/self.chemistry_data_oxides["Cl2O"])**(-1),
            "Cs2O": (2*self.chemistry_data["Cs"]/self.chemistry_data_oxides["Cs2O"])**(-1),
            "CuO": (self.chemistry_data["Cu"]/self.chemistry_data_oxides["CuO"])**(-1),
            "PbO": (self.chemistry_data["Pb"]/self.chemistry_data_oxides["PbO"])**(-1),
            "SO3": (self.chemistry_data["S"]/self.chemistry_data_oxides["SO3"])**(-1),
            "Sb2O3": (2*self.chemistry_data["Sb"]/self.chemistry_data_oxides["Sb2O3"])**(-1),
            "SrO": (self.chemistry_data["Sr"]/self.chemistry_data_oxides["SrO"])**(-1),
            "WO3": (self.chemistry_data["W"]/self.chemistry_data_oxides["WO3"])**(-1),
            "ZnO": (self.chemistry_data["Zn"]/self.chemistry_data_oxides["ZnO"])**(-1),
            "MgO": (self.chemistry_data["Mg"]/self.chemistry_data_oxides["MgO"])**(-1),
            "K2O": (2*self.chemistry_data["K"]/self.chemistry_data_oxides["K2O"])**(-1),
            "SnO2": (self.chemistry_data["Sn"]/self.chemistry_data_oxides["SnO2"])**(-1),
            "Ag2O": (2*self.chemistry_data["Ag"]/self.chemistry_data_oxides["Ag2O"])**(-1),
            "Bi2O5": (2*self.chemistry_data["Bi"]/self.chemistry_data_oxides["Bi2O5"])**(-1),
            "CO2": (self.chemistry_data["C"]/self.chemistry_data_oxides["CO2"])**(-1),
            "CdO": (self.chemistry_data["Cd"]/self.chemistry_data_oxides["CdO"])**(-1),
            "Ce2O3": (2*self.chemistry_data["Ce"]/self.chemistry_data_oxides["Ce2O3"])**(-1),
            "CeO2": (self.chemistry_data["Ce"]/self.chemistry_data_oxides["CeO2"])**(-1),
            "CoO": (self.chemistry_data["Co"]/self.chemistry_data_oxides["CoO"])**(-1),
            "CrO": (self.chemistry_data["Cr"]/self.chemistry_data_oxides["CrO"])**(-1),
            "Cr2O3": (2*self.chemistry_data["Cr"]/self.chemistry_data_oxides["Cr2O3"])**(-1),
            "Dy2O3": (2*self.chemistry_data["Dy"]/self.chemistry_data_oxides["Dy2O3"])**(-1),
            "Er2O3": (2*self.chemistry_data["Er"]/self.chemistry_data_oxides["Er2O3"])**(-1),
            "Eu2O3": (2*self.chemistry_data["Eu"]/self.chemistry_data_oxides["Eu2O3"])**(-1),
            "Gd2O3": (2*self.chemistry_data["Gd"]/self.chemistry_data_oxides["Gd2O3"])**(-1),
            "HfO2": (self.chemistry_data["Hf"]/self.chemistry_data_oxides["HfO2"])**(-1),
            "HgO": (self.chemistry_data["Hg"]/self.chemistry_data_oxides["HgO"])**(-1),
            "Ho2O3": (2*self.chemistry_data["Ho"]/self.chemistry_data_oxides["Ho2O3"])**(-1),
            "In2O3": (2*self.chemistry_data["In"]/self.chemistry_data_oxides["In2O3"])**(-1),
            "IrO": (self.chemistry_data["Ir"]/self.chemistry_data_oxides["IrO"])**(-1),
            "La2O3": (2*self.chemistry_data["La"]/self.chemistry_data_oxides["La2O3"])**(-1),
            "Lu2O3": (2*self.chemistry_data["Lu"]/self.chemistry_data_oxides["Lu2O3"])**(-1),
            "MnO2": (self.chemistry_data["Mn"]/self.chemistry_data_oxides["MnO2"])**(-1),
            "MoO3": (self.chemistry_data["Mo"]/self.chemistry_data_oxides["MoO3"])**(-1),
            "N2O5": (2*self.chemistry_data["N"]/self.chemistry_data_oxides["N2O5"])**(-1),
            "Nb2O5": (2*self.chemistry_data["Nb"]/self.chemistry_data_oxides["Nb2O5"])**(-1),
            "Nd2O3": (2*self.chemistry_data["Nd"]/self.chemistry_data_oxides["Nd2O3"])**(-1),
            "NiO": (self.chemistry_data["Ni"]/self.chemistry_data_oxides["NiO"])**(-1),
            "OsO": (self.chemistry_data["Os"]/self.chemistry_data_oxides["OsO"])**(-1),
            "P2O5": (2*self.chemistry_data["P"]/self.chemistry_data_oxides["P2O5"])**(-1),
            "PbO2": (self.chemistry_data["Pb"]/self.chemistry_data_oxides["PbO2"])**(-1),
            "PdO": (self.chemistry_data["Pd"]/self.chemistry_data_oxides["PdO"])**(-1),
            "Pr2O3": (2*self.chemistry_data["Pr"]/self.chemistry_data_oxides["Pr2O3"])**(-1),
            "Pr6O11": (6*self.chemistry_data["Pr"]/self.chemistry_data_oxides["Pr6O11"])**(-1),
            "PtO": (self.chemistry_data["Pt"]/self.chemistry_data_oxides["PtO"])**(-1),
            "ReO": (self.chemistry_data["Re"]/self.chemistry_data_oxides["ReO"])**(-1),
            "RhO": (self.chemistry_data["Rh"]/self.chemistry_data_oxides["RhO"])**(-1),
            "RuO": (self.chemistry_data["Ru"]/self.chemistry_data_oxides["RuO"])**(-1),
            "SO4": (self.chemistry_data["S"]/self.chemistry_data_oxides["SO4"])**(-1),
            "Sb2O5": (2*self.chemistry_data["Sb"]/self.chemistry_data_oxides["Sb2O5"])**(-1),
            "Sc2O3": (2*self.chemistry_data["Sc"]/self.chemistry_data_oxides["Sc2O3"])**(-1),
            "SeO3": (self.chemistry_data["Se"]/self.chemistry_data_oxides["SeO3"])**(-1),
            "Sm2O3": (2*self.chemistry_data["Sm"]/self.chemistry_data_oxides["Sm2O3"])**(-1),
            "Ta2O5": (2*self.chemistry_data["Ta"]/self.chemistry_data_oxides["Ta2O5"])**(-1),
            "Tb2O3": (2*self.chemistry_data["Tb"]/self.chemistry_data_oxides["Tb2O3"])**(-1),
            "Tb4O7": (4*self.chemistry_data["Tb"]/self.chemistry_data_oxides["Tb4O7"])**(-1),
            "TeO3": (self.chemistry_data["Te"]/self.chemistry_data_oxides["TeO3"])**(-1),
            "ThO2": (self.chemistry_data["Th"]/self.chemistry_data_oxides["ThO2"])**(-1),
            "Tl2O3": (2*self.chemistry_data["Tl"]/self.chemistry_data_oxides["Tl2O3"])**(-1),
            "Tm2O3": (2*self.chemistry_data["Tm"]/self.chemistry_data_oxides["Tm2O3"])**(-1),
            "UO2": (self.chemistry_data["U"]/self.chemistry_data_oxides["UO2"])**(-1),
            "UO3": (self.chemistry_data["U"]/self.chemistry_data_oxides["UO3"])**(-1),
            "U3O8": (3*self.chemistry_data["U"]/self.chemistry_data_oxides["U3O8"])**(-1),
            "V2O5": (2*self.chemistry_data["V"]/self.chemistry_data_oxides["V2O5"])**(-1),
            "Y2O3": (2*self.chemistry_data["Y"]/self.chemistry_data_oxides["Y2O3"])**(-1),
            "Yb2O3": (2*self.chemistry_data["Yb"]/self.chemistry_data_oxides["Yb2O3"])**(-1),
            "ZrO2": (self.chemistry_data["Zr"]/self.chemistry_data_oxides["ZrO2"])**(-1),
            "I2O4": (2*self.chemistry_data["I"]/self.chemistry_data_oxides["I2O4"])**(-1),
            "I2O5": (2*self.chemistry_data["I"]/self.chemistry_data_oxides["I2O5"])**(-1),
            "I4O9": (4*self.chemistry_data["I"]/self.chemistry_data_oxides["I4O9"])**(-1),
            "I2O": (2*self.chemistry_data["I"]/self.chemistry_data_oxides["I2O"])**(-1),
            "Co2O3": (2*self.chemistry_data["Co"]/self.chemistry_data_oxides["Co2O3"])**(-1),
            "Ni2O3": (2*self.chemistry_data["Ni"]/self.chemistry_data_oxides["Ni2O3"])**(-1),
            "H2O": (2*self.chemistry_data["H"]/self.chemistry_data_oxides["H2O"])**(-1),
            "SO2": (self.chemistry_data["S"]/self.chemistry_data_oxides["SO2"])**(-1)}

        self.chemistry_oxides_sorted = {
            "H": ["H2O"], "Li": ["Li2O"], "Be": ["BeO"], "B": ["B2O3"], "C": ["CO", "CO2"],
            "N": ["NO", "N2O3", "NO2", "N2O5"], "Na": ["Na2O"], "Mg": ["MgO"], "Al": ["Al2O3"], "Si": ["SiO2"],
            "P": ["P2O3", "P2O5"], "S": ["SO", "SO2", "SO3"], "Cl": ["Cl2O", "ClO2", "Cl2O3", "Cl2O5", "Cl2O7"],
            "K": ["K2O"], "Ca": ["CaO"], "Sc": ["Sc2O3"], "Ti": ["Ti2O3", "TiO2"], "V": ["VO", "V2O3", "VO2", "V2O5"],
            "Cr": ["CrO", "Cr2O3", "CrO3"], "Mn": ["MnO", "Mn2O3", "MnO2", "MnO3", "Mn2O7"],
            "Fe": ["FeO", "Fe2O3", "FeO3"], "Co": ["CoO", "Co2O3"], "Ni": ["NiO", "Ni2O3"], "Cu": ["Cu2O", "CuO"],
            "Zn": ["ZnO"], "Ga": ["Ga2O3"], "Ge": ["GeO2"], "As": ["As2O3", "As2O5"], "Se": ["SeO2", "SiO3"],
            "Br": ["Br2O", "Br2O3", "Br2O5", "Br2O7"], "Kr": ["KrO"], "Rb": ["Rb2O"], "Sr": ["SrO"], "Y": ["Y2O3"],
            "Zr": ["ZrO2"], "Nb": ["Nb2O3", "Nb2O5"], "Mo": ["MoO", "Mo2O3", "MoO2", "Mo2O5", "MoO3"], "Tc": ["Tc2O7"],
            "Ru": ["RuO", "Ru2O3", "RuO2", "RuO3", "RuO4"], "Rh": ["Rh2O", "RhO", "Rh2O3", "RhO2", "Rh2O5"],
            "Pd": ["PdO", "PdO2"], "Ag": ["Ag2O", "AgO"], "Cd": ["CdO"], "In": ["In2O3"], "Sn": ["SnO", "SnO2"],
            "Sb": ["Sb2O3", "Sb2O5"], "Te": ["TeO2", "TeO3"], "I": ["I2O", "I2O4", "I2O5", "I4O9"],
            "Xe": ["XeO", "XeO2", "XeO3"], "Cs": ["Cs2O"], "Ba": ["BaO"], "La": ["La2O3"], "Ce": ["Ce2O3", "CeO2"],
            "Pr": ["Pr2O3", "PrO2"], "Nd": ["Nd2O3"], "Pm": ["Pm2O3"], "Sm": ["SmO", "Sm2O3"], "Eu": ["EuO", "Eu2O3"],
            "Gd": ["Gd2O3"], "Tb": ["Tb2O3", "TbO2"], "Dy": ["Dy2O3"], "Ho": ["Ho2O3"], "Er": ["Er2O3"],
            "Tm": ["TmO", "Tm2O3"], "Yb": ["YbO", "Yb2O3"], "Lu": ["Lu2O3"], "Hf": ["HfO2"], "Ta": ["Ta2O5"],
            "W": ["WO", "WO2O3", "WO2", "W2O5", "WO3"], "Re": ["ReO", "ReO2", "ReO3", "Re2O7"],
            "Os": ["OsO", "Os2O3", "OsO2", "OsO3", "OsO4"], "Ir": ["Ir2O", "IrO", "Ir2O3", "IrO2", "IrO3"],
            "Pt": ["PtO", "PtO2"], "Au": ["Au2O", "Au2O3"], "Hg": ["Hg2O", "HgO"], "Tl": ["Tl2O", "Tl2O3"],
            "Pb": ["PbO", "PbO2"], "Bi": ["Bi2O3", "B2O5"], "Po": ["PoO", "PoO2", "PoO3"],
            "At": ["At2O", "At2O3", "At2O5", "At2O7"], "Rn": ["RnO"], "Fr": ["Fr2O"], "Ra": ["RaO"], "Ac": ["Ac2O3"],
            "Th": ["ThO2"], "Pa": ["PaO2", "Pa2O5"], "U": ["U2O3", "UO2", "U2O5", "UO3"],
            "Np": ["Np2O3", "NpO2", "Np2O5", "NpO3"], "Pu": ["Pu2O3", "PuO2", "Pu2O5", "PuO3"],
            "Am": ["Am2O3", "AmO2", "Am2O5", "AmO3"], "Cm": ["Cm2O3", "CmO2"], "Bk": ["Bk2O3", "BkO2"],
            "Cf": ["Cf2O3", "CfO2"], "Es": ["Es2O3"], "Fm": ["Fm2O3"], "Md": ["Md2O3"], "No": ["NoO", "No2O3"],
            "Lr": ["Lr2O3"]}


        # Stratigraphy
        #
        ### General Settings
        self.parent = parent
        self.parent.title("GebPy")
        #
        var_geometry = ""
        var_window_width = int(round(0.95*int(var_screen_width), -2))
        if var_window_width > 1800:
            var_window_width = 1800
        var_geometry += str(var_window_width)
        var_geometry += "x"
        var_window_height = int(round(0.85*int(var_screen_height), -2))
        if var_window_height > 975:
            var_window_height = 975
        var_geometry += str(var_window_height)
        var_geometry += "+0+0"
        self.parent.geometry(var_geometry)
        #self.parent.geometry("1800x975+0+0")
        self.parent.resizable(False, False)
        self.parent["bg"] = self.colors_gebpy["Background"]
        #
        ## Geometry and Layout
        # window_width = 1800
        # window_heigth = 975
        # #
        # row_min = 15
        # n_rows = int(window_heigth / row_min)
        # #
        # column_min = 10
        # n_columns = int(window_width / column_min)
        #
        window_width = var_window_width
        window_heigth = var_window_height
        #
        row_min = 15
        n_rows = int(window_heigth/row_min)
        #
        column_min = 10
        n_columns = int(window_width/column_min)
        #
        self.n_rows = n_rows
        self.n_columns = n_columns
        #
        for x in range(n_columns):
            tk.Grid.columnconfigure(self.parent, x, weight=1)
        for y in range(n_rows):
            tk.Grid.rowconfigure(self.parent, y, weight=1)
        #
        # Rows
        for i in range(0, n_rows):
            self.parent.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            self.parent.grid_columnconfigure(i, minsize=column_min)
        #
        ## Navigation Bar
        frm_navigation = SimpleElements(
            parent=self.parent, row_id=0, column_id=0, n_rows=100, n_columns=32, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Navigation"]).create_frame()
        frm_accent = SimpleElements(
            parent=self.parent, row_id=0, column_id=33, n_rows=100, n_columns=1, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Accent"]).create_frame()
        #
        self.gui_elements["Static"]["Frame"].extend([frm_navigation, frm_accent])
        #
        ## Logo
        try:
            gebpy_logo = tk.PhotoImage(file=self.path_gebpy + "/documents/readme_images/GebPy_Logo_new.png")
            gebpy_logo = gebpy_logo.subsample(5, 5)
            img = tk.Label(self.parent, image=gebpy_logo, bg=self.colors_gebpy["Navigation"])
            img.image = gebpy_logo
            img.grid(row=0, column=0, rowspan=5, columnspan=32, sticky="nesw")

            ## Icon
            #pysills_icon = tk.PhotoImage(file=self.path_pysills + str("/documentation/images/PySILLS_Icon.png"))
            gebpy_icon = tk.PhotoImage(file=self.path_gebpy + str("/GebPy_Logo.png"))
            self.parent.iconphoto(False, gebpy_icon)
        except:
            print("Image not found!")
        #
        ### Menu Bar
        menubar = tk.Menu(self.parent)
        self.parent.config(menu=menubar)
        ## Project
        project_menu = tk.Menu(menubar, tearoff=0)
        #
        project_menu.add_command(
            label="New", command=self.restart_gebpy)
        project_menu.add_separator()
        project_menu.add_command(
            label="Exit", command=self.parent.destroy)
        #
        menubar.add_cascade(
            label="Project",
            menu=project_menu)
        #
        ## Mineralogy
        mineralogy_menu = tk.Menu(menubar, tearoff=0)
        #
        sub_mineral_groups = tk.Menu(mineralogy_menu, tearoff=0)
        mineral_groups = [
            "Oxides", "Sulfides", "Sulfates", "Halides", "Tectosilicates", "Nesosilicates", "Sorosilicates",
            "Cyclosilicates", "Inosilicates", "Phyllosilicates", "Carbonates", "Phosphates", "Phospides",
            "Miscellaneous"]
        mineral_groups.sort()
        for mineral_group in mineral_groups:
            if mineral_group == "Oxides":
                sub_oxides = tk.Menu(sub_mineral_groups, tearoff=0)
                self.oxide_minerals = [
                    "Quartz", "Magnetite", "Hematite", "Al-Spinel", "Ilmenite", "Cassiterite", "Chromite", "Corundum",
                    "Rutile", "Pyrolusite", "Magnesiochromite", "Zincochromite", "Cr-Spinel", "Cuprospinel",
                    "Jacobsite", "Magnesioferrite", "Trevorite", "Franklinite", "Ulvospinel", "Fe-Spinel", "Uraninite",
                    "Litharge", "Massicot", "Minium", "Plattnerite", "Scrutinyite", "Zincite", "Columbite", "Tantalite",
                    "Coltan", "Crocoite", "Wulfenite", "Goethite", "Wolframite", "Huebnerite", "Ferberite", "Boehmite",
                    "Gibbsite", "Au(III)-Oxide", "Brucite", "Valentinite", "Senarmontite", "Ferberite-Huebnerite",
                    "Spinel", "Diaspore", "Cuprite", "Brookite", "Anatase", "Manganite", "Groutite", "Pyrophanite",
                    "Geikielite", "Claudetite", "Arsenolite", "Bismite", "Sphaerobismite", "Manganochromite",
                    "Nichromite", "Cochromite"]
                self.oxide_minerals.sort()
                for mineral in self.oxide_minerals:
                    sub_oxides.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Oxides",
                    menu=sub_oxides)
            elif mineral_group == "Carbonates":
                sub_carbonates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.carbonate_minerals = [
                    "Calcite", "Dolomite", "Magnesite", "Siderite", "Rhodochrosite", "Aragonite", "Cerrusite",
                    "Ankerite", "Azurite", "Malachite", "Ikaite", "Smithsonite"]
                self.carbonate_minerals.sort()
                for mineral in self.carbonate_minerals:
                    sub_carbonates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Carbonates",
                    menu=sub_carbonates)
            elif mineral_group == "Sulfides":
                sub_sulfides = tk.Menu(sub_mineral_groups, tearoff=0)
                self.sulfide_minerals = [
                    "Pyrite", "Chalcopyrite", "Galena", "Acanthite", "Chalcocite", "Bornite", "Sphalerite",
                    "Pyrrhotite", "Millerite", "Pentlandite", "Covellite", "Cinnabar", "Realgar", "Orpiment",
                    "Stibnite", "Marcasite", "Molybdenite", "Fahlore", "Gallite", "Roquesite", "Lenaite", "Laforetite",
                    "Vaesite", "Cattierite", "Cobaltite", "Marmatite"]
                self.sulfide_minerals.sort()
                for mineral in self.sulfide_minerals:
                    sub_sulfides.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Sulfides",
                    menu=sub_sulfides)
            elif mineral_group == "Sulfates":
                sub_sulfates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.sulfate_minerals = [
                    "Barite", "Celestite", "Anglesite", "Anhydrite", "Hanksite", "Gypsum", "Alunite", "Jarosite",
                    "Chalcanthite", "Kieserite", "Scheelite", "Hexahydrite", "Kainite"]
                self.sulfate_minerals.sort()
                for mineral in self.sulfate_minerals:
                    sub_sulfates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Sulfates",
                    menu=sub_sulfates)
            elif mineral_group == "Halides":
                sub_halides = tk.Menu(sub_mineral_groups, tearoff=0)
                self.halide_minerals = ["Halite", "Fluorite", "Sylvite", "Carnallite"]
                self.halide_minerals.sort()
                for mineral in self.halide_minerals:
                    sub_halides.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Halides",
                    menu=sub_halides)
            elif mineral_group == "Phosphates":
                sub_phosphates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.phosphate_minerals = ["Apatite", "Fluoroapatite", "Hydroxyapatite", "Chloroapatite"]
                self.phosphate_minerals.sort()
                for mineral in self.phosphate_minerals:
                    sub_phosphates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Phosphates",
                    menu=sub_phosphates)
            elif mineral_group == "Phospides":
                sub_phosphides = tk.Menu(sub_mineral_groups, tearoff=0)
                self.phospide_minerals = ["Allabogdanite"]
                self.phospide_minerals.sort()
                for mineral in self.phospide_minerals:
                    sub_phosphides.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Phospides",
                    menu=sub_phosphides)
            elif mineral_group == "Tectosilicates":
                sub_tectosilicates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.tectosilicate_minerals = [
                    "Alkaline Feldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline", "Orthoclase"]
                self.tectosilicate_minerals.sort()
                for mineral in self.tectosilicate_minerals:
                    sub_tectosilicates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Tectosilicates",
                    menu=sub_tectosilicates)
            elif mineral_group == "Sorosilicates":
                sub_sorosilicates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.sorosilicate_minerals = ["Epidote", "Zoisite", "Gehlenite"]
                self.sorosilicate_minerals.sort()
                for mineral in self.sorosilicate_minerals:
                    sub_sorosilicates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Sorosilicates",
                    menu=sub_sorosilicates)
            elif mineral_group == "Inosilicates":
                sub_inosilicates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.inosilicate_minerals = [
                    "Enstatite", "Ferrosilite", "Diopside", "Jadeite", "Aegirine", "Spodumene", "Wollastonite",
                    "Tremolite", "Actinolite", "Glaucophane", "Augite", "Arfvedsonite", "Ca-Amphibole", "Na-Amphibole",
                    "Mg-Fe-Pyroxene", "Ca-Pyroxene", "Na-Pyroxene", "Donpeacorite", "Orthopyroxene", "Clinopyroxene"]
                self.inosilicate_minerals.sort()
                for mineral in self.inosilicate_minerals:
                    sub_inosilicates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Inosilicates",
                    menu=sub_inosilicates)
            elif mineral_group == "Nesosilicates":
                sub_nesosilicates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.nesosilicate_minerals = [
                    "Zircon", "Titanite", "Thorite", "Andalusite", "Kyanite", "Sillimanite", "Topaz", "Staurolite",
                    "Fayalite", "Forsterite", "Tephroite", "Ca-Olivine", "Liebenbergite", "Olivine", "Pyrope",
                    "Almandine", "Grossular", "Anhadrite", "Uvarovite", "Al-Garnet", "Ca-Garnet"]
                self.nesosilicate_minerals.sort()
                for mineral in self.nesosilicate_minerals:
                    sub_nesosilicates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Nesosilicates",
                    menu=sub_nesosilicates)
            elif mineral_group == "Phyllosilicates":
                sub_phyllosilicates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.phyllosilicate_minerals = [
                    "Illite", "Kaolinite", "Montmorillonite", "Chamosite", "Clinochlore", "Pennantite", "Nimite",
                    "Chlorite", "Vermiculite", "Annite", "Phlogopite", "Eastonite", "Siderophyllite", "Biotite",
                    "Muscovite", "Glauconite", "Nontronite", "Saponite", "Talc", "Chrysotile", "Antigorite",
                    "Pyrophyllite"]
                self.phyllosilicate_minerals.sort()
                for mineral in self.phyllosilicate_minerals:
                    sub_phyllosilicates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Phyllosilicates",
                    menu=sub_phyllosilicates)
            elif mineral_group == "Cyclosilicates":
                sub_cyclosilicates = tk.Menu(sub_mineral_groups, tearoff=0)
                self.cyclosilicate_minerals = [
                    "Beryl", "Benitoite", "Cordierite", "Sekaninaite", "Schorl", "Elbaite", "Liddicoatite"]
                self.cyclosilicate_minerals.sort()
                for mineral in self.cyclosilicate_minerals:
                    sub_cyclosilicates.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Cyclosilicates",
                    menu=sub_cyclosilicates)
            elif mineral_group == "Miscellaneous":
                sub_miscellaneous = tk.Menu(sub_mineral_groups, tearoff=0)
                self.miscellaneous_minerals = ["Organic Matter"]
                self.miscellaneous_minerals.sort()
                for mineral in self.miscellaneous_minerals:
                    sub_miscellaneous.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Miscellaneous",
                    menu=sub_miscellaneous)
            else:
                sub_mineral_groups.add_command(
                    label=mineral_group)
        #
        sub_special_groups = tk.Menu(mineralogy_menu, tearoff=0)
        special_groups = [
            "Spinel Group", "Hematite Group", "Rutile Group", "Periclase Group", "Scheelite Group"]
        special_groups.sort()
        for special_group in special_groups:
            sub_special_groups.add_command(
                label=special_group, command=lambda var_name=special_group: self.select_special_group(var_name))
        #
        sub_ore_minerals = tk.Menu(mineralogy_menu, tearoff=0)
        ore_minerals = [
            "Fe Ores", "Pb Ores", "Zn Ores", "Cu Ores", "Al Ores", "Sn Ores"]
        ore_minerals.sort()
        for ore_mineral in ore_minerals:
            sub_ore_minerals.add_command(
                label=ore_mineral)
        #
        mineralogy_menu.add_cascade(
            label="Mineral Groups",
            menu=sub_mineral_groups)
        mineralogy_menu.add_cascade(
            label="Special Groups",
            menu=sub_special_groups)
        mineralogy_menu.add_cascade(
            label="Ore Minerals",
            menu=sub_ore_minerals)
        menubar.add_cascade(
            label="Mineralogy",
            menu=mineralogy_menu)
        #
        ## Petrology
        petrology_menu = tk.Menu(menubar, tearoff=0)
        #
        sub_rock_groups = tk.Menu(mineralogy_menu, tearoff=0)
        rock_groups = [
            "Sedimentary Rocks", "Igneous Rocks", "Metamorphic Rocks", "Evaporite Rocks", "Ore Rocks"]
        rock_groups.sort()
        for rock_group in rock_groups:
            if rock_group == "Sedimentary Rocks":
                sub_sedimentary = tk.Menu(petrology_menu, tearoff=0)
                sedimentary_rocks = {
                    "Siliciclastic": ["Sandstone", "Shale", "Mudstone", "Conglomerate", "Greywacke (Huckenholz)",
                                      "Siltstone"],
                    "Carbonate": ["Limestone", "Dolostone", "Marl"]}
                sedimentary_rocks = collections.OrderedDict(sorted(sedimentary_rocks.items()))
                i = 1
                n = len(sedimentary_rocks)
                for rock_group, rock_list in sedimentary_rocks.items():
                    rock_list.sort()
                    for rock in rock_list:
                        sub_sedimentary.add_command(
                            label=rock, command=lambda var_name=rock: self.select_rock(var_name))
                    if i < n:
                        sub_sedimentary.add_separator()
                        i += 1
                #
                sub_rock_groups.add_cascade(
                    label="Sedimentary Rocks",
                    menu=sub_sedimentary)
            elif rock_group == "Igneous Rocks":
                sub_igneous = tk.Menu(petrology_menu, tearoff=0)
                igneous_rocks = {
                    "Plutonic": ["Granite (Streckeisen)", "Granodiorite (Streckeisen)", "Tonalite (Streckeisen)",
                                 "Gabbro (Streckeisen)", "Diorite (Streckeisen)", "Monzonite (Streckeisen)",
                                 "Syenite (Streckeisen)", "Granitoid (Streckeisen)", "Quarzolite (Streckeisen)",
                                 "Foid-bearing Syenite (Streckeisen)", "Foid-bearing Monzonite (Streckeisen)",
                                 "Foid-bearing Monzodiorite (Streckeisen)", "Foid-bearing Monzogabbro (Streckeisen)",
                                 "Foid Monzosyenite (Streckeisen)", "Foid Monzodiorite (Streckeisen)",
                                 "Foid Monzogabbro (Streckeisen)", "Foidolite (Streckeisen)", "Norite (Streckeisen)",
                                 "Monzodiorite (Streckeisen)", "Monzogabbro (Streckeisen)"],
                    "Volcanic": ["Rhyolite (Streckeisen)", "Dacite (Streckeisen)", "Trachyte (Streckeisen)",
                                 "Latite (Streckeisen)", "Andesite (Streckeisen)", "Basalt (Streckeisen)",
                                 "Phonolite (Streckeisen)", "Tephrite (Streckeisen)", "Foidite (Streckeisen)",
                                 "Foid-bearing Trachyte (Streckeisen)", "Foid-bearing Latite (Streckeisen)",
                                 "Foid-bearing Andesite (Streckeisen)", "Foid-bearing Basalt (Streckeisen)",
                                 "Basalt (Yoder & Tilley)"],
                    "Ultramafic": ["Orthopyroxenite", "Clinopyroxenite", "Dunite", "Harzburgite", "Wehrlite",
                                   "Websterite", "Lherzolite", "Olivine-Websterite", "Olivine-Orthopyroxenite",
                                   "Olivine-Clinopyroxenite", "Peridotite", "Pyroxenite"]}
                #
                igneous_rocks = collections.OrderedDict(sorted(igneous_rocks.items()))
                i = 1
                n = len(igneous_rocks)
                for rock_group, rock_list in igneous_rocks.items():
                    rock_list.sort()
                    for rock in rock_list:
                        sub_igneous.add_command(
                            label=rock, command=lambda var_name=rock: self.select_rock(var_name))
                    if i < n:
                        sub_igneous.add_separator()
                        i += 1
                #
                sub_rock_groups.add_cascade(
                    label="Igneous Rocks",
                    menu=sub_igneous)
            #
            elif rock_group == "Ore Rocks":
                sub_ore = tk.Menu(petrology_menu, tearoff=0)
                ore_rocks = {
                    "Fe-ore": ["Itabirite", "Compact Hematite", "Friable Hematite", "Goethite Hematite",
                               "Al-rich Itabirite", "Compact Quartz Itabirite", "Friable Quartz Itabirite",
                               "Goethite Itabirite", "Banded Iron Formation"],
                    "Al-ore": ["Bauxite"]}
                #
                ore_rocks = collections.OrderedDict(sorted(ore_rocks.items()))
                i = 1
                n = len(ore_rocks)
                for rock_group, rock_list in ore_rocks.items():
                    rock_list.sort()
                    for rock in rock_list:
                        sub_ore.add_command(
                            label=rock, command=lambda var_name=rock: self.select_rock(var_name))
                    if i < n:
                        sub_ore.add_separator()
                        i += 1
                #
                sub_rock_groups.add_cascade(
                    label="Ore Rocks",
                    menu=sub_ore)
            #
            elif rock_group == "Metamorphic Rocks":
                sub_metamorphic = tk.Menu(petrology_menu, tearoff=0)
                metamorphic_rocks = {
                    "Granulite-Facies": ["Felsic Granulite", "Mafic Granulite"],
                    "Greenschist-Facies": ["Basaltic Greenschist", "Ultramafic Greenschist", "Pelitic Greenschist"],
                    "Amphibolite-Facies": ["Ortho-Amphibolite"]}
                #
                metamorphic_rocks = collections.OrderedDict(sorted(metamorphic_rocks.items()))
                i = 1
                n = len(metamorphic_rocks)
                for rock_group, rock_list in metamorphic_rocks.items():
                    rock_list.sort()
                    for rock in rock_list:
                        sub_metamorphic.add_command(
                            label=rock, command=lambda var_name=rock: self.select_rock(var_name))
                    if i < n:
                        sub_metamorphic.add_separator()
                        i += 1
                #
                sub_rock_groups.add_cascade(
                    label="Metamorphic Rocks",
                    menu=sub_metamorphic)
            else:
                sub_rock_groups.add_command(
                    label=rock_group)
        #
        petrology_menu.add_cascade(
            label="Rock Groups",
            menu=sub_rock_groups)
        #
        petrology_menu.add_separator()
        petrology_menu.add_command(
            label="Rock Analysis", command=self.rock_analysis)
        petrology_menu.add_command(
            label="Rock Comparison", command=self.rock_comparison)
        petrology_menu.add_command(
            label="Rock Builder", command=self.rock_builder)
        petrology_menu.add_command(
            label="Rock Builder (sedimentary)", command=self.rock_builder_sedimentary)
        #
        menubar.add_cascade(
            label="Petrology",
            menu=petrology_menu)
        #
        ## STRATIGRAPHY
        stratigraphy_menu = tk.Menu(menubar, tearoff=0)
        #
        sub_permian = tk.Menu(stratigraphy_menu, tearoff=0)
        permian_units = ["Zechstein", "Rotliegendes"]
        for unit in permian_units:
            sub_permian.add_command(label=unit, command=lambda var_unit=unit: self.real_sequences(var_unit))
        #
        sub_triassic = tk.Menu(stratigraphy_menu, tearoff=0)
        triassic_units = ["Keuper", "Muschelkalk"]
        for unit in triassic_units:
            sub_triassic.add_command(label=unit, command=lambda var_unit=unit: self.real_sequences(var_unit))
        #
        sub_buntstandstein = tk.Menu(stratigraphy_menu, tearoff=0)
        buntsandstein_units = ["Buntsandstein", "Upper Buntsandstein", "Medium Buntsandstein", "Lower Buntsandstein"]
        for unit in buntsandstein_units:
            sub_buntstandstein.add_command(label=unit, command=lambda var_unit=unit: self.real_sequences(var_unit))
            if unit == "Buntsandstein":
                sub_buntstandstein.add_separator()
        #
        # Real Sequences
        stratigraphy_menu.add_cascade(
            label="Triassic Units",
            menu=sub_triassic)
        stratigraphy_menu.add_cascade(
            label="Permian Units",
            menu=sub_permian)
        #
        sub_triassic.add_cascade(
            label="Buntsandstein",
            menu=sub_buntstandstein)
        #
        stratigraphy_menu.add_separator()
        stratigraphy_menu.add_command(
            label="Create Sequences")
        #
        menubar.add_cascade(
            label="Stratigraphy",
            menu=stratigraphy_menu)

        ## EXPLORATION
        exploration_menu = tk.Menu(menubar, tearoff=0)
        exploration_menu.add_command(
            label="Borehole data", command=lambda mode="Borehole data": self.select_exploration(mode))

        menubar.add_cascade(
            label="Exploration",
            menu=exploration_menu)

        ## DATABASE
        database_menu = tk.Menu(menubar, tearoff=0)
        database_menu.add_command(
            label="Elements")
        database_menu.add_command(
            label="Minerals")
        #
        menubar.add_cascade(
            label="Database",
            menu=database_menu)
        #
        ## Help
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(
            label="About")

        menubar.add_cascade(
            label="Help",
            menu=help_menu)

        ## Labels
        lbl_title = SimpleElements(
            parent=self.parent, row_id=n_rows - 2, column_id=1, n_rows=1, n_columns=30,
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Option"]).create_label(
            text=self.val_version, font_option="sans 11 bold", relief=tk.FLAT)

        ## Buttons
        btn_quit = SimpleElements(
            parent=self.parent, row_id=n_rows - 6, column_id=16, n_rows=3, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Quit GebPy", command=self.parent.quit)
        btn_restart = SimpleElements(
            parent=self.parent, row_id=n_rows - 6, column_id=1, n_rows=3, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Restart GebPy", command=self.restart_gebpy)

    ####################################################################################################################
    ## EXPLORATION #####################################################################################################
    ####################################################################################################################
    def select_exploration(self, mode):
        ExplorationInterface(parent=self.parent).create_subwindow_borehole_data()
        if mode == "Borehole data":
            ## Labels
            lbl_title = SimpleElements(
                parent=self.parent, row_id=5, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
                fg=self.colors_gebpy["Navigation"]).create_label(
                text="Borehole data", font_option="sans 14 bold", relief=tk.FLAT)
            lbl_datapoints = SimpleElements(
                parent=self.parent, row_id=7, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Number of data points", font_option="sans 10 bold", relief=tk.FLAT)
            lbl_units = SimpleElements(
                parent=self.parent, row_id=9, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Number of units", font_option="sans 10 bold", relief=tk.FLAT)



    #########################
    ## M i n e r a l o g y ##
    #########################
    #
    def select_mineral(self, name):
        ## Cleaning
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Static"][category]:
                gui_element.grid_remove()
        #
        for key_01 in self.gui_elements_sub["Mineral Chemistry"].keys():
            for key_02 in self.gui_elements_sub["Mineral Chemistry"][key_01].keys():
                for gui_item in self.gui_elements_sub["Mineral Chemistry"][key_01][key_02]:
                    gui_item.grid_remove()
                #
                self.gui_elements_sub["Mineral Chemistry"][key_01][key_02].clear()
        #
        ## Initialization
        self.gui_variables["Entry"]["Number Samples"] = tk.IntVar()
        self.gui_variables["Entry"]["Number Samples"].set(100)
        self.gui_variables["Radiobutton"]["Analysis Mode"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.gui_variables["Radiobutton"]["Trace Elements"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Trace Elements"].set(0)
        self.gui_variables["Radiobutton"]["Oxidation State"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Oxidation State"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Elements"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Elements"].set(0)
        self.gui_variables["Option Menu"]["Trace Elements"] = tk.StringVar()
        self.gui_variables["Option Menu"]["Trace Elements"].set("Select Trace Element")
        self.gui_variables["Checkbox"]["Trace Elements"] = {}
        self.gui_variables["Entry"]["Minimum"] = {}
        self.gui_variables["Entry"]["Maximum"] = {}
        self.gui_variables["Entry"]["Mean"] = {}
        self.gui_variables["Entry"]["Error"] = {}
        self.gui_variables["Entry"]["Trace Elements"] = {}
        self.gui_variables["Entry"]["Trace Elements"]["Minimum"] = {}
        self.gui_variables["Entry"]["Trace Elements"]["Maximum"] = {}
        self.trace_elements_all = {"All": []}
        categories_short = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U"]
        for category in categories_short:
            self.gui_variables["Entry"]["Minimum"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][category].set(0.0)
            self.gui_variables["Entry"]["Maximum"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][category].set(0.0)
            self.gui_variables["Entry"]["Mean"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][category].set(0.0)
            self.gui_variables["Entry"]["Error"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][category].set(0.0)
        self.btn_traces = None
        self.subwindow_traces = False
        #
        ## Labels
        lbl_mineralogy = SimpleElements(
            parent=self.parent, row_id=5, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Mineralogy", font_option="sans 14 bold", relief=tk.FLAT)
        lbl_name = SimpleElements(
            parent=self.parent, row_id=7, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text=name, font_option="sans 12 bold", relief=tk.FLAT)
        lbl_samples = SimpleElements(
            parent=self.parent, row_id=9, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Number of Datapoints", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_traces = SimpleElements(
            parent=self.parent, row_id=11, column_id=0, n_rows=6, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Trace Elements", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_addsetup = SimpleElements(
            parent=self.parent, row_id=31, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Additional Settings", font_option="sans 14 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].extend(
            [lbl_mineralogy, lbl_name, lbl_samples, lbl_traces, lbl_addsetup])
        #
        ## Entries
        entr_samples = SimpleElements(
            parent=self.parent, row_id=9, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Number Samples"])
        #
        self.gui_elements["Static"]["Entry"].extend([entr_samples])
        #
        ## Radiobuttons
        rb_trace_without = SimpleElements(
            parent=self.parent, row_id=11, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Without Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=0,
            color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
        self.gui_elements["Static"]["Radiobutton"].extend([rb_trace_without])
        if name in ["Quartz", "Orthoclase"]:
            rb_trace_with = SimpleElements(
                parent=self.parent, row_id=13, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="With Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=1,
                color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
            self.gui_elements["Static"]["Radiobutton"].extend([rb_trace_with])

        ## Buttons
        btn_generate_data = SimpleElements(
            parent=self.parent, row_id=18, column_id=16, n_rows=2, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Run Simulation", command=lambda var_name=name: self.run_simulation_mineralogy(var_name))
        #
        self.gui_elements["Static"]["Button"].append(btn_generate_data)
    #
    ####################
    ## SPECIAL GROUPS ##################################################################################################
    ####################
    #
    def select_special_group(self, var_name):
        ## Mineral Data
        #
        if var_name == "Spinel Group":
            mineral_names = np.sort(["Spinel", "Chihmingite", "Chromite", "Cochromite", "Coulsonite", "Cuprospinel",
                                     "Dellagiustaite", "Franklinite", "Gahnite", "Galaxite", "Hercynite", "Jacobsite",
                                     "Magnesiochromite", "Magnesiocoulsonite", "Magnesioferrite", "Magnetite",
                                     "Manganochromite", "Nichromite", "Thermaerogenite", "Trevorite", "Vuorelainenite",
                                     "Zincochromite"])
            mineral_group = {}
            #
            for mineral in mineral_names:
                try:
                    data_mineral = Oxides(mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=1)
                    mineral_group[mineral] = data_mineral
                except:
                    pass
            #
        elif var_name == "Hematite Group":
            mineral_names = np.sort(["Corundum", "Eskolaite", "Hematite", "Karelianite", "Tistarite"])
            mineral_group = {}
            #
            for mineral in mineral_names:
                try:
                    data_mineral = Oxides(mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=1)
                    mineral_group[mineral] = data_mineral
                except:
                    pass
            #
        elif var_name == "Rutile Group":
            mineral_names = np.sort(["Argutite", "Cassiterite", "Paratellurite", "Plattnerite", "Pyrolusite", "Rutile",
                                     "Stishovite"])
            mineral_group = {}
            #
            for mineral in mineral_names:
                try:
                    data_mineral = Oxides(mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=1)
                    mineral_group[mineral] = data_mineral
                except:
                    pass
            #
        elif var_name == "Periclase Group":
            mineral_names = np.sort(["Bunsenite", "Manganosite", "Periclase", "Wuestite"])
            mineral_group = {}
            #
            for mineral in mineral_names:
                try:
                    data_mineral = Oxides(mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=1)
                    mineral_group[mineral] = data_mineral
                except:
                    pass
            #
        elif var_name == "Scheelite Group":
            mineral_names = np.sort(["Powellite", "Scheelite", "Stolzite", "Wulfenite"])
            mineral_group = {}
            #
            for mineral in mineral_names:
                try:
                    data_mineral = Oxides(mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=1)
                    mineral_group[mineral] = data_mineral
                except:
                    pass
        #
        ## TABLE
        #
        start_column = 35
        list_categories = ["Mineral",
            "M (kg/mol)", "V (\u00C5\u00B3/mol)", "rho (kg/m\u00B3)", "vP (m/s)", "vS (m/s)",
            "K (GPa)", "G (GPa)", "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
        list_categories_short = ["Mineral", "M", "V", "rho", "vP", "vS", "K", "G", "E", "nu", "GR", "PE"]
        list_minerals = list(mineral_group.keys())
        list_width = list(80*np.ones(len(list_categories)))
        list_width = [int(item) for item in list_width]
        list_width[0] = 120
        #
        n_columns = 100
        tv_ma_results = SimpleElements(
            parent=self.parent, row_id=0, column_id=start_column, n_rows=20, n_columns=n_columns,
            fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
            n_categories=len(list_categories), text_n=list_categories,
            width_n=list_width, individual=True)
        #
        scb_v = ttk.Scrollbar(self.parent, orient="vertical")
        scb_h = ttk.Scrollbar(self.parent, orient="horizontal")
        tv_ma_results.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
        scb_v.config(command=tv_ma_results.yview)
        scb_h.config(command=tv_ma_results.xview)
        scb_v.grid(row=0, column=start_column + n_columns, rowspan=20, columnspan=1, sticky="ns")
        scb_h.grid(row=20, column=start_column, rowspan=1, columnspan=n_columns, sticky="ew")
        #
        for index, mineral in enumerate(list_minerals):
            entries = [mineral]
            #
            for category_parameter in list_categories_short[1:]:
                try:
                    entries.append(round(np.mean(mineral_group[mineral][category_parameter]), 3))
                except:
                    entries.append(mineral_group[mineral][category_parameter])
            #
            tv_ma_results.insert("", tk.END, values=entries)
        #
        ## DIAGRAMS
        #
        self.mineral_group_scatter(mineral_data=mineral_group)
        #
    def mineral_group_scatter(self, mineral_data):
        fig_mineral_group = Figure(
            figsize=(6, 6), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
        ax_mineral_group = fig_mineral_group.subplots(nrows=3, ncols=3)
        #
        categories = [["M", "V", "rho"], ["vP", "vS", "vP/vS"], ["GR", "PE", "nu"]]
        labels = [["M (kg/mol)", "V (A$^3$/mol", "rho (g/cm$^3$"], ["vP (m/s)", "vS (m/s)", "vP/vS (1)"],
                  ["GR (API)", "PE (barns/e$^-$)", "nu (1)"]]
        #
        helper_x = []
        helper_y = {}
        for i, subcategories in enumerate(categories):
            for j, subcategory in enumerate(subcategories):
                helper_y[subcategory] = []
                for mineral, dataset in mineral_data.items():
                    helper_x.append(dataset["rho"][-1]/1000)
                    helper_y[subcategory].append(dataset[subcategory][-1])
                    #
                    ax_mineral_group[i][j].scatter(dataset["rho"][-1]/1000, dataset[subcategory][-1])
                    #
        x_min = 0.9*min(helper_x)
        x_max = 1.1*max(helper_x)
        #
        for i, subcategories in enumerate(categories):
            for j, subcategory in enumerate(subcategories):
                y_min = 0.9*min(helper_y[subcategory])
                y_max = 1.1*max(helper_y[subcategory])
                #
                ax_mineral_group[i][j].set_xlim(left=x_min, right=x_max)
                ax_mineral_group[i][j].set_ylim(bottom=y_min, top=y_max)
                ax_mineral_group[i][j].set_xticks(np.around(
                    np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), 2))
                ax_mineral_group[i][j].set_yticks(np.around(
                    np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 2))
                ax_mineral_group[i][j].xaxis.set_tick_params(labelsize=8)
                ax_mineral_group[i][j].yaxis.set_tick_params(labelsize=8)
                ax_mineral_group[i][j].set_xlabel("Density - kg/m$^3$", fontsize=8)
                ax_mineral_group[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=8)
                ax_mineral_group[i][j].grid(True)
                ax_mineral_group[i][j].set_axisbelow(True)
        #
        canvas_mineral_group = FigureCanvasTkAgg(fig_mineral_group, master=self.parent)
        canvas_mineral_group.get_tk_widget().grid(
            row=21, column=35, rowspan=int(self.n_rows - 21), columnspan=int(2*(self.n_rows - 21)), sticky="nesw")

    def closest_number(self, n, m):
        # Find the quotient
        q = int(n/m)
        # 1st possible closest number
        n1 = m*q
        # 2nd possible closest number
        if ((n*m) > 0):
            n2 = (m*(q + 1))
        else:
            n2 = (m*(q - 1))
        # if true, then n1 is the required closest number
        if (abs(n - n1) < abs(n - n2)):
            return n1
        # else n2 is the required closest number
        return n2

    def myround(self, x, prec=2, base=.05):
        return round(base*round(float(x)/base), prec)

    def improve_xlimits(self, x_min, x_max):
        if x_min > 1:
            if (x_max - x_min) < 0.5:
                x_min = self.myround(x=x_min - 0.005, base=0.005)
                x_max = self.myround(x=x_max + 0.005, base=0.005)
            elif 0.5 <= (x_max - x_min) < 5:
                x_min = self.myround(x=x_min - 2, base=5)
                x_max = self.myround(x=x_max + 2, base=5)
            elif 5 <= (x_max - x_min) < 10:
                x_min = self.myround(x=x_min - 2, base=5)
                x_max = self.myround(x=x_max + 2, base=5)
            elif 10 <= (x_max - x_min) < 15:
                x_min = self.myround(x=x_min - 3, base=2.5)
                x_max = self.myround(x=x_max + 3, base=2.5)
            elif 15 <= (x_max - x_min) < 25:
                x_min = self.myround(x=x_min - 3.5, base=5)
                x_max = self.myround(x=x_max + 3.5, base=5)
            elif 25 <= (x_max - x_min) < 50:
                x_min = self.myround(x=x_min - 5, base=5)
                x_max = self.myround(x=x_max + 5, base=5)
            elif 50 <= (x_max - x_min) < 100:
                x_min = self.myround(x=x_min - 25, base=5)
                x_max = self.myround(x=x_max + 25, base=5)
            else:
                x_min = self.myround(x=x_min - 50, base=5)
                x_max = self.myround(x=x_max + 50, base=5)
        else:
            if 5 <= (x_max - x_min) < 10:
                x_min = self.myround(x=x_min - 5, base=5)
                x_max = self.myround(x=x_max + 5, base=5)
            elif 1 <= (x_max - x_min) < 5:
                x_min = self.myround(x=x_min - 0.5, base=0.5)
                x_max = self.myround(x=x_max + 0.5, base=0.5)
            elif 0.05 <= (x_max - x_min) < 0.5:
                x_min = self.myround(x=x_min - 0.025, base=0.025)
                x_max = self.myround(x=x_max + 0.025, base=0.025)
            elif 0.01 <= (x_max - x_min) < 0.05:
                x_min = self.myround(x=x_min - 0.005, base=0.005)
                x_max = self.myround(x=x_max + 0.005, base=0.005)
            else:
                delta_x = x_max - x_min
                x_min = self.myround(x=x_min - 0.0005, base=0.0005)
                x_max = self.myround(x=x_max + 0.0005, base=0.0005)

        return x_min, x_max

    def improve_ylimits(self, y_min, y_max):
        if (y_max - y_min) > 50:
            y_min = self.myround(x=y_min - 10, base=5)
            y_max = self.myround(x=y_max + 10, base=5)
        elif 50 >= (y_max - y_min) > 25:
            y_min = self.myround(x=y_min - 5, base=5)
            y_max = self.myround(x=y_max + 5, base=5)
        elif 25 >= (y_max - y_min) > 10:
            y_min = self.myround(x=y_min - 2.5, base=2.5)
            y_max = self.myround(x=y_max + 2.5, base=2.5)
        else:
            if y_min > 5:
                y_min = self.myround(x=y_min - 5, base=0.5)
                y_max = self.myround(x=y_max + 5, base=0.5)
            else:
                y_min = self.myround(x=y_min - 0.0005, base=0.005)
                y_max = self.myround(x=y_max + 0.0005, base=0.005)

        if y_min < 0:
            y_min = 0

        return y_min, y_max

    def change_rb_diagram(self):    # RB DIAGRAM MINERALOGY
        if self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 0 \
                and self.last_rb_analysis_mineral.get() not in [1, 2]:   # MINERAL PHYSICS
            ## Diagram
            if self.gui_variables["Radiobutton"]["Diagram Type Mineral"].get() == 0:    # HISTOGRAM
                if self.last_rb_diagram_mineral.get() != 0:
                    if "Mineralogy" in self.gui_elements["Temporary"]["Canvas"]:
                        categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS",
                                      "Mineral Chemistry Histogram", "Mineral Chemistry Scatter"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                    #
                    if "Mineral Physics Histogram" not in self.gui_elements["Temporary"]["Axis"]:
                        fig_mineralogy = Figure(
                            figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                        ax_mp_histo = fig_mineralogy.subplots(nrows=3, ncols=3)
                        #
                        categories = [["M", "V", "rho"], ["vP", "vS", "vP/vS"], ["GR", "PE", "nu"]]
                        labels = [["M (kg/mol)", "V (A$^3$/mol)", "rho (kg/m$^3$)"],
                                  ["vP (m/s)", "vS (m/s)", "vP/vS (1)"], ["GR (API)", "PE (barns/e$^-$)", "nu (1)"]]

                        for i, subcategories in enumerate(categories):
                            for j, key in enumerate(subcategories):
                                dataset_x = self.data_mineral[key]
                                y, x, _ = ax_mp_histo[i][j].hist(
                                    x=dataset_x, color=self.colors_gebpy["Option"], edgecolor="black",
                                    bins=12)
                                #
                                x_min = min(x)
                                x_max = max(x)
                                delta_x = round(x_max - x_min, 4)
                                y_min = 0
                                y_max = round(1.05*max(y), 2)
                                #
                                if delta_x < 1:
                                    n_digits = 3
                                elif 1 <= delta_x < 5:
                                    n_digits = 2
                                elif delta_x >= 5:
                                    n_digits = 0

                                x_min, x_max = self.improve_xlimits(x_min=x_min, x_max=x_max)
                                y_min, y_max = self.improve_ylimits(y_min=y_min, y_max=y_max)

                                if key != "nu":
                                    if x_min < 0:
                                        x_min = 0

                                if x_max - x_min == 0:
                                    pass
                                else:
                                    ax_mp_histo[i][j].set_xlim(left=x_min, right=x_max)
                                ax_mp_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                                if x_max - x_min == 0:
                                    pass
                                else:
                                    ax_mp_histo[i][j].set_xticks(np.around(
                                        np.linspace(x_min, x_max, 6, dtype=float, endpoint=True), n_digits))
                                ax_mp_histo[i][j].set_yticks(np.around(
                                    np.linspace(y_min, y_max, 6, dtype=float, endpoint=True), 1))
                                ax_mp_histo[i][j].xaxis.set_tick_params(labelsize="x-small")
                                ax_mp_histo[i][j].yaxis.set_tick_params(labelsize="x-small")
                                ax_mp_histo[i][j].set_xlabel(labels[i][j], fontsize="x-small")
                                ax_mp_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize="x-small")

                                ax_mp_histo[i][j].grid(which="major", axis="both", linestyle="-")
                                ax_mp_histo[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                ax_mp_histo[i][j].minorticks_on()
                                ax_mp_histo[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
                                ax_mp_histo[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
                                ax_mp_histo[i][j].set_axisbelow(True)

                            canvas_mineralogy = FigureCanvasTkAgg(fig_mineralogy, master=self.parent)
                            canvas_mineralogy.get_tk_widget().grid(
                                row=0, column=81, rowspan=int(self.n_rows - 3), columnspan=int(self.n_columns - 81),
                                sticky="nesw")
                            #
                            self.gui_elements["Temporary"]["Axis"]["Mineral Physics Histogram"] = ax_mp_histo
                            self.gui_elements["Temporary"]["Figure"]["Mineralogy"] = fig_mineralogy
                            self.gui_elements["Temporary"]["Canvas"]["Mineralogy"] = canvas_mineralogy
                        #
                    else:
                        ## Cleaning
                        for gui_axes in self.gui_elements["Temporary"]["Axis"]["Mineral Physics Histogram"]:
                            for gui_axis in gui_axes:
                                gui_axis.axis("on")
                                gui_axis.set_visible(True)
                        #
                        categories = ["Mineral Physics Scatter", "LA ICP MS", "Mineral Chemistry Histogram",
                                      "Mineral Chemistry Scatter"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                    #
                else:
                    pass
                #
            elif self.gui_variables["Radiobutton"]["Diagram Type Mineral"].get() == 1:  # SCATTER
                if self.last_rb_diagram_mineral.get() != 1:
                    if "Mineralogy" in self.gui_elements["Temporary"]["Canvas"]:
                        categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS",
                                      "Mineral Chemistry Histogram", "Mineral Chemistry Scatter"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                    #
                    if "Mineral Physics Scatter" not in self.gui_elements["Temporary"]["Axis"]:
                        if "Mineralogy" not in self.gui_elements["Temporary"]["Figure"]:
                            fig_mineralogy = Figure(
                                figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                        else:
                            fig_mineralogy = self.gui_elements["Temporary"]["Figure"]["Mineralogy"]
                        #
                        ax_mp_scatter = fig_mineralogy.subplots(nrows=3, ncols=3)
                        #
                        categories = [["M", "V", "rho"], ["vP", "vS", "vP/vS"], ["GR", "PE", "nu"]]
                        labels = [["M (kg/mol)", "V (A$^3$/mol)", "rho (kg/m$^3$)"],
                                  ["vP (m/s)", "vS (m/s)", "vP/vS (1)"], ["GR (API)", "PE (barns/e$^-$)", "nu (1)"]]
                        #
                        dataset_x = self.data_mineral["rho"]
                        #
                        for i, subcategories in enumerate(categories):
                            for j, key in enumerate(subcategories):
                                dataset_y = self.data_mineral[key]
                                ax_mp_scatter[i][j].scatter(
                                    dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black",
                                    alpha=0.5)

                                x_min = min(dataset_x)
                                x_max = max(dataset_x)
                                delta_x = round(x_max - x_min, 4)
                                y_min = min(dataset_y)
                                y_max = max(dataset_y)
                                delta_y = round(y_max - y_min, 4)

                                if delta_x < 1:
                                    n_digits_x = 3
                                    x_factor = 0.9
                                elif 1 <= delta_x < 5:
                                    n_digits_x = 2
                                    x_factor = 0.5
                                elif delta_x >= 5:
                                    n_digits_x = 0
                                    x_factor = 0.1

                                if delta_y < 1:
                                    n_digits_y = 3
                                    y_factor = 0.9
                                    var_dtype = float
                                elif 1 <= delta_y < 5:
                                    n_digits_y = 2
                                    y_factor = 0.5
                                    var_dtype = float
                                elif delta_y >= 5:
                                    n_digits_y = 0
                                    y_factor = 0.1
                                    var_dtype = int

                                # x_min, x_max = self.improve_xlimits(x_min=x_min, x_max=x_max)
                                # y_min, y_max = self.improve_ylimits(y_min=y_min, y_max=y_max)
                                #
                                # if key != "nu":
                                #     if x_min < 0:
                                #         x_min = 0

                                x_min = round(x_min - x_factor*delta_x, n_digits_x)
                                x_max = round(x_max + x_factor*delta_x, n_digits_x)
                                y_min = round(y_min - y_factor*delta_y, n_digits_y)
                                y_max = round(y_max + y_factor*delta_y, n_digits_y)

                                if key in ["rho"]:
                                    x_min = round(round(min(dataset_x) - 25, -1) + 15, n_digits_x)
                                    x_max = round(round(max(dataset_x) + 25, -1) - 15, n_digits_x)
                                if key in ["vP", "vS"]:
                                    y_min = round(round(min(dataset_y) - 150, -2), n_digits_y)
                                    y_max = round(round(max(dataset_y) + 150, -2), n_digits_y)
                                elif key in ["K", "G", "E", "GR"]:
                                    y_min = round(round(min(dataset_y) - 15, -1) + 10, n_digits_y)
                                    y_max = round(round(max(dataset_y) + 15, -1) - 10, n_digits_y)
                                elif key in ["phi"]:
                                    y_min = round(round(min(dataset_y) - 7, -1) + 5, n_digits_y)
                                    y_max = round(round(max(dataset_y) + 7, -1) - 5, n_digits_y)
                                elif key in ["PE", "vPvS", "nu"]:
                                    y_min = round(round(min(dataset_y) - 0.005, 3), n_digits_y)
                                    y_max = round(round(max(dataset_y) + 0.005, 3), n_digits_y)

                                if key != "nu":
                                    if x_min < 0:
                                        x_min = 0
                                if key != "nu":
                                    if y_min < 0:
                                        y_min = 0

                                step_x = (x_max - x_min)/4
                                step_x = int(round(step_x, n_digits_x))

                                if step_x == 0:
                                    x_value = (x_min + x_max)/2
                                    x_min = 0.9*x_value
                                    x_max = 1.1*x_value
                                    x_ticks = np.linspace(x_min, x_max, 5)
                                else:
                                    x_ticks = np.arange(x_min, x_max + step_x, step_x)

                                step_y = (y_max - y_min)/4

                                if var_dtype == float:
                                    step_y = round(step_y, n_digits_y)
                                else:
                                    step_y = int(round(step_y, n_digits_y))

                                if step_y > 0:
                                    y_ticks = np.arange(y_min, y_max + step_y, step_y)
                                else:
                                    y_value = (y_min + y_max)/2
                                    y_min = 0.9*y_value
                                    y_max = 1.1*y_value
                                    y_ticks = np.linspace(y_min, y_max, 5)

                                ax_mp_scatter[i][j].set_xlim(left=x_min, right=x_max)
                                ax_mp_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                                ax_mp_scatter[i][j].set_xticks(x_ticks)
                                ax_mp_scatter[i][j].set_yticks(y_ticks)
                                ax_mp_scatter[i][j].xaxis.set_tick_params(labelsize="x-small")
                                ax_mp_scatter[i][j].yaxis.set_tick_params(labelsize="x-small")
                                ax_mp_scatter[i][j].set_xlabel("Density - kg/m$^3$", fontsize="x-small")
                                ax_mp_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize="x-small")

                                ax_mp_scatter[i][j].grid(which="major", axis="both", linestyle="-")
                                ax_mp_scatter[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                ax_mp_scatter[i][j].minorticks_on()
                                ax_mp_scatter[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
                                ax_mp_scatter[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
                                ax_mp_scatter[i][j].set_axisbelow(True)

                        if "Mineralogy" not in self.gui_elements["Temporary"]["Canvas"]:
                            canvas_mineralogy = FigureCanvasTkAgg(fig_mineralogy, master=self.parent)
                            #
                        else:
                            canvas_mineralogy = self.gui_elements["Temporary"]["Canvas"]["Mineralogy"]
                        #
                        canvas_mineralogy.draw()
                        #
                        self.gui_elements["Temporary"]["Axis"]["Mineral Physics Scatter"] = ax_mp_scatter
                        self.gui_elements["Temporary"]["Canvas"]["Mineralogy"] = canvas_mineralogy
                        #
                    else:
                        ## Cleaning
                        for gui_axes in self.gui_elements["Temporary"]["Axis"]["Mineral Physics Scatter"]:
                            for gui_axis in gui_axes:
                                gui_axis.axis("on")
                                gui_axis.set_visible(True)
                        #
                        categories = ["Mineral Physics Histogram", "LA ICP MS", "Mineral Chemistry Histogram",
                                      "Mineral Chemistry Scatter"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                    #
                else:
                    pass
        elif self.last_rb_analysis_mineral.get() != 0 and self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 0:
            if self.gui_variables["Radiobutton"]["Diagram Type Mineral"].get() == 0:
                key = "Mineral Physics Histogram"
            else:
                key = "Mineral Physics Scatter"
            #
            categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS",
                          "Mineral Chemistry Histogram", "Mineral Chemistry Scatter"]
            for category in categories:
                if category in self.gui_elements["Temporary"]["Axis"]:
                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                        for gui_axis in gui_axes:
                            if category == key:
                                gui_axis.axis("on")
                                gui_axis.set_visible(True)
                            else:
                                gui_axis.axis("off")
                                gui_axis.set_visible(False)
            #
            self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
            #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 1: # MINERAL CHEMISTRY
            if self.gui_variables["Radiobutton"]["Diagram Type Elements"].get() == 0:   # HISTOGRAM
                if "Mineralogy" in self.gui_elements["Temporary"]["Canvas"]:
                    categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS",
                                  "Mineral Chemistry Histogram", "Mineral Chemistry Scatter"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                #
                if "Mineralogy" not in self.gui_elements["Temporary"]["Figure"]:
                    fig_mineralogy = Figure(
                        figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                else:
                    fig_mineralogy = self.gui_elements["Temporary"]["Figure"]["Mineralogy"]
                #
                ax_mc_histo = fig_mineralogy.subplots(nrows=3, ncols=3)
                #
                categories = []
                labels = []
                if self.gui_variables["Radiobutton"]["Concentration Type"].get() == 0:
                    var_key = "chemistry"
                else:
                    var_key = "compounds"
                #
                for index, (element, values) in enumerate(self.data_mineral[var_key].items()):
                    if index in [0, 3, 6]:
                        categories.append([])
                        labels.append([])
                    #
                    if 0 <= index < 3:
                        categories[0].append(element)
                        labels[0].append(str(element)+" (wt.%)")
                    elif 3 <= index < 6:
                        categories[1].append(element)
                        labels[1].append(str(element)+" (wt.%)")
                    elif 6 <= index < 9:
                        categories[2].append(element)
                        labels[2].append(str(element)+" (wt.%)")
                #
                for i, subcategories in enumerate(categories):
                    for j, key in enumerate(subcategories):
                        if self.gui_variables["Radiobutton"]["Concentration Type"].get() == 0:
                            dataset_x = np.array(self.data_mineral["chemistry"][key])*10**2
                        else:
                            dataset_x = np.array(self.data_mineral["compounds"][key])*10**2
                        #
                        y, x, _ = ax_mc_histo[i][j].hist(
                            x=dataset_x, color=self.colors_gebpy["Option"], edgecolor="black", bins=12)
                        #
                        x_min = min(x)
                        x_max = max(x)
                        delta_x = round(x_max - x_min, 4)
                        y_min = 0
                        y_max = round(1.05*max(y), 2)
                        #
                        if delta_x < 1:
                            n_digits = 2
                        elif 1 <= delta_x < 10:
                            n_digits = 1
                        elif delta_x >= 10:
                            n_digits = 0
                        #
                        x_min = round(x_min - 0.1*delta_x, n_digits)
                        x_max = round(x_max + 0.1*delta_x, n_digits)
                        #
                        if x_min < 0:
                            x_min = 0
                        #
                        ax_mc_histo[i][j].set_xlim(left=x_min, right=x_max)
                        ax_mc_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                        ax_mc_histo[i][j].set_xticks(np.around(
                            np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits))
                        ax_mc_histo[i][j].set_yticks(np.around(
                            np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 1))
                        ax_mc_histo[i][j].xaxis.set_tick_params(labelsize=8)
                        ax_mc_histo[i][j].yaxis.set_tick_params(labelsize=8)
                        ax_mc_histo[i][j].set_xlabel(labels[i][j], fontsize=8)
                        ax_mc_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=8)

                        ax_mc_histo[i][j].grid(which="major", axis="both", linestyle="-")
                        ax_mc_histo[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                        ax_mc_histo[i][j].minorticks_on()
                        ax_mc_histo[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(3))
                        ax_mc_histo[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(3))
                        ax_mc_histo[i][j].set_axisbelow(True)
                if "Mineralogy" not in self.gui_elements["Temporary"]["Canvas"]:
                    canvas_mineralogy = FigureCanvasTkAgg(fig_mineralogy, master=self.parent)
                    #
                else:
                    canvas_mineralogy = self.gui_elements["Temporary"]["Canvas"]["Mineralogy"]
                #
                canvas_mineralogy.draw()
                #
                self.gui_elements["Temporary"]["Axis"]["Mineral Chemistry Histogram"] = ax_mc_histo
                self.gui_elements["Temporary"]["Canvas"]["Mineralogy"] = canvas_mineralogy
                #
            elif self.gui_variables["Radiobutton"]["Diagram Type Elements"].get() == 1: # SCATTER
                if "Mineralogy" in self.gui_elements["Temporary"]["Canvas"]:
                    categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS",
                                  "Mineral Chemistry Histogram", "Mineral Chemistry Scatter"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                #
                if "Mineralogy" not in self.gui_elements["Temporary"]["Figure"]:
                    fig_mineralogy = Figure(
                        figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                else:
                    fig_mineralogy = self.gui_elements["Temporary"]["Figure"]["Mineralogy"]
                #
                ax_mc_scatter = fig_mineralogy.subplots(nrows=3, ncols=3)
                #
                categories = []
                labels = []
                ref_key = None
                ref_mean = 0
                if self.gui_variables["Radiobutton"]["Concentration Type"].get() == 0:
                    var_key = "chemistry"
                else:
                    var_key = "compounds"
                #
                for index, (element, values) in enumerate(self.data_mineral[var_key].items()):
                    mean_element = round(np.mean(values), 3)
                    if mean_element > ref_mean:
                        ref_key = element
                        ref_mean = mean_element
                    #
                    if index in [0, 3, 6]:
                        categories.append([])
                        labels.append([])
                    #
                    if 0 <= index < 3:
                        categories[0].append(element)
                        labels[0].append(str(element) + " (wt.%)")
                    elif 3 <= index < 6:
                        categories[1].append(element)
                        labels[1].append(str(element) + " (wt.%)")
                    elif 6 <= index < 9:
                        categories[2].append(element)
                        labels[2].append(str(element) + " (wt.%)")
                #
                dataset_x = np.array(self.data_mineral[var_key][ref_key])*10**2
                #
                for i, subcategories in enumerate(categories):
                    for j, key in enumerate(subcategories):
                        dataset_y = np.array(self.data_mineral[var_key][key])*10**2
                        ax_mc_scatter[i][j].scatter(
                            dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black", alpha=0.5)
                        #
                        x_min = min(dataset_x)
                        x_max = max(dataset_x)
                        delta_x = round(x_max - x_min, 4)
                        y_min = min(dataset_y)
                        y_max = max(dataset_y)
                        delta_y = round(y_max - y_min, 4)
                        #
                        if delta_x < 1:
                            n_digits_x = 2
                        elif 1 <= delta_x < 10:
                            n_digits_x = 1
                        elif delta_x >= 10:
                            n_digits_x = 0
                        if delta_y < 1:
                            n_digits_y = 2
                        elif 1 <= delta_y < 10:
                            n_digits_y = 1
                        elif delta_y >= 10:
                            n_digits_y = 0
                        #
                        x_min = round(x_min - 0.1*delta_x, n_digits_x)
                        x_max = round(x_max + 0.1*delta_x, n_digits_x)
                        y_min = round(y_min - 0.1*delta_y, n_digits_y)
                        y_max = round(y_max + 0.1*delta_y, n_digits_y)
                        #
                        if x_min < 0:
                            x_min = 0
                        if x_max > 100:
                            x_max = 100
                        if y_min < 0:
                            y_min = 0
                        if y_max > 100:
                            y_max = 100

                        if x_min == x_max:
                            x_value = (x_min + x_max)/2
                            x_min = 0.9*x_value
                            x_max = 1.1*x_value

                        if y_min == y_max:
                            y_value = (y_min + y_max)/2
                            y_min = 0.9*y_value
                            y_max = 1.1*y_value

                        ax_mc_scatter[i][j].set_xlim(left=x_min, right=x_max)
                        ax_mc_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                        ax_mc_scatter[i][j].set_xticks(
                            np.around(np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits_x))
                        ax_mc_scatter[i][j].set_yticks(
                            np.around(np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), n_digits_y))
                        ax_mc_scatter[i][j].xaxis.set_tick_params(labelsize=8)
                        ax_mc_scatter[i][j].yaxis.set_tick_params(labelsize=8)
                        ax_mc_scatter[i][j].set_xlabel(str(ref_key) + " (wt.%)", fontsize=8)
                        ax_mc_scatter[i][j].set_ylabel(labels[i][j], fontsize=8)

                        ax_mc_scatter[i][j].grid(which="major", axis="both", linestyle="-")
                        ax_mc_scatter[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                        ax_mc_scatter[i][j].minorticks_on()
                        ax_mc_scatter[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(3))
                        ax_mc_scatter[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(3))
                        ax_mc_scatter[i][j].set_axisbelow(True)
                #
                if "Mineralogy" not in self.gui_elements["Temporary"]["Canvas"]:
                    canvas_mineralogy = FigureCanvasTkAgg(fig_mineralogy, master=self.parent)
                    #
                else:
                    canvas_mineralogy = self.gui_elements["Temporary"]["Canvas"]["Mineralogy"]
                #
                canvas_mineralogy.draw()
                #
                self.gui_elements["Temporary"]["Axis"]["Mineral Chemistry Scatter"] = ax_mc_scatter
                self.gui_elements["Temporary"]["Canvas"]["Mineralogy"] = canvas_mineralogy
                #
        #
        self.last_rb_diagram_mineral.set(self.gui_variables["Radiobutton"]["Diagram Type Mineral"].get())
    #
    def change_rb_analysis(self):
        #
        start_row = 0
        start_column = 35
        #
        if self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 0:   # Mineral Physics
            if self.last_rb_analysis_mineral.get() != 0:
                ## Cleaning
                if "Mineralogy" in self.gui_elements["Temporary"]["Canvas"]:
                    categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                #
                categories = ["Label", "Radiobutton", "Entry", "Option Menu"]
                for category in categories:
                    for gui_item in self.gui_elements["Temporary"][category]:
                        gui_item.grid_remove()
                #
                ## Labels
                lbl_diagram_type = SimpleElements(
                    parent=self.parent, row_id=33, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
                #
                self.gui_elements["Temporary"]["Label"].extend([lbl_diagram_type])
                #
                ## Radiobuttons
                rb_diagram_type_01 = SimpleElements(
                    parent=self.parent, row_id=34, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Histogram", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Mineral"], value_rb=0,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram)
                rb_diagram_type_02 = SimpleElements(
                    parent=self.parent, row_id=36, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Scatter", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Mineral"], value_rb=1,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram)
                #
                self.gui_elements["Temporary"]["Radiobutton"].extend(
                    [rb_diagram_type_01, rb_diagram_type_02])
                #
                ## TREE VIEW
                categories = [
                    "M (kg/mol)", "V (\u00C5\u00B3/mol)", "rho (kg/m\u00B3)", "vP (m/s)", "vS (m/s)", "vP/vS (1)",
                    "K (GPa)", "G (GPa)", "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)", "U (barns/cm\u00B3)"]
                categories_short = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U"]
                list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
                list_width = list(75*np.ones(len(list_categories)))
                list_width = [int(item) for item in list_width]
                list_width[0] = 90
                list_width[-1] = 150

                self.tv_ma_results = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=30, n_columns=45,
                    fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
                    n_categories=len(list_categories), text_n=list_categories,
                    width_n=list_width, individual=True)

                scb_v = ttk.Scrollbar(self.parent, orient="vertical")
                scb_h = ttk.Scrollbar(self.parent, orient="horizontal")
                self.tv_ma_results.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
                scb_v.config(command=self.tv_ma_results.yview)
                scb_h.config(command=self.tv_ma_results.xview)
                scb_v.grid(row=0, column=start_column + 45, rowspan=30, columnspan=1, sticky="ns")
                scb_h.grid(row=30, column=start_column, rowspan=1, columnspan=45, sticky="ew")

                for index, category in enumerate(categories):
                    entries = [category]

                    n_digits = 3
                    var_entr_min = round(min(self.data_mineral[categories_short[index]]), n_digits)
                    var_entr_max = round(max(self.data_mineral[categories_short[index]]), n_digits)
                    var_entr_mean = round(np.mean(self.data_mineral[categories_short[index]]), n_digits)
                    var_entr_error = round(np.std(self.data_mineral[categories_short[index]], ddof=1), n_digits)

                    entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                    self.tv_ma_results.insert("", tk.END, values=entries)

                entries = ["-", "-", "-", "-", "-"]
                self.tv_ma_results.insert("", tk.END, values=entries)

                for element, dataset in self.data_mineral["chemistry"].items():
                    try:
                        if element in self.data_mineral["major elements"]:
                            entries = [str(element)+str(" (%)")]
                            var_factor = 100
                            n_digits = 2
                        else:
                            entries = [str(element) + str(" (ppm)")]
                            var_factor = 1000000
                            n_digits = 0
                    except:
                        entries = [str(element) + str(" (%)")]
                        var_factor = 100
                        n_digits = 2

                    #n_digits = 2
                    #var_factor = 100

                    var_entr_min = round(var_factor*min(dataset), n_digits)
                    var_entr_max = round(var_factor*max(dataset), n_digits)
                    var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
                    var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)

                    entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                    self.tv_ma_results.insert("", tk.END, values=entries)

                if "compounds" in self.data_mineral:
                    entries = ["-", "-", "-", "-", "-"]
                    self.tv_ma_results.insert("", tk.END, values=entries)
                    for compound, dataset in self.data_mineral["compounds"].items():
                        try:
                            if compound in self.data_mineral["major compounds"]:
                                entries = [str(compound) + str(" (%)")]
                                var_factor = 100
                            else:
                                entries = [str(compound) + str(" (ppm)")]
                                var_factor = 1000000
                        except:
                            entries = [str(compound) + str(" (%)")]
                            var_factor = 100

                        n_digits = 2
                        #var_factor = 100

                        var_entr_min = round(var_factor*min(dataset), n_digits)
                        var_entr_max = round(var_factor*max(dataset), n_digits)
                        var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
                        var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)

                        entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                        self.tv_ma_results.insert("", tk.END, values=entries)

                self.change_rb_diagram()

            else:
                pass
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 1:   # MINERAL CHEMISTRY
            if self.last_rb_analysis_mineral.get() != 1:
                ## Cleaning
                categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS"]
                for category in categories:
                    if category in self.gui_elements["Temporary"]["Axis"]:
                        for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                            for gui_axis in gui_axes:
                                gui_axis.axis("off")
                                gui_axis.set_visible(False)
                #
                self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                #
                categories = ["Label", "Radiobutton", "Entry", "Option Menu"]
                for category in categories:
                    for gui_item in self.gui_elements["Temporary"][category]:
                        gui_item.grid_remove()
                #
                ## Labels
                lbl_diagram_type = SimpleElements(
                    parent=self.parent, row_id=33, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
                lbl_concentration_setup = SimpleElements(
                    parent=self.parent, row_id=38, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Concentration Setup", font_option="sans 10 bold", relief=tk.FLAT)
                lbl_element = SimpleElements(
                    parent=self.parent, row_id=42, column_id=0, n_rows=2, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Element Selection", font_option="sans 10 bold", relief=tk.FLAT)
                #
                if self.oxides_present == True:
                    lbl_compound = SimpleElements(
                        parent=self.parent, row_id=44, column_id=0, n_rows=2, n_columns=14,
                        bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                        text="Oxide Selection", font_option="sans 10 bold", relief=tk.FLAT)
                    self.gui_variables["Option Menu"]["Amount Compound"].set("Select Oxide")
                    compound_list = ["Select Oxide"]
                    #
                if self.sulfides_present == True:
                    lbl_compound = SimpleElements(
                        parent=self.parent, row_id=44, column_id=0, n_rows=2, n_columns=14,
                        bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                        text="Sulfide Selection", font_option="sans 10 bold", relief=tk.FLAT)
                    self.gui_variables["Option Menu"]["Amount Compound"].set("Select Sulfide")
                    compound_list = ["Select Sulfide"]
                    #
                if self.phospides_present == True:
                    lbl_compound = SimpleElements(
                        parent=self.parent, row_id=44, column_id=0, n_rows=2, n_columns=14,
                        bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                        text="Phospide Selection", font_option="sans 10 bold", relief=tk.FLAT)
                    self.gui_variables["Option Menu"]["Amount Compound"].set("Select Phospide")
                    compound_list = ["Select Phospide"]
                #
                if "compounds" in self.data_mineral:
                    compound_list.extend(list(self.data_mineral["compounds"].keys()))
                #
                self.gui_elements["Temporary"]["Label"].extend(
                    [lbl_diagram_type, lbl_element, lbl_compound, lbl_concentration_setup, lbl_element])
                #
                ## Radiobuttons
                rb_diagram_type_01 = SimpleElements(
                    parent=self.parent, row_id=34, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Histogram", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Elements"], value_rb=0,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram)
                rb_diagram_type_02 = SimpleElements(
                    parent=self.parent, row_id=36, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Scatter", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Elements"], value_rb=1,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram)
                rb_concentration_type_01 = SimpleElements(
                    parent=self.parent, row_id=38, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Element Concentrations", var_rb=self.gui_variables["Radiobutton"]["Concentration Type"],
                    value_rb=0, color_bg=self.colors_gebpy["Background"],
                    command=lambda var_rb=self.gui_variables["Radiobutton"]["Concentration Type"]:
                    self.change_chemical_composition(var_rb))
                rb_concentration_type_02 = SimpleElements(
                    parent=self.parent, row_id=40, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Oxide Concentrations", var_rb=self.gui_variables["Radiobutton"]["Concentration Type"],
                    value_rb=1, color_bg=self.colors_gebpy["Background"],
                    command=lambda var_rb=self.gui_variables["Radiobutton"]["Concentration Type"]:
                    self.change_chemical_composition(var_rb))
                #
                self.gui_elements["Temporary"]["Radiobutton"].extend(
                    [rb_diagram_type_01, rb_diagram_type_02, rb_concentration_type_01, rb_concentration_type_02])
                #
                ## Option Menu
                self.list_elements.sort()
                opt_element = SimpleElements(
                    parent=self.parent, row_id=42, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_option_menu(
                    var_opt=self.gui_variables["Option Menu"]["Amount Element"],
                    var_opt_set=self.gui_variables["Option Menu"]["Amount Element"].get(), opt_list=self.list_elements,
                    active_bg=self.colors_gebpy["Accent"])
                opt_compound = SimpleElements(
                    parent=self.parent, row_id=44, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_option_menu(
                    var_opt=self.gui_variables["Option Menu"]["Amount Compound"],
                    var_opt_set=self.gui_variables["Option Menu"]["Amount Compound"].get(), opt_list=compound_list,
                    active_bg=self.colors_gebpy["Accent"])
                #
                self.gui_elements["Temporary"]["Option Menu"].extend([opt_element, opt_compound])
                #
                self.change_rb_diagram()
                #
            else:
                pass
        #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 2:   # SYNTHETIC LA-ICP-MS
            if self.last_rb_analysis_mineral.get() != 2:
                ## Cleaning
                categories = ["Mineral Physics Histogram", "Mineral Physics Scatter", "LA ICP MS",
                              "Mineral Chemistry Histogram", "Mineral Chemistry Scatter"]
                for category in categories:
                    if category in self.gui_elements["Temporary"]["Axis"]:
                        for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                            for gui_axis in gui_axes:
                                gui_axis.axis("off")
                                gui_axis.set_visible(False)
                #
                self.gui_elements["Temporary"]["Canvas"]["Mineralogy"].draw()
                #
                categories = ["Label", "Radiobutton", "Entry", "Option Menu"]
                for category in categories:
                    for gui_item in self.gui_elements["Temporary"][category]:
                        gui_item.grid_remove()
                #
                ## LA-ICP-MS Experiment Simulation
                time_data, intensity_data = self.simulate_laicpms_experiment()
                #
                ## Diagram
                if "LA ICP MS" not in self.gui_elements["Temporary"]["Canvas"]:
                    if "Mineralogy" not in self.gui_elements["Temporary"]["Figure"]:
                        fig_mineralogy = Figure(
                            dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                    else:
                        fig_mineralogy = self.gui_elements["Temporary"]["Figure"]["Mineralogy"]
                    #
                    ax_laicpms = [[fig_mineralogy.subplots(nrows=1, ncols=1)]]
                    #
                    element_laicpms = list(intensity_data.keys())
                    for element in element_laicpms:
                        ax_laicpms[0][0].plot(time_data, intensity_data[element], label=element, linewidth=1)
                        #
                    ax_laicpms[0][0].set_xlabel("Time (s)", fontsize=9)
                    ax_laicpms[0][0].set_ylabel("Intensity (cps)", labelpad=0.5, fontsize=9)
                    ax_laicpms[0][0].set_xlim(0, 60)
                    ax_laicpms[0][0].set_ylim(1, 5*10**9)
                    ax_laicpms[0][0].set_xticks(np.arange(0, 60 + 5, 5))
                    ax_laicpms[0][0].set_yscale("log")
                    ax_laicpms[0][0].grid(True)
                    #
                    ax_laicpms[0][0].grid(which="major", axis="both", linestyle="-")
                    ax_laicpms[0][0].minorticks_on()
                    ax_laicpms[0][0].grid(which="minor", axis="both", linestyle="--", alpha=0.25)
                    #
                    ax_laicpms[0][0].set_axisbelow(True)
                    #
                    ax_laicpms[0][0].legend(fontsize="x-small", framealpha=1.0)
                    #
                    canvas_mineralogy = self.gui_elements["Temporary"]["Canvas"]["Mineralogy"]
                    #
                    canvas_mineralogy.draw()
                    canvas_mineralogy.get_tk_widget().grid(
                        row=0, column=81, rowspan=int(self.n_rows - 3), columnspan=int(self.n_columns - 81), sticky="nesw")
                    #
                    self.gui_elements["Temporary"]["Axis"]["LA ICP MS"] = ax_laicpms
                    self.gui_elements["Temporary"]["Canvas"]["Mineralogy"] = canvas_mineralogy
                else:
                    self.gui_elements["Temporary"]["Canvas"]["LA ICP MS"].get_tk_widget().grid()
                #
            else:
                pass
        #
        self.last_rb_analysis_mineral.set(self.gui_variables["Radiobutton"]["Analysis Mode"].get())
    #
    def simulate_srm(self, var_srm="NIST 610 (GeoRem)", var_is="Si"):
        data_srm = {}
        if var_srm == "NIST 610 (GeoRem)":
            concentration_data = {
                "Li": 485, "Be": 466, "B": 356, "C": 0.0, "Na": 102970, "Mg": 465, "Al": 10791, "Si": 327091, "P": 409,
                "S": 570, "Cl": 470, "K": 465, "Ca": 81833, "Sc": 452, "Ti": 460, "V": 442, "Cr": 405, "Mn": 433,
                "Fe": 457, "Co": 405, "Ni": 459, "Cu": 430, "Zn": 456, "Ga": 438, "Ge": 426, "As": 317, "Se": 112,
                "Br": 33, "Rb": 426, "Sr": 516, "Y": 458, "Zr": 437, "Nb": 485, "Mo": 410, "Rh": 1, "Pd": 1, "Ag": 239,
                "Cd": 259, "In": 441, "Sn": 427, "Sb": 405, "Te": 327, "I": 0.0, "Cs": 357, "Ba": 454, "La": 440, "Ce": 458,
                "Pr": 443, "Nd": 437, "Sm": 453, "Eu": 444, "Gd": 456, "Tb": 440, "Dy": 436, "Ho": 440, "Er": 456,
                "Tm": 423, "Yb": 455, "Lu": 440, "Hf": 421, "Ta": 482, "W": 447, "Re": 50, "Pt": 3, "Au": 23, "Tl": 61,
                "Pb": 426, "Bi": 358, "Th": 457, "U": 462}

            normalized_sensitivity_data = {
                "Li": 2718.7, "Be": 348, "B": 600.7, "C": 0.0, "Na": 4011.6, "Mg": 52.3, "Al": 2282.0, "Si": 59.7,
                "Cl": 0.4, "K": 4938.7, "Ca": 126.2, "Ti": 13.3, "V": 5949.7, "Cr": 4862.2, "Mn": 2028.0, "Fe": 28.1,
                "Co": 5592.8, "Ni": 1180.1, "Cu": 3006.2, "Zn": 517.9, "Br": 19.6, "Rb": 7764.4, "Sr": 9145.0,
                "Sn": 3439, "I": 1298.6, "Cs": 11484.1, "Ba": 1198.2, "Pb": 3823.9}

            analytical_sensitivity_data = {
                "Li": 32.96, "Be": 4.49, "B": 5.27, "Na": 69.66, "Al": 49.34, "Si": 1.00, "K": 2.25, "Ti": 5.95,
                "Mn": 106.83, "Fe": 4.38, "Sn": 44.35}

            intensity_ratio_data = {}

            s_is = normalized_sensitivity_data[var_is]
            concentration_is = concentration_data[var_is]
            intensity_is = s_is*concentration_is
            for key, value in normalized_sensitivity_data.items():
                s_i = value
                xi_i = round(s_i/s_is, 3)
                analytical_sensitivity_data[key] = xi_i

                concentration_i = concentration_data[key]
                intensity_i = s_i*concentration_i
                ratio_i = intensity_i/intensity_is
                intensity_ratio_data[key] = ratio_i

        data_srm["Concentration"] = concentration_data
        data_srm["Normalized Sensitivity"] = normalized_sensitivity_data
        data_srm["Analytical Sensitivity"] = analytical_sensitivity_data
        data_srm["Intensity Ratio"] = intensity_ratio_data

        return data_srm

    def change_chemical_composition(self, var_rb):
        if var_rb.get() == 0:   # Element Composition
            ## Cleaning
            for gui_category in ["Label", "Entry"]:
                for gui_item in self.gui_elements_sub["Mineral Chemistry"]["Compound Concentration"][gui_category]:
                    gui_item.grid_remove()
            #
            if len(self.gui_elements_sub["Mineral Chemistry"]["Element Concentration"]["Label"]) == 0:
                for index, element in enumerate(self.list_elements):
                    lbl_element = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=element, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Temporary"]["Label"].append(lbl_element)
                    self.gui_elements_sub["Mineral Chemistry"]["Element Concentration"]["Label"].append(lbl_element)
                    #
                    ## Entries
                    var_entr_min = int(min(self.data_mineral["chemistry"][element])*10**6)
                    var_entr_max = int(max(self.data_mineral["chemistry"][element])*10**6)
                    var_entr_mean = int(np.mean(self.data_mineral["chemistry"][element])*10**6)
                    var_entr_error = int(np.std(self.data_mineral["chemistry"][element], ddof=1)*10**6)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 9, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][element], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 18, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][element], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 27, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][element], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 36, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][element], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                    self.gui_elements_sub["Mineral Chemistry"]["Element Concentration"]["Entry"].extend(
                        [entr_min, entr_max, entr_mean, entr_error])
            else:
                ## Reconstruction
                for gui_category in ["Label", "Entry"]:
                    for gui_item in self.gui_elements_sub["Mineral Chemistry"]["Element Concentration"][gui_category]:
                        gui_item.grid()
            #
        elif var_rb.get() == 1: # Oxide/Sulfide/Phospide Composition
            ## Cleaning
            for gui_category in ["Label", "Entry"]:
                for gui_item in self.gui_elements_sub["Mineral Chemistry"]["Element Concentration"][gui_category]:
                    gui_item.grid_remove()
            #
            if len(self.gui_elements_sub["Mineral Chemistry"]["Compound Concentration"]["Label"]) == 0:
                for index, compound in enumerate(self.list_compounds):
                    lbl_element = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=compound, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Temporary"]["Label"].append(lbl_element)
                    self.gui_elements_sub["Mineral Chemistry"]["Compound Concentration"]["Label"].append(lbl_element)
                    #
                    ## Entries
                    var_entr_min = int(min(self.data_mineral["compounds"][compound])*10**6)
                    var_entr_max = int(max(self.data_mineral["compounds"][compound])*10**6)
                    var_entr_mean = int(np.mean(self.data_mineral["compounds"][compound])*10**6)
                    var_entr_error = int(np.std(self.data_mineral["compounds"][compound], ddof=1)*10**6)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 9, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][compound], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 18, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][compound], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 27, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][compound], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=35 + 36, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][compound], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                    self.gui_elements_sub["Mineral Chemistry"]["Compound Concentration"]["Entry"].extend(
                        [entr_min, entr_max, entr_mean, entr_error])
            else:
                ## Reconstruction
                for gui_category in ["Label", "Entry"]:
                    for gui_item in self.gui_elements_sub["Mineral Chemistry"]["Compound Concentration"][gui_category]:
                        gui_item.grid()
            #
    #
    def simulate_laicpms_experiment(self):
        total_ppm = rd.randint(10**8, 10**9)
        time_step = 0.1
        time_data = list(np.around(np.arange(0, 60 + time_step, time_step), 1))
        #
        index_lower_start = time_data.index(10.1)
        index_upper_start = time_data.index(11.0)
        n_values_start = index_upper_start - index_lower_start
        amount_start = np.around(np.linspace(0.5, 1.0, n_values_start + 1, endpoint=True), 4)
        #
        index_lower_sig = time_data.index(11.1)
        index_upper_sig = time_data.index(50.0)
        n_values_sig = index_upper_sig - index_lower_sig
        amount_sig = np.around(np.linspace(1.0, 0.33, n_values_sig + 1, endpoint=True), 4)
        #
        index_lower_end = time_data.index(50.1)
        index_upper_end = time_data.index(52.0)
        n_values_end = index_upper_end - index_lower_end
        amount_end = np.around(np.linspace(1.0, 0.05, n_values_end + 1, endpoint=True), 4)
        #
        index_lower_end2 = time_data.index(52.1)
        index_upper_end2 = time_data.index(54.0)
        n_values_end2 = index_upper_end2 - index_lower_end2
        amount_end2 = np.around(np.geomspace(0.05, 0.0001, n_values_end2 + 1, endpoint=True), 4)
        #
        data_srm = self.simulate_srm()
        #
        intensity_data = {}
        if self.data_mineral["mineral"] in ["Qz", "Kfs", "Pl", "Scp", "Dnb", "Nph"]:
            target_id = self.list_elements.index("Si")
            self.list_elements.insert(0, self.list_elements.pop(target_id))
            var_list_element = self.list_elements
        else:
            var_list_element = self.list_elements
        #
        for element in var_list_element:
            index_start = 0
            index_sig = 0
            index_end = 0
            index_end2 = 0
            if element != "O":
                intensity_data[element] = []

                mean_bg = np.random.randint(1, 100)
                error_bg = np.random.uniform(0.1, 0.5)*mean_bg
                if self.data_mineral["mineral"] == "Qz":
                    concentration_is = np.mean(self.data_mineral["chemistry"]["Si"])
                    concentration_i = np.mean(self.data_mineral["chemistry"][element])
                    if element == "Si":
                        mean_sig_is = np.mean(self.data_mineral["chemistry"][element])*10**6*self.data_mineral[
                            "LA-ICP-MS"][element]
                    else:
                        mean_sig = (concentration_i/concentration_is)*mean_sig_is
                elif self.data_mineral["mineral"] in ["Kfs", "Pl", "Scp", "Nph", "Dnb"]:
                    oxide_is = self.data_mineral["LA-ICP-MS"]["Oxides"]["Si"]
                    oxide_i = self.data_mineral["LA-ICP-MS"]["Oxides"][element]
                    #concentration_is = np.mean(self.data_mineral["compounds"][oxide_is])
                    concentration_is = np.mean(self.data_mineral["chemistry"]["Si"])
                    sensitivity_i = data_srm["Analytical Sensitivity"][element]
                    s_i = data_srm["Normalized Sensitivity"][element]
                    s_is = data_srm["Normalized Sensitivity"]["Si"]
                    concentration_std_is = data_srm["Concentration"]["Si"]
                    intensity_std_is = s_is*concentration_std_is

                    if oxide_i not in ["Cl", "F", "Br", "I", "At"]:
                        #concentration_i = np.mean(self.data_mineral["compounds"][oxide_i])
                        concentration_i = np.mean(self.data_mineral["chemistry"][element])
                    else:
                        concentration_i = np.mean(self.data_mineral["chemistry"][element])

                    if element == "Si":
                        mean_sig_is = np.mean(self.data_mineral["chemistry"][element])*10**6*self.data_mineral[
                            "LA-ICP-MS"][element]
                        mean_sig = mean_sig_is
                    elif element in ["H", "C"]:
                        value = int(np.random.normal(loc=mean_bg, scale=error_bg, size=1)[0])
                        if value < 0:
                            value = 0
                        mean_sig = value
                    else:
                        #mean_sig = (concentration_i/concentration_is)*mean_sig_is*sensitivity_i
                        #mean_sig = (concentration_i/concentration_is)*mean_sig_is
                        mean_sig = (mean_sig_is/intensity_std_is)*(concentration_std_is/concentration_is)*concentration_i*s_i
                else:
                    mean_sig = np.mean(self.data_mineral["chemistry"][element])*total_ppm

                # error_sig = np.random.uniform(0, 0.025)*mean_sig

                for time_value in time_data:
                    if time_value <= 10:
                        value = int(np.random.normal(loc=mean_bg, scale=error_bg, size=1)[0])
                        if value < 0:
                            value = 0
                        value_bg = value
                    elif 10 < time_value <= 11:
                        error_sig = np.random.uniform(0.025, 0.075)*mean_sig
                        value = int(amount_start[index_start]*(
                                np.random.normal(loc=mean_sig, scale=error_sig, size=1)[0] + value_bg))
                        index_start += 1
                    elif 11 < time_value <= 50:
                        error_sig = np.random.uniform(0.05, 0.15)*mean_sig
                        value = int(amount_sig[index_sig]*(
                                np.random.normal(loc=mean_sig, scale=error_sig, size=1)[0] + value_bg))
                        last_value_sig = value
                        index_sig += 1
                    elif 50 < time_value <= 51:
                        error_sig = np.random.uniform(0.025, 0.075)*last_value_sig
                        value = int(amount_end[index_end]*(np.random.normal(
                            loc=last_value_sig, scale=error_sig, size=1)[0] + value_bg))
                        index_end += 1
                    elif time_value > 51:
                        value = int(np.random.normal(loc=mean_bg, scale=error_bg, size=1)[0])
                        if value < 0:
                            value = 0

                    intensity_data[element].append(value)

        return time_data, intensity_data

    def change_rb_traces(self, var_name):
        if self.gui_variables["Radiobutton"]["Trace Elements"].get() == 0:
            ## Cleaning
            for gui_item in self.gui_elements["Temporary"]["Button"]:
                gui_item.grid_remove()
        elif self.gui_variables["Radiobutton"]["Trace Elements"].get() == 1:
            ## Reconstruction
            for gui_item in self.gui_elements["Temporary"]["Button"]:
                if gui_item == self.btn_traces:
                    gui_item.grid()
            if var_name in self.oxide_minerals:
                data_mineral = Oxides(mineral=var_name, data_type=True).get_data()
                self.trace_elements_all = data_mineral["trace elements"]
            elif var_name in self.sulfide_minerals:
                data_mineral = Sulfides(mineral=var_name, data_type=True).get_data()
                self.trace_elements_all = data_mineral["trace elements"]
            elif var_name in self.tectosilicate_minerals:
                data_mineral = Tectosilicates(mineral=var_name, data_type=True).get_data()
                self.trace_elements_all = data_mineral["trace elements"]

            if self.btn_traces == None:
                self.btn_traces = SimpleElements(
                    parent=self.parent, row_id=15, column_id=16, n_rows=2, n_columns=15,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                    text="Trace Elements",
                    command=lambda var_traces=self.trace_elements_all:
                    self.select_trace_elements(var_traces))

            self.gui_elements["Temporary"]["Button"].append(self.btn_traces)

    def select_trace_elements(self, var_traces):
        self.window_trace_elements = tk.Toplevel(self.parent)
        self.window_trace_elements.title("Trace Elements")
        self.window_trace_elements.geometry("750x450")
        self.window_trace_elements.resizable(False, False)
        self.window_trace_elements["bg"] = self.colors_gebpy["Background"]
        #
        ## Cleaning
        categories = ["Frame", "Label", "Radiobutton", "Entry", "Checkbox"]
        priorities = ["Static", "Temporary"]
        for priority in priorities:
            for category in categories:
                if len(self.gui_elements_sub["Trace Elements"][priority][category]) > 0:
                    self.gui_elements_sub["Trace Elements"][priority][category].clear()
        #
        ## Geometry and Layout
        window_width = 750
        window_heigth = 450
        row_min = 15
        n_rows = int(window_heigth / row_min)
        column_min = 15
        n_columns = int(window_width / column_min)
        #
        for x in range(n_columns):
            tk.Grid.columnconfigure(self.window_trace_elements, x, weight=1)
        for y in range(n_rows):
            tk.Grid.rowconfigure(self.window_trace_elements, y, weight=1)
        #
        # Rows
        for i in range(0, n_rows):
            self.window_trace_elements.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            self.window_trace_elements.grid_columnconfigure(i, minsize=column_min)
        #
        ## Frames
        frm_navigation = SimpleElements(
            parent=self.window_trace_elements, row_id=0, column_id=0, n_rows=40, n_columns=12,
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Navigation"]).create_frame()
        frm_accent = SimpleElements(
            parent=self.window_trace_elements, row_id=0, column_id=13, n_rows=40, n_columns=1,
            bg=self.colors_gebpy["Accent Blue"], fg=self.colors_gebpy["Accent Blue"]).create_frame()
        #
        self.gui_elements_sub["Trace Elements"]["Static"]["Frame"].extend([frm_navigation, frm_accent])
        #
        ## Labels
        lbl_title = SimpleElements(
            parent=self.window_trace_elements, row_id=0, column_id=0, n_rows=2, n_columns=12,
            bg=self.colors_gebpy["Accent Blue"], fg=self.colors_gebpy["Navigation"]).create_label(
            text="Oxidation State", font_option="sans 12 bold", relief=tk.FLAT)
        #
        self.gui_elements_sub["Trace Elements"]["Static"]["Label"].append(lbl_title)

        self.oxidation_states = list(var_traces.keys())
        self.oxidation_states.remove("All")
        self.oxidation_states.sort(reverse=True)
        for index, oxidation_state in enumerate(self.oxidation_states):
            dict_traces = var_traces
            rb_oxidation_state = SimpleElements(
                parent=self.window_trace_elements, row_id=2*index + 2, column_id=0, n_rows=2, n_columns=12,
                bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Mainly "+str(oxidation_state), var_rb=self.gui_variables["Radiobutton"]["Oxidation State"],
                value_rb=index, color_bg=self.colors_gebpy["Navigation"],
                command=lambda var_traces=dict_traces: self.change_oxidation_state(var_traces))
            #
            self.gui_elements_sub["Trace Elements"]["Static"]["Radiobutton"].append(rb_oxidation_state)

        self.change_oxidation_state(var_traces=var_traces, first_run=True)
    #
    def change_oxidation_state(self, var_traces, first_run=False):
        if first_run == False:
            ## Cleaning
            categories = ["Label", "Entry", "Checkbox"]
            for category in categories:
                if len(self.gui_elements_sub["Trace Elements"]["Temporary"][category]) > 0:
                    for gui_item in self.gui_elements_sub["Trace Elements"]["Temporary"][category]:
                        gui_item.grid_remove()
                    self.gui_elements_sub["Trace Elements"]["Temporary"][category].clear()

        oxidation_state = self.oxidation_states[self.gui_variables["Radiobutton"]["Oxidation State"].get()]

        for index, trace_element in enumerate(var_traces[oxidation_state]):
            if index == 0:
                ## Labels
                lbl_mineral = SimpleElements(
                    parent=self.window_trace_elements, row_id=0, column_id=15, n_rows=2, n_columns=7,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Trace Element", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.window_trace_elements, row_id=0, column_id=22, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Min", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.window_trace_elements, row_id=0, column_id=26, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Max", font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements_sub["Trace Elements"]["Temporary"]["Label"].extend(
                    [lbl_mineral, lbl_min, lbl_max])
                #
            elif index == 12:
                ## Labels
                lbl_mineral = SimpleElements(
                    parent=self.window_trace_elements, row_id=0, column_id=31, n_rows=2, n_columns=7,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Trace Element", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.window_trace_elements, row_id=0, column_id=38, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Min", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.window_trace_elements, row_id=0, column_id=42, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Max", font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements_sub["Trace Elements"]["Temporary"]["Label"].extend(
                    [lbl_mineral, lbl_min, lbl_max])
            #
            if index < 12:
                ## Labels
                lbl_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1), column_id=15, n_rows=2, n_columns=5,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=trace_element, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if trace_element not in self.gui_variables["Checkbox"]["Trace Elements"]:
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element].set(0)
                cb_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1), column_id=20, n_rows=2, n_columns=2,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Trace Elements"][trace_element])
                #
                ## Entries
                if trace_element not in self.gui_variables["Entry"]["Trace Elements"]["Minimum"]:
                    self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element] = tk.IntVar()
                    self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element].set(0)
                    self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element] = tk.IntVar()
                    self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element].set(0)
                #
                entr_min = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * index + 2, column_id=22, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element])
                entr_max = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * index + 2, column_id=26, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element])
            #
            elif 12 <= index < 24:
                ## Labels
                lbl_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1), column_id=31, n_rows=2, n_columns=5,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=trace_element, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if trace_element not in self.gui_variables["Checkbox"]["Trace Elements"]:
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element].set(0)
                cb_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1), column_id=36, n_rows=2, n_columns=2,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Trace Elements"][trace_element])
                #
                ## Entries
                if trace_element not in self.gui_variables["Entry"]["Trace Elements"]["Minimum"]:
                    self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element] = tk.IntVar()
                    self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element].set(0)
                    self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element] = tk.IntVar()
                    self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element].set(0)
                #
                entr_min = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * index + 2, column_id=38, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element])
                entr_max = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * index + 2, column_id=42, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element])
            #
            self.gui_elements_sub["Trace Elements"]["Temporary"]["Label"].append(lbl_trace)
            self.gui_elements_sub["Trace Elements"]["Temporary"]["Checkbox"].append(cb_trace)
            self.gui_elements_sub["Trace Elements"]["Temporary"]["Entry"].extend([entr_min, entr_max])
    #
    def run_simulation_mineralogy(self, var_name):
        ## Cleaning
        for key_01 in self.gui_elements_sub["Mineral Chemistry"].keys():
            for key_02 in self.gui_elements_sub["Mineral Chemistry"][key_01].keys():
                for gui_item in self.gui_elements_sub["Mineral Chemistry"][key_01][key_02]:
                    gui_item.grid_remove()
                #
                self.gui_elements_sub["Mineral Chemistry"][key_01][key_02].clear()
        #
        self.last_rb_analysis_mineral.set(42)
        self.last_rb_diagram_mineral.set(42)
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"].set(0)
        #
        self.gui_variables["Radiobutton"]["Concentration Type"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Concentration Type"].set(0)
        self.gui_variables["Option Menu"]["Amount Element"] = tk.StringVar()
        self.gui_variables["Option Menu"]["Amount Element"].set("Select Element")
        self.gui_variables["Option Menu"]["Amount Compound"] = tk.StringVar()
        self.gui_variables["Option Menu"]["Amount Compound"].set("Select Compound")
        #
        self.trace_elements = {}
        for trace_element in self.trace_elements_all["All"]:
            if trace_element in self.gui_variables["Checkbox"]["Trace Elements"]:
                if self.gui_variables["Checkbox"]["Trace Elements"][trace_element].get() == 1:
                    self.trace_elements[trace_element] = {
                        "Min": self.gui_variables["Entry"]["Trace Elements"]["Minimum"][trace_element].get(),
                        "Max": self.gui_variables["Entry"]["Trace Elements"]["Maximum"][trace_element].get()}
        #
        self.traces_list = []
        for key, var_cb in self.gui_variables["Checkbox"]["Trace Elements"].items():
            if var_cb.get() == 1:
                self.traces_list.append(key)
        #
        self.oxides_present = False
        self.sulfides_present = False
        self.phospides_present = False
        #
        time_start = datetime.now()
        if var_name in self.oxide_minerals:         # Oxides
            self.data_mineral = Oxides(
                mineral=var_name, data_type=True, traces_list=self.trace_elements).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.carbonate_minerals:   # Carbonates
            self.data_mineral = Carbonates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.sulfide_minerals:   # Sulfides
            self.data_mineral = Sulfides(
                mineral=var_name, data_type=True, traces_list=self.trace_elements).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.sulfides_present = True
        elif var_name in self.sulfate_minerals:   # Sulfates
            self.data_mineral = Sulfates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.halide_minerals:   # Halides
            self.data_mineral = Halides(
                mineral=var_name, dict=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.phospide_minerals:   # Phospides
            self.data_mineral = Phospides(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.phospides_present = True
        elif var_name in self.phosphate_minerals:   # Phosphates
            self.data_mineral = Phosphates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.tectosilicate_minerals:   # Tectosilicates
            self.data_mineral = Tectosilicates(
                mineral=var_name, data_type=True, traces_list=self.trace_elements).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.nesosilicate_minerals:   # Nesosilicates
            self.data_mineral = Nesosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.sorosilicate_minerals:   # Sorosilicates
            self.data_mineral = Sorosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.inosilicate_minerals:   # Inosilicates
            self.data_mineral = Inosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.phyllosilicate_minerals:   # Phyllosilicates
            self.data_mineral = Phyllosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.cyclosilicate_minerals:   # Cyclosilicates
            self.data_mineral = Cyclosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        elif var_name in self.miscellaneous_minerals:   # Miscellaneous
            self.data_mineral = Organics(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
            self.oxides_present = True
        time_end = datetime.now()
        time_delta = (time_end - time_start)*1000
        print(f"Taken time: {time_delta.total_seconds()} ms")
        #
        self.list_compounds = []
        for key, dataset in self.data_mineral.items():
            if key == "chemistry":
                self.list_elements = list(dataset.keys())
            #
            if key == "compounds":
                self.list_compounds = list(dataset.keys())
        #
        for element in self.list_elements:
            self.gui_variables["Entry"]["Minimum"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][element].set(0.0)
            self.gui_variables["Entry"]["Maximum"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][element].set(0.0)
            self.gui_variables["Entry"]["Mean"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][element].set(0.0)
            self.gui_variables["Entry"]["Error"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][element].set(0.0)
        #
        for compound in self.list_compounds:
            self.gui_variables["Entry"]["Minimum"][compound] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][compound].set(0.0)
            self.gui_variables["Entry"]["Maximum"][compound] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][compound].set(0.0)
            self.gui_variables["Entry"]["Mean"][compound] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][compound].set(0.0)
            self.gui_variables["Entry"]["Error"][compound] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][compound].set(0.0)
        #
        ## Labels
        lbl_analysis = SimpleElements(
            parent=self.parent, row_id=21, column_id=0, n_rows=6, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Analysis Mode", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].append(lbl_analysis)
        #
        ## Radiobutton
        rb_geophysics = SimpleElements(
            parent=self.parent, row_id=21, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Physics", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=0,
            color_bg=self.colors_gebpy["Background"], command=self.change_rb_analysis)
        rb_geochemistry = SimpleElements(
            parent=self.parent, row_id=23, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Chemistry", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis)
        # rb_laicpms = SimpleElements(
        #     parent=self.parent, row_id=25, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
        #     fg=self.colors_gebpy["Background"]).create_radiobutton(
        #     text="Synthetic LA-ICP-MS", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=2,
        #     color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis)
        #
        self.gui_elements["Static"]["Radiobutton"].extend([rb_geophysics, rb_geochemistry])
        #
        ## Button
        name_mineral = var_name
        btn_export = SimpleElements(
            parent=self.parent, row_id=28, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Export Data", command=lambda var_dataset=self.data_mineral, var_name=name_mineral:
            self.export_mineral_data(var_dataset, var_name))
        #
        self.gui_elements["Static"]["Button"].append(btn_export)
        #
        for key, gui_element in self.gui_elements["Temporary"].items():
            if key not in ["Canvas", "Button"]:
                if type(gui_element) == list:
                    for gui_item in gui_element:
                        gui_item.grid_remove()
            elif key == "Canvas":
                for key_2, gui_item in gui_element.items():
                    gui_item.get_tk_widget().grid_remove()
            gui_element.clear()
        #
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.change_rb_analysis()
    #
    #######################
    ## P e t r o l o g y ##
    #######################
    #
    def select_rock(self, var_name):
        ## Cleaning
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Static"][category]:
                gui_element.grid_remove()
        #
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Rockbuilder Static"][category]:
                gui_element.grid_remove()
        #
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Rockbuilder Temporary"][category]:
                gui_element.grid_remove()
        #
        for key, gui_items in self.gui_elements["Rockbuilder Temporary"].items():
            if len(gui_items) > 0:
                if key not in ["Canvas"]:
                    if type(gui_items) == list:
                        for gui_item in gui_items:
                            gui_item.grid_remove()
                        gui_items.clear()
                else:
                    for key_2, gui_item in gui_items.items():
                        gui_item.get_tk_widget().grid_remove()
        #
        ## Initialization
        if var_name in ["Limestone"]:
            phi_max = 40
        elif var_name in ["Sandstone", "Dolostone", "Marl", "Mudstone"]:
            phi_max = 30
        else:
            phi_max = 10
        #
        self.gui_variables["Entry"]["Number Datapoints"] = tk.IntVar()
        self.gui_variables["Entry"]["Number Datapoints"].set(100)
        self.gui_variables["Entry"]["Porosity Min"] = tk.IntVar()
        self.gui_variables["Entry"]["Porosity Min"].set(0)
        self.gui_variables["Entry"]["Porosity Max"] = tk.IntVar()
        self.gui_variables["Entry"]["Porosity Max"].set(phi_max)
        #
        ## Labels
        lbl_title = SimpleElements(
            parent=self.parent, row_id=5, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Petrology", font_option="sans 14 bold", relief=tk.FLAT)
        lbl_name = SimpleElements(
            parent=self.parent, row_id=7, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text=var_name, font_option="sans 12 bold", relief=tk.FLAT)
        lbl_samples = SimpleElements(
            parent=self.parent, row_id=9, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Number of Datapoints", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_phi_min = SimpleElements(
            parent=self.parent, row_id=12, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Porosity Minimum", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_phi_max = SimpleElements(
            parent=self.parent, row_id=14, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Porosity Maximum", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].extend([lbl_title, lbl_name, lbl_samples, lbl_phi_min, lbl_phi_max])
        #
        ## Entries
        entr_samples = SimpleElements(
            parent=self.parent, row_id=9, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Number Datapoints"])
        entr_phi_min = SimpleElements(
            parent=self.parent, row_id=12, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Porosity Min"])
        entr_phi_max = SimpleElements(
            parent=self.parent, row_id=14, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Porosity Max"])
        #
        self.gui_elements["Static"]["Entry"].extend([entr_samples, entr_phi_min, entr_phi_max])
        #
        ## Buttons
        name_rock = var_name
        btn_simulation = SimpleElements(
            parent=self.parent, row_id=17, column_id=16, n_rows=2, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Run Simulation", command=lambda var_name=name_rock: self.run_simulation_petrology(var_name))

        self.gui_elements["Static"]["Button"].extend([btn_simulation])

        ## TREEVIEWS
        list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
        list_width = list(75*np.ones(len(list_categories)))
        list_width = [int(item) for item in list_width]
        list_width[0] = 90
        list_width[-1] = 150
        #
        self.tv_petrology_results = SimpleElements(
            parent=self.parent, row_id=0, column_id=35, n_rows=30, n_columns=45,
            fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
            n_categories=len(list_categories), text_n=list_categories,
            width_n=list_width, individual=True)
        #
        scb_v = ttk.Scrollbar(self.parent, orient="vertical")
        scb_h = ttk.Scrollbar(self.parent, orient="horizontal")
        self.tv_petrology_results.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
        scb_v.config(command=self.tv_petrology_results.yview)
        scb_h.config(command=self.tv_petrology_results.xview)
        scb_v.grid(row=0, column=35 + 45, rowspan=30, columnspan=1, sticky="ns")
        scb_h.grid(row=30, column=35, rowspan=1, columnspan=45, sticky="ew")
    #
    def run_simulation_petrology(self, var_name):
        try:
            self.canvas_scatter.get_tk_widget().grid_forget()
        except:
            pass
        #
        ## Initialization
        self.gui_variables["Radiobutton"]["Analysis Mode"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Rock"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Rock"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Element"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Element"].set(0)
        self.gui_variables["Radiobutton"]["Concentration Type"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Concentration Type"].set(0)
        self.gui_variables["Radiobutton"]["Amount Mineral"] = tk.StringVar()
        self.gui_variables["Radiobutton"]["Amount Mineral"].set("Select Mineral")
        self.gui_variables["Radiobutton"]["Amount Element"] = tk.StringVar()
        self.gui_variables["Radiobutton"]["Amount Element"].set("Select Element")
        self.gui_variables["Radiobutton"]["Amount Oxide"] = tk.StringVar()
        self.gui_variables["Radiobutton"]["Amount Oxide"].set("Select Oxide")
        #
        self.selected_minerals = {}
        n_digits = 8
        for key, container in self.gui_elements["Temporary"].items():
            if len(container) > 0:
                container.clear()
        for key, item in self.last_rb_setting.items():
            for key2, item2 in item.items():
                item2.set(42)
        #
        self.gui_variables["Entry"]["Minimum"] = {}
        self.gui_variables["Entry"]["Maximum"] = {}
        self.gui_variables["Entry"]["Mean"] = {}
        self.gui_variables["Entry"]["Error"] = {}
        categories_short = ["rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi"]
        for category in categories_short:
            self.gui_variables["Entry"]["Minimum"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][category].set(0.0)
            self.gui_variables["Entry"]["Maximum"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][category].set(0.0)
            self.gui_variables["Entry"]["Mean"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][category].set(0.0)
            self.gui_variables["Entry"]["Error"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][category].set(0.0)
        #
        ## Labels
        lbl_analysis = SimpleElements(
            parent=self.parent, row_id=22, column_id=0, n_rows=6, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Analysis Mode", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].append(lbl_analysis)

        time_start = datetime.now()
        ## Siliciclastic Rocks
        if var_name == "Sandstone":
            data = SiliciclasticRocks(fluid="water", actualThickness=0).create_sandstone(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Conglomerate":
            data = SiliciclasticRocks(fluid="water", actualThickness=0).create_conglomerate(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Siltstone":
            data = SiliciclasticRocks(fluid="water", actualThickness=0).create_siltstone(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Mudstone":
            data = SiliciclasticRocks(fluid="water", actualThickness=0).create_mudstone_alt(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Shale":
            data = SiliciclasticRocks(fluid="water", actualThickness=0).create_shale_alt(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Greywacke (Huckenholz)":
            data = SiliciclasticRocks(fluid="water", actualThickness=0).create_greywacke_huckenholz(
                rock="Greywacke", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        #
        ## Carbonate Rocks
        elif var_name == "Limestone":
            data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Dolostone":
            data = CarbonateRocks(fluid="water", actualThickness=0).create_dolostone(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Marl":
            data = CarbonateRocks(fluid="water", actualThickness=0).create_marl(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        #
        ## Igneous Rocks (Plutonic)
        elif var_name == "Granite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Granite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Granodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Granodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Tonalite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Tonalite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Gabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Gabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Ca",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Norite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Norite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl=[0.1, 0.5],
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Diorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Diorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Na",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Monzodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Monzodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Na",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Monzogabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Monzogabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Ca",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Monzonite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Monzonite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Syenite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Syenite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Granitoid (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Granitoid", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Quarzolite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Quarzolite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid-bearing Syenite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Syenite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False,
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid-bearing Monzonite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Monzonite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False,
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid-bearing Monzodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Monzodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Na",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid-bearing Monzogabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Monzogabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Ca",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid Monzosyenite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid Monzosyenite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False,
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid Monzodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid Monzodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Na",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foid Monzogabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid Monzogabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Ca",
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Foidolite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foidolite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False,
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        #
        ## Igneous Rocks (Volcanic)
        elif var_name == "Rhyolite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Rhyolite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Dacite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Dacite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Trachyte (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Trachyte", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Latite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Latite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Andesite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Andesite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Na")
        elif var_name == "Basalt (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Basalt", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Ca")
        elif var_name == "Foid-bearing Trachyte (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Foid-bearing Trachyte", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
        elif var_name == "Foid-bearing Latite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Foid-bearing Latite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
        elif var_name == "Foid-bearing Andesite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Foid-bearing Andesite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                enrichment_pl="Na", upper_streckeisen=False)
        elif var_name == "Foid-bearing Basalt (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Foid-bearing Basalt", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                enrichment_pl="Ca", upper_streckeisen=False)
        elif var_name == "Phonolite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Phonolite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
        elif var_name == "Tephrite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Tephrite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), upper_streckeisen=False)
        elif var_name == "Foidite (Streckeisen)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_volcanic_rock_streckeisen(
                rock="Foidite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), upper_streckeisen=False)
        elif var_name == "Basalt (Yoder & Tilley)":
            data = Volcanic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_basaltic_rock_yoder_tilley(
                rock="Basalt", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        #
        ## Igneous Rocks (Ultramafic)
        elif var_name == "Orthopyroxenite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Orthopyroxenite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Clinpyroxenite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Clinpyroxenite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Dunite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Dunite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Harzburgite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Harzburgite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Wehrlite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Wehrlite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Websterite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Websterite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Lherzolite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Lherzolite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Olivine-Websterite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Olivine-Websterite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Olivine-Orthopyroxenite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Olivine-Orthopyroxenite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Olivine-Clinopyroxenite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Olivine-Clinopyroxenite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Peridotite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Peridotite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Pyroxenite":
            data = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_ultramafic_rock(
                rock="Pyroxenite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        #
        ## Ore Rocks
        # Fe-Ore
        elif var_name == "Itabirite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock="Itabirite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Compact Hematite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock="Compact Hematite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                classification=var_name)
        elif var_name == "Banded Iron Formation":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_bandedironformation(
                rock="Banded Iron Formation", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Friable Hematite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock="Friable Hematite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                classification=var_name)
        elif var_name == "Goethite Hematite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock=var_name, number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification=var_name)
        elif var_name == "Al-rich Itabirite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock=var_name, number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification=var_name)
        elif var_name == "Compact Quartz Itabirite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock=var_name, number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification=var_name)
        elif var_name == "Friable Quartz Itabirite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock=var_name, number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification=var_name)
        elif var_name == "Goethite Itabirite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_siliciclastic_itabirite(
                rock=var_name, number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification=var_name)
        # Al-ore
        elif var_name == "Bauxite":
            data = OreRocks(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_bauxite(
                rock="Bauxite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        #
        ## Metamorphic Rocks
        # Granulite-Facies
        elif var_name == "Felsic Granulite":
            data = GranuliteFacies(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_granulite(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification="felsic")
        elif var_name == "Mafic Granulite":
            data = GranuliteFacies(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_granulite(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(), classification="mafic")
        # Greenschist-Facies
        elif var_name == "Basaltic Greenschist":
            data = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_greenschist_basaltic_alt(
                number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Ultramafic Greenschist":
            data = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_greenschist_ultramafic_alt(
                number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Pelitic Greenschist":
            data = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_greenschist_pelitic_alt(
                number=self.gui_variables["Entry"]["Number Datapoints"].get())
        # Amphibolite-Facies
        elif var_name == "Ortho-Amphibolite":
            data = AmphiboliteFacies(
                fluid="water", actual_thickness=0, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_amphibolite_ortho(
                number=self.gui_variables["Entry"]["Number Datapoints"].get())

        time_end = datetime.now()
        time_delta = (time_end - time_start)*1000
        print(f"Taken time: {time_delta.total_seconds()} ms")

        self.data_rock = {}
        self.rock_data = {}

        categories = ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi", "phi_true", "fluid",
                      "mineralogy", "chemistry"]

        if "compounds" in data:
            categories.append("compounds")

        for category in categories:
            if category in ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]:
                self.data_rock[category] = data[category]
                self.rock_data[category] = data[category]
            elif category in ["mineralogy", "chemistry", "compounds"]:
                self.data_rock[category] = data[category]
                self.rock_data[category] = data[category]
            elif category in ["phi", "phi_true"]:
                if category in data:
                    try:
                        self.data_rock[category] = list(np.array(data[category])*100)
                        self.rock_data[category] = list(np.array(data[category])*100)
                    except:
                        self.data_rock[category] = [data[category]*100]
                        self.rock_data[category] = [data[category]*100]

        self.list_elements_rock = list(self.data_rock["chemistry"].keys())
        self.list_minerals_rock = list(self.data_rock["mineralogy"].keys())

        if "compounds" in data:
            self.list_compounds_rock = list(self.data_rock["compounds"].keys())

        ## Radiobuttons
        rb_geophysics = SimpleElements(
            parent=self.parent, row_id=22, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Physical Data", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=0,
            color_bg=self.colors_gebpy["Background"], command=self.change_rb_analysis_rocks)
        rb_geochemistry = SimpleElements(
            parent=self.parent, row_id=24, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Composition", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis_rocks)
        rb_geochemistry2 = SimpleElements(
            parent=self.parent, row_id=26, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Element Composition", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=2,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis_rocks)

        self.gui_elements["Static"]["Radiobutton"].extend([rb_geophysics, rb_geochemistry, rb_geochemistry2])

        ## Button
        name_rock = var_name
        btn_export = SimpleElements(
            parent=self.parent, row_id=29, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Export Data", command=lambda var_dataset=self.data_rock, var_name=name_rock:
            self.export_rock_data(var_dataset, var_name))

        self.gui_elements["Static"]["Button"].append(btn_export)

        for key, gui_element in self.gui_elements["Temporary"].items():
            if key not in ["Canvas", "Button"]:
                if type(gui_element) == list:
                    for gui_item in gui_element:
                        gui_item.grid_remove()
            elif key == "Canvas":
                for key_2, gui_item in gui_element.items():
                    gui_item.get_tk_widget().grid_remove()
            gui_element.clear()

        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.last_rb_analysis_rock.set(42)
        self.change_rb_analysis_rocks()

        for mineral in self.list_minerals_rock:
            self.gui_variables["Entry"]["Minimum"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][mineral].set(0.0)
            self.gui_variables["Entry"]["Maximum"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][mineral].set(0.0)
            self.gui_variables["Entry"]["Mean"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][mineral].set(0.0)
            self.gui_variables["Entry"]["Error"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][mineral].set(0.0)

        for element in self.list_elements_rock:
            self.gui_variables["Entry"]["Minimum"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][element].set(0.0)
            self.gui_variables["Entry"]["Maximum"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][element].set(0.0)
            self.gui_variables["Entry"]["Mean"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][element].set(0.0)
            self.gui_variables["Entry"]["Error"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][element].set(0.0)

        if len(self.tv_petrology_results.get_children()) > 0:
            for item in self.tv_petrology_results.get_children():
                self.tv_petrology_results.delete(item)

        for property, dataset in self.rock_data.items():
            entries = [property]

            if property not in ["mineralogy", "chemistry", "compounds"]:
                n_digits = 2

                var_entr_min = round(min(dataset), n_digits)
                var_entr_max = round(max(dataset), n_digits)
                var_entr_mean = round(np.mean(dataset), n_digits)
                var_entr_error = round(np.std(dataset, ddof=1), n_digits)

                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                self.tv_petrology_results.insert("", tk.END, values=entries)

        entries = ["-", "-", "-", "-", "-"]
        self.tv_petrology_results.insert("", tk.END, values=entries)

        for mineral in np.sort(self.list_minerals_rock):
            dataset = self.rock_data["mineralogy"][mineral]
            entries = [str(mineral)+str(" (%)")]

            n_digits = 2
            var_factor = 100

            var_entr_min = round(var_factor*min(dataset), n_digits)
            var_entr_max = round(var_factor*max(dataset), n_digits)
            var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
            var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)

            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

            self.tv_petrology_results.insert("", tk.END, values=entries)

        entries = ["-", "-", "-", "-", "-"]
        self.tv_petrology_results.insert("", tk.END, values=entries)

        for element in np.sort(self.list_elements_rock):
            dataset = self.rock_data["chemistry"][element]
            entries = [str(element)+str(" (%)")]
            n_digits = 2
            var_factor = 100
            #
            var_entr_min = round(var_factor*min(dataset), n_digits)
            var_entr_max = round(var_factor*max(dataset), n_digits)
            var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
            var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)
            #
            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            #
            self.tv_petrology_results.insert("", tk.END, values=entries)
        #
        entries = ["-", "-", "-", "-", "-"]
        self.tv_petrology_results.insert("", tk.END, values=entries)

        if "compounds" in data:
            for element in np.sort(self.list_compounds_rock):
                dataset = self.rock_data["compounds"][element]
                entries = [str(element) + str(" (%)")]
                n_digits = 2
                var_factor = 100

                var_entr_min = round(var_factor*min(dataset), n_digits)
                var_entr_max = round(var_factor*max(dataset), n_digits)
                var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
                var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)

                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                self.tv_petrology_results.insert("", tk.END, values=entries)

    def rock_analysis(self):
        rock_analysis = tk.Toplevel(self.parent)
        rock_analysis.title("GebPy - Rock analysis")
        rock_analysis.geometry("1350x825")
        rock_analysis.resizable(False, False)
        rock_analysis["bg"] = self.colors_gebpy["Background"]

        ## Cleaning
        categories = ["Frame", "Label", "Radiobutton", "Entry", "Checkbox"]
        priorities = ["Static", "Temporary"]
        for priority in priorities:
            for category in categories:
                if len(self.gui_elements_sub["Mineralogy"][priority][category]) > 0:
                    self.gui_elements_sub["Mineralogy"][priority][category].clear()

        ## Geometry and Layout
        window_width = 1200
        window_heigth = 800
        row_min = 20
        n_rows = int(window_heigth/row_min)
        column_min = 20
        n_columns = int(window_width/column_min)

        for x in range(n_columns):
            tk.Grid.columnconfigure(rock_analysis, x, weight=1)
        for y in range(n_rows):
            tk.Grid.rowconfigure(rock_analysis, y, weight=1)

        # Rows
        for i in range(0, n_rows):
            rock_analysis.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            rock_analysis.grid_columnconfigure(i, minsize=column_min)

        ## FRAMES
        SimpleElements(parent=rock_analysis, row_id=0, column_id=0, n_rows=n_rows, n_columns=15,
                       bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Option"]).create_frame()

        ## LABELS
        lbl_01 = SimpleElements(
            parent=rock_analysis, row_id=0, column_id=0, n_rows=1, n_columns=15, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Rock analysis", font_option="sans 14 bold", relief=tk.FLAT)
        lbl_02 = SimpleElements(
            parent=rock_analysis, row_id=2, column_id=1, n_rows=1, n_columns=10, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Accent"]).create_label(
            text="Rock definition", font_option="sans 14 bold", relief=tk.FLAT, anchor_option=tk.W)
        lbl_03 = SimpleElements(
            parent=rock_analysis, row_id=6, column_id=1, n_rows=1, n_columns=10, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Accent"]).create_label(
            text="External data", font_option="sans 14 bold", relief=tk.FLAT, anchor_option=tk.W)
        lbl_04 = SimpleElements(
            parent=rock_analysis, row_id=10, column_id=1, n_rows=1, n_columns=10, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Accent"]).create_label(
            text="Plotting options", font_option="sans 14 bold", relief=tk.FLAT, anchor_option=tk.W)
        lbl_04a = SimpleElements(
            parent=rock_analysis, row_id=11, column_id=1, n_rows=1, n_columns=7, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Define x-axis", font_option="sans 12 bold", relief=tk.FLAT, anchor_option=tk.W)
        lbl_04b = SimpleElements(
            parent=rock_analysis, row_id=12, column_id=1, n_rows=1, n_columns=7, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Define y-axis", font_option="sans 12 bold", relief=tk.FLAT, anchor_option=tk.W)

        ## RADIOBUTTONS
        var_rb_02 = self.var_rb["Rock analysis"]["Rock definition"]
        var_rb_02.set(0)

        rb_02a = SimpleElements(
            parent=rock_analysis, row_id=3, column_id=1, n_rows=1, n_columns=9, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_radiobutton(
            text="Predefined rock", var_rb=var_rb_02, value_rb=0, color_bg=self.colors_gebpy["Accent"])
        rb_02b = SimpleElements(
            parent=rock_analysis, row_id=4, column_id=1, n_rows=1, n_columns=9, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_radiobutton(
            text="Custom rock", var_rb=var_rb_02, value_rb=1, color_bg=self.colors_gebpy["Accent"])

        rb_02a.configure(
            activebackground=self.colors_gebpy["Navigation"], activeforeground=self.colors_gebpy["Accent"],
            font="sans 12 bold")
        rb_02b.configure(
            activebackground=self.colors_gebpy["Navigation"], activeforeground=self.colors_gebpy["Accent"],
            font="sans 12 bold")

        var_rb_03 = self.var_rb["Rock analysis"]["External data"]
        var_rb_03.set(0)

        rb_03a = SimpleElements(
            parent=rock_analysis, row_id=7, column_id=1, n_rows=1, n_columns=9, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_radiobutton(
            text="Import data", var_rb=var_rb_03, value_rb=0, color_bg=self.colors_gebpy["Accent"])
        rb_03b = SimpleElements(
            parent=rock_analysis, row_id=8, column_id=1, n_rows=1, n_columns=9, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_radiobutton(
            text="Custom input", var_rb=var_rb_03, value_rb=1, color_bg=self.colors_gebpy["Accent"])

        rb_03a.configure(
            activebackground=self.colors_gebpy["Navigation"], activeforeground=self.colors_gebpy["Accent"],
            font="sans 12 bold")
        rb_03b.configure(
            activebackground=self.colors_gebpy["Navigation"], activeforeground=self.colors_gebpy["Accent"],
            font="sans 12 bold")

        ## BUTTONS
        btn_02 = SimpleElements(
            parent=rock_analysis, row_id=3, column_id=9, n_rows=2, n_columns=5, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(text="Setup")
        btn_03 = SimpleElements(
            parent=rock_analysis, row_id=7, column_id=9, n_rows=2, n_columns=5, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(text="Setup")
        btn_04 = SimpleElements(
            parent=rock_analysis, row_id=n_rows - 2, column_id=1, n_rows=1, n_columns=13, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(text="Export data")

        btn_02.configure(activebackground=self.colors_gebpy["Accent"], font="sans 12 bold")
        btn_03.configure(activebackground=self.colors_gebpy["Accent"], font="sans 12 bold")
        btn_04.configure(activebackground=self.colors_gebpy["Accent"], font="sans 12 bold")

        ## OPTION MENUS
        var_opt_04a = self.var_rb["Rock analysis"]["x-axis"]
        var_opt_04a.set("Select parameter")
        list_opt_04a = ["rho", "vP", "vS", "GR", "PE"]

        opt_04a = SimpleElements(
            parent=rock_analysis, row_id=11, column_id=9, n_rows=1, n_columns=5, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_option_menu(
            var_opt=var_opt_04a, var_opt_set=var_opt_04a.get(), opt_list=list_opt_04a,
            active_bg=self.colors_gebpy["Accent"])

        var_opt_04b = self.var_rb["Rock analysis"]["y-axis"]
        var_opt_04b.set("Select parameter")
        list_opt_04b = ["rho", "vP", "vS", "GR", "PE"]

        opt_04b = SimpleElements(
            parent=rock_analysis, row_id=12, column_id=9, n_rows=1, n_columns=5, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_option_menu(
            var_opt=var_opt_04b, var_opt_set=var_opt_04b.get(), opt_list=list_opt_04b,
            active_bg=self.colors_gebpy["Accent"])

        opt_04a.configure(font="sans 10 bold")
        opt_04b.configure(font="sans 10 bold")

    def rock_comparison(self):
        ## Cleaning
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Static"][category]:
                gui_element.grid_remove()
            for gui_element in self.gui_elements["Temporary"][category]:
                gui_element.grid_remove()
        #
        ## Labels
        lbl_title = SimpleElements(
            parent=self.parent, row_id=5, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Rock Comparison", font_option="sans 14 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].extend([lbl_title])

    def rock_builder_sedimentary(self):
        ## Window Settings
        window_width = 1800
        window_height = 900
        var_geometry = str(window_width) + "x" + str(window_height) + "+" + str(0) + "+" + str(0)

        row_min = 25
        n_rows = int(window_height/row_min)
        column_min = 20
        n_columns = int(window_width/column_min)
        self.n_rows_rbs = n_rows
        self.n_columns_rbs = n_columns

        str_title_window = "GebPy - Rock Builder (sedimentary)"
        self.subwindow_rockbuilder_sedimentary = tk.Toplevel(self.parent)
        self.subwindow_rockbuilder_sedimentary.title(str_title_window)
        self.subwindow_rockbuilder_sedimentary.geometry(var_geometry)
        self.subwindow_rockbuilder_sedimentary.resizable(False, False)
        self.subwindow_rockbuilder_sedimentary["bg"] = self.colors_gebpy["Navigation"]

        for x in range(n_columns):
            tk.Grid.columnconfigure(self.subwindow_rockbuilder_sedimentary, x, weight=1)
        for y in range(self.n_rows):
            tk.Grid.rowconfigure(self.subwindow_rockbuilder_sedimentary, y, weight=1)

        # Rows
        for i in range(0, self.n_rows):
            self.subwindow_rockbuilder_sedimentary.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            self.subwindow_rockbuilder_sedimentary.grid_columnconfigure(i, minsize=column_min)

        self.start_row = 0
        self.n_columns_setup = 11

        ## Frames
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row, column_id=self.n_columns_setup + 1,
            n_rows=self.n_rows, n_columns=n_columns - self.n_columns_setup - 1, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Background"]).create_frame()
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row,
            column_id=3*self.n_columns_setup + 6, n_rows=1, n_columns=self.n_columns_rbs - (3*self.n_columns_setup + 6),
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Navigation"]).create_frame()

        ## Labels
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Settings", relief=tk.FLAT, font_option="sans 12 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 1, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Qz", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 3, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Fsp (Kfs+Pl)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 5, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Kfs/(Kfs+Pl) ratio", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 7, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Cal+Dol", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 9, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Cal/Dol ratio", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 11, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Ilt+Mnt", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 13, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Ilt/Mnt ratio", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 15, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Anh+Gp", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 17, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Anh/Gp ratio", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 19, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Org", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 21, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Total amount of Py+Sd", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 23, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Py/Sd ratio", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 25, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Range of porosity", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 27, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["White"]).create_label(
            text="Number of datapoints", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row, column_id=self.n_columns_setup + 2,
            n_rows=1, n_columns=2*self.n_columns_setup + 3, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Black"]).create_label(
            text="Results", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row, column_id=3*self.n_columns_setup + 6,
            n_rows=1, n_columns=self.n_columns_rbs - (3*self.n_columns_setup + 6), bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Black"]).create_label(
            text="Diagrams", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)

        ## RADIOBUTTONS
        self.rb_01a = SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 1,
            column_id=3*self.n_columns_setup + 6, n_rows=1, n_columns=5, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Black"]).create_radiobutton(
            text="Histogram", var_rb=self.container_variables["Radiobuttons"]["Category Diagram"], value_rb=0,
            color_bg=self.colors_gebpy["Accent"],
            command=lambda var_rb=self.container_variables["Radiobuttons"]["Category Diagram"], key="Category Diagram":
            self.change_diagram_category(var_rb, key))
        self.rb_01b = SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 1,
            column_id=3*self.n_columns_setup + 11, n_rows=1, n_columns=5, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Black"]).create_radiobutton(
            text="Scatter", var_rb=self.container_variables["Radiobuttons"]["Category Diagram"], value_rb=1,
            color_bg=self.colors_gebpy["Accent"],
            command=lambda var_rb=self.container_variables["Radiobuttons"]["Category Diagram"], key="Category Diagram":
            self.change_diagram_category(var_rb, key))
        self.rb_02a = SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 1,
            column_id=3*self.n_columns_setup + 16, n_rows=1, n_columns=8, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Black"]).create_radiobutton(
            text="Geophysical data", var_rb=self.container_variables["Radiobuttons"]["Category Data"], value_rb=0,
            color_bg=self.colors_gebpy["Accent"],
            command=lambda var_rb=self.container_variables["Radiobuttons"]["Category Data"], key="Category Data":
            self.change_diagram_category(var_rb, key))
        self.rb_02b = SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 1,
            column_id=3*self.n_columns_setup + 24, n_rows=1, n_columns=8, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Black"]).create_radiobutton(
            text="Geochemical data", var_rb=self.container_variables["Radiobuttons"]["Category Data"], value_rb=1,
            color_bg=self.colors_gebpy["Accent"],
            command=lambda var_rb=self.container_variables["Radiobuttons"]["Category Data"], key="Category Data":
            self.change_diagram_category(var_rb, key))

        self.rb_01a.configure(state="disabled")
        self.rb_01b.configure(state="disabled")
        self.rb_02a.configure(state="disabled")
        self.rb_02b.configure(state="disabled")

        ## ENTRIES
        var_entr_01a = tk.StringVar()
        var_entr_01a.set("90.0")
        var_entr_01b = tk.StringVar()
        var_entr_01b.set("100.0")
        var_entr_02a = tk.StringVar()
        var_entr_02a.set("0.0")
        var_entr_02b = tk.StringVar()
        var_entr_02b.set("0.0")
        var_entr_03 = tk.StringVar()
        var_entr_03.set("0.0")
        var_entr_04a = tk.StringVar()
        var_entr_04a.set("0.0")
        var_entr_04b = tk.StringVar()
        var_entr_04b.set("0.0")
        var_entr_05 = tk.StringVar()
        var_entr_05.set("0.0")
        var_entr_06a = tk.StringVar()
        var_entr_06a.set("0.0")
        var_entr_06b = tk.StringVar()
        var_entr_06b.set("0.0")
        var_entr_07 = tk.StringVar()
        var_entr_07.set("0.0")
        var_entr_08a = tk.StringVar()
        var_entr_08a.set("0.0")
        var_entr_08b = tk.StringVar()
        var_entr_08b.set("0.0")
        var_entr_09 = tk.StringVar()
        var_entr_09.set("0.0")
        var_entr_10a = tk.StringVar()
        var_entr_10a.set("0.0")
        var_entr_10b = tk.StringVar()
        var_entr_10b.set("0.0")
        var_entr_11a = tk.StringVar()
        var_entr_11a.set("0.0")
        var_entr_11b = tk.StringVar()
        var_entr_11b.set("0.0")
        var_entr_12 = tk.StringVar()
        var_entr_12.set("0.0")
        var_entr_13a = tk.IntVar()
        var_entr_13a.set(0)
        var_entr_13b = tk.IntVar()
        var_entr_13b.set(10)
        var_entr_14 = tk.IntVar()
        var_entr_14.set(100)

        self.helper_rockbuilder_sedimentary_variables = {
            "w(Qz,min)": var_entr_01a, "w(Qz,max)": var_entr_01b, "w(Kfs+Pl,min)": var_entr_02a,
            "w(Kfs+Pl,max)": var_entr_02b, "w(Kfs)/w(Kfs+Pl)": var_entr_03, "w(Cal+Dol,min)": var_entr_04a,
            "w(Cal+Dol,max)": var_entr_04b, "w(Cal)/w(Cal+Dol)": var_entr_05, "w(Ilt+Mnt,min)": var_entr_06a,
            "w(Ilt+Mnt,max)": var_entr_06b, "w(Ilt)/w(Ilt+Mnt)": var_entr_07, "w(Anh+Gp,min)": var_entr_08a,
            "w(Anh+Gp,max)": var_entr_08b, "w(Anh)/w(Anh+Gp)": var_entr_09, "w(org,min)": var_entr_10a,
            "w(org,max)": var_entr_10b, "w(Py+Sd,min)": var_entr_11a, "w(Py+Sd,max)": var_entr_11b,
            "w(Py)/w(Py+Sd)": var_entr_12, "phi(min)": var_entr_13a, "phi(max)": var_entr_13b,
            "n(datapoints)": var_entr_14}

        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 2, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_01a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 2,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_01b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 4, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_02a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 4,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_02b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 6, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_03)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 8, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_04a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 8,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_04b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 10, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_05)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 12, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_06a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 12,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_06b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 14, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_07)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 16, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_08a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 16,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_08b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 18, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_09)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 20, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_10a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 20,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_10b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 22, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_11a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 22,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_11b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 24, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_12)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 26, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_13a)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 26,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_13b)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 28, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors_gebpy["White"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=var_entr_14)

        ## BUTTONS
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 30, column_id=1, n_rows=1,
            n_columns=int((self.n_columns_setup - 1)/2), bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Run calculation", command=self.run_calculation_rockbuilder_sedimentary)
        SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 30,
            column_id=int((self.n_columns_setup - 1)/2) + 1, n_rows=1, n_columns=int((self.n_columns_setup - 1)/2),
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(text="Export results")

        ## TREEVIEWS
        list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
        list_width = list(80*np.ones(len(list_categories)))
        list_width = [int(item) for item in list_width]
        list_width[0] = 100
        list_width[-1] = 135

        self.tv_rockbuilder_01 = SimpleElements(
            parent=self.subwindow_rockbuilder_sedimentary, row_id=self.start_row + 1,
            column_id=self.n_columns_setup + 2, n_rows=n_rows - 2, n_columns=2*self.n_columns_setup + 2,
            fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
            n_categories=len(list_categories), text_n=list_categories, width_n=list_width, individual=True)

        scb_v = ttk.Scrollbar(self.subwindow_rockbuilder_sedimentary, orient="vertical")
        scb_h = ttk.Scrollbar(self.subwindow_rockbuilder_sedimentary, orient="horizontal")
        self.tv_rockbuilder_01.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
        scb_v.config(command=self.tv_rockbuilder_01.yview)
        scb_h.config(command=self.tv_rockbuilder_01.xview)
        scb_v.grid(
            row=self.start_row + 1, column=3*self.n_columns_setup + 4, rowspan=n_rows - 2, columnspan=1, sticky="ns")
        scb_h.grid(
            row=n_rows - 1, column=self.n_columns_setup + 2, rowspan=1, columnspan=2*self.n_columns_setup + 2,
            sticky="ew")

    def run_calculation_rockbuilder_sedimentary(self):
        self.rb_01a.configure(state="normal")
        self.rb_01b.configure(state="normal")
        self.rb_02a.configure(state="normal")
        self.rb_02b.configure(state="normal")

        self.helper_results_geophysics = {
            "rho": [], "rho(solid)": [], "V(molar)": [], "vP": [], "vS": [], "vP/vS": [], "GR": [], "PE": [], "phi": [],
            "K": [], "G": [], "E": [], "M": [], "lame": [], "poisson": []}
        self.helper_results_minerals = {
            "Qz": [], "Kfs": [], "Pl": [], "Cal": [], "Dol": [], "Ilt": [], "Mnt": [], "Gp": [], "Anh": [],
            "Org": [], "Py": [], "Sd": []}

        for key, item in self.helper_rockbuilder_sedimentary_variables.items():
            if key == "w(Kfs)/w(Kfs+Pl)":
                val_ratio_KfsPl = float(item.get())
            elif key == "w(Cal)/w(Cal+Dol)":
                val_ratio_CalDol = float(item.get())
            elif key == "w(Ilt)/w(Ilt+Mnt)":
                val_ratio_IltMnt = float(item.get())
            elif key == "w(Anh)/w(Anh+Gp)":
                val_ratio_AnhGp = float(item.get())
            elif key == "w(Py)/w(Py+Sd)":
                val_ratio_PySd = float(item.get())
            elif key == "phi(min)":
                val_phiMin = float(item.get())
            elif key == "phi(max)":
                val_phiMax = float(item.get())
            elif key == "n(datapoints)":
                val_nDatapoints = int(item.get())

        ## Results
        data_quartz = Oxides(mineral="Quartz", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_calcite = Carbonates(mineral="Calcite", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_dolomite = Carbonates(mineral="Dolomite", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_illite = Phyllosilicates(mineral="Illite", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_anhydrite = Sulfates(mineral="Anhydrite", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_gypsum = Sulfates(mineral="Gypsum", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_pyrite = Sulfides(mineral="Pyrite", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_siderite = Carbonates(mineral="Siderite", data_type=True, traces_list=[]).generate_dataset(number=1)
        data_water = Water.water("")

        helper_results = {
            "rho(solid)": [], "rho": [], "\u03C6": [], "V(molar)": [], "K": [], "G": [], "E": [], "\u03BB": [],
            "\u03BD": [], "M": [], "vP": [], "vS": [], "vP/vS": [], "GR": [], "PE": []}
        helper_results_minerals = {}
        helper_results_chemistry = {}
        helper_results_elements = {
            "O": [], "Si": [], "Na": [], "Al": [], "K": [], "Ca": [], "C": [], "Mg": [], "H": [], "S": [], "N": [],
            "Fe": []}
        helper_results_oxides = {
            "SiO2": [], "Na2O": [], "Al2O3": [], "K2O": [], "CaO": [], "CO2": [], "MgO": [], "H2O": [], "SO2": [],
            "NO2": [], "Fe2O3": []}

        for index in range(val_nDatapoints):
            phi = round(np.random.uniform(val_phiMin, val_phiMax)/100, 4)
            helper_rho_s = []
            helper_vol = []
            helper_bulkmod = []
            helper_shearmod = []
            helper_poisson = []
            helper_gr = []
            helper_pe = []
            helper_u = []
            helper_results_chemistry = {}
            ## Mineralogy
            # Amounts
            amounts_minerals = self.find_mineral_amounts()
            w_qz = amounts_minerals["Qz"]/100
            w_fsp = round(amounts_minerals["Fsp"]/100, 4)
            w_carb = round(amounts_minerals["Carb"]/100, 4)
            w_clay = round(amounts_minerals["Clay"]/100, 4)
            w_sulfat = round(amounts_minerals["Sulfat"]/100, 4)
            w_org = round(amounts_minerals["Org"]/100, 4)
            w_sulfid = round(amounts_minerals["Sulfid"]/100, 4)
            # Quartz
            helper_rho_s.append(w_qz*data_quartz["rho"][0])
            helper_vol.append(w_qz*data_quartz["V"][0])
            helper_bulkmod.append(w_qz*data_quartz["K"][0])
            helper_shearmod.append(w_qz*data_quartz["G"][0])
            helper_poisson.append(w_qz*data_quartz["nu"][0])
            helper_gr.append(w_qz*data_quartz["GR"][0])
            helper_pe.append(w_qz*data_quartz["PE"][0])
            helper_u.append(w_qz*data_quartz["U"][0])
            if "Qz" not in helper_results_minerals:
                helper_results_minerals["Qz"] = []
            helper_results_minerals["Qz"].append(100*w_qz)
            self.helper_results_minerals["Qz"].append(100*w_qz)
            for element, value in data_quartz["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_qz > 0:
                    helper_results_chemistry[element].append(w_qz*value[0])
                    if round(w_qz, 4) == 1.0:
                        print(w_qz, w_qz*value[0])
            # Feldspar minerals
            # Alkaline feldspar
            w_kfs = round(val_ratio_KfsPl*w_fsp, 4)
            data_alkaline_feldspar = Tectosilicates(
                mineral="Alkali Feldspar", data_type=True, traces_list=[]).generate_dataset(number=1)
            helper_rho_s.append(w_kfs*data_alkaline_feldspar["rho"][0])
            helper_vol.append(w_kfs*data_alkaline_feldspar["V"][0])
            helper_bulkmod.append(w_kfs*data_alkaline_feldspar["K"][0])
            helper_shearmod.append(w_kfs*data_alkaline_feldspar["G"][0])
            helper_poisson.append(w_kfs*data_alkaline_feldspar["nu"][0])
            helper_gr.append(w_kfs*data_alkaline_feldspar["GR"][0])
            helper_pe.append(w_kfs*data_alkaline_feldspar["PE"][0])
            helper_u.append(w_kfs*data_alkaline_feldspar["U"][0])
            if "Kfs" not in helper_results_minerals:
                helper_results_minerals["Kfs"] = []
            helper_results_minerals["Kfs"].append(100*w_kfs)
            self.helper_results_minerals["Kfs"].append(100*w_kfs)
            for element, value in data_alkaline_feldspar["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_kfs > 0:
                    helper_results_chemistry[element].append(w_kfs*value[0])
            # Plagioclase
            w_pl = round(w_fsp - w_kfs, 4)
            data_plagioclase = Tectosilicates(
                mineral="Plagioclase", data_type=True, traces_list=[]).generate_dataset(number=1)
            helper_rho_s.append(w_pl*data_plagioclase["rho"][0])
            helper_vol.append(w_pl*data_plagioclase["V"][0])
            helper_bulkmod.append(w_pl*data_plagioclase["K"][0])
            helper_shearmod.append(w_pl*data_plagioclase["G"][0])
            helper_poisson.append(w_pl*data_plagioclase["nu"][0])
            helper_gr.append(w_pl*data_plagioclase["GR"][0])
            helper_pe.append(w_pl*data_plagioclase["PE"][0])
            helper_u.append(w_pl*data_plagioclase["U"][0])
            if "Pl" not in helper_results_minerals:
                helper_results_minerals["Pl"] = []
            helper_results_minerals["Pl"].append(100*w_pl)
            self.helper_results_minerals["Pl"].append(100*w_pl)
            for element, value in data_plagioclase["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_pl > 0:
                    helper_results_chemistry[element].append(w_pl*value[0])
            # Carbonate minerals
            # Calcite
            w_cal = round(val_ratio_CalDol*w_carb, 4)
            helper_rho_s.append(w_cal*data_calcite["rho"][0])
            helper_vol.append(w_cal*data_calcite["V"][0])
            helper_bulkmod.append(w_cal*data_calcite["K"][0])
            helper_shearmod.append(w_cal*data_calcite["G"][0])
            helper_poisson.append(w_cal*data_calcite["nu"][0])
            helper_gr.append(w_cal*data_calcite["GR"][0])
            helper_pe.append(w_cal*data_calcite["PE"][0])
            helper_u.append(w_cal*data_calcite["U"][0])
            if "Cal" not in helper_results_minerals:
                helper_results_minerals["Cal"] = []
            helper_results_minerals["Cal"].append(100*w_cal)
            self.helper_results_minerals["Cal"].append(100*w_cal)
            for element, value in data_calcite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_cal > 0:
                    helper_results_chemistry[element].append(w_cal*value[0])
            # Dolomite
            w_dol = round(w_carb - w_cal, 4)
            helper_rho_s.append(w_dol*data_dolomite["rho"][0])
            helper_vol.append(w_dol*data_dolomite["V"][0])
            helper_bulkmod.append(w_dol*data_dolomite["K"][0])
            helper_shearmod.append(w_dol*data_dolomite["G"][0])
            helper_poisson.append(w_dol*data_dolomite["nu"][0])
            helper_gr.append(w_dol*data_dolomite["GR"][0])
            helper_pe.append(w_dol*data_dolomite["PE"][0])
            helper_u.append(w_dol*data_dolomite["U"][0])
            if "Dol" not in helper_results_minerals:
                helper_results_minerals["Dol"] = []
            helper_results_minerals["Dol"].append(100*w_dol)
            self.helper_results_minerals["Dol"].append(100*w_dol)
            for element, value in data_dolomite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_dol > 0:
                    helper_results_chemistry[element].append(w_dol*value[0])
            # Clay minerals
            # Illite
            w_ilt = round(val_ratio_IltMnt*w_clay, 4)
            helper_rho_s.append(w_ilt*data_illite["rho"][0])
            helper_vol.append(w_ilt*data_illite["V"][0])
            helper_bulkmod.append(w_ilt*data_illite["K"][0])
            helper_shearmod.append(w_ilt*data_illite["G"][0])
            helper_poisson.append(w_ilt*data_illite["nu"][0])
            helper_gr.append(w_ilt*data_illite["GR"][0])
            helper_pe.append(w_ilt*data_illite["PE"][0])
            helper_u.append(w_ilt*data_illite["U"][0])
            if "Ilt" not in helper_results_minerals:
                helper_results_minerals["Ilt"] = []
            helper_results_minerals["Ilt"].append(100*w_ilt)
            self.helper_results_minerals["Ilt"].append(100*w_ilt)
            for element, value in data_illite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_ilt > 0:
                    helper_results_chemistry[element].append(w_ilt*value[0])
            # Montmorillonite
            w_mnt = round(w_clay - w_ilt, 4)
            data_montmorillonite = Phyllosilicates(
                mineral="Montmorillonite", data_type=True, traces_list=[]).generate_dataset(number=1)
            helper_rho_s.append(w_mnt*data_montmorillonite["rho"][0])
            helper_vol.append(w_mnt*data_montmorillonite["V"][0])
            helper_bulkmod.append(w_mnt*data_montmorillonite["K"][0])
            helper_shearmod.append(w_mnt*data_montmorillonite["G"][0])
            helper_poisson.append(w_mnt*data_montmorillonite["nu"][0])
            helper_gr.append(w_mnt*data_montmorillonite["GR"][0])
            helper_pe.append(w_mnt*data_montmorillonite["PE"][0])
            helper_u.append(w_mnt*data_montmorillonite["U"][0])
            if "Mnt" not in helper_results_minerals:
                helper_results_minerals["Mnt"] = []
            helper_results_minerals["Mnt"].append(100*w_mnt)
            self.helper_results_minerals["Mnt"].append(100*w_mnt)
            for element, value in data_montmorillonite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_mnt > 0:
                    helper_results_chemistry[element].append(w_mnt*value[0])
            # Sulfate minerals
            # Anhydrite
            w_anh = round(val_ratio_AnhGp*w_sulfat, 4)
            helper_rho_s.append(w_anh*data_anhydrite["rho"][0])
            helper_vol.append(w_anh*data_anhydrite["V"][0])
            helper_bulkmod.append(w_anh*data_anhydrite["K"][0])
            helper_shearmod.append(w_anh*data_anhydrite["G"][0])
            helper_poisson.append(w_anh*data_anhydrite["nu"][0])
            helper_gr.append(w_anh*data_anhydrite["GR"][0])
            helper_pe.append(w_anh*data_anhydrite["PE"][0])
            helper_u.append(w_anh*data_anhydrite["U"][0])
            if "Anh" not in helper_results_minerals:
                helper_results_minerals["Anh"] = []
            helper_results_minerals["Anh"].append(100*w_anh)
            self.helper_results_minerals["Anh"].append(100*w_anh)
            for element, value in data_anhydrite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_anh > 0:
                    helper_results_chemistry[element].append(w_anh*value[0])
            # Gypsum
            w_gp = round(w_sulfat - w_anh, 4)
            helper_rho_s.append(w_gp*data_gypsum["rho"][0])
            helper_vol.append(w_gp*data_gypsum["V"][0])
            helper_bulkmod.append(w_gp*data_gypsum["K"][0])
            helper_shearmod.append(w_gp*data_gypsum["G"][0])
            helper_poisson.append(w_gp*data_gypsum["nu"][0])
            helper_gr.append(w_gp*data_gypsum["GR"][0])
            helper_pe.append(w_gp*data_gypsum["PE"][0])
            helper_u.append(w_gp*data_gypsum["U"][0])
            if "Gp" not in helper_results_minerals:
                helper_results_minerals["Gp"] = []
            helper_results_minerals["Gp"].append(100*w_gp)
            self.helper_results_minerals["Gp"].append(100*w_gp)
            for element, value in data_gypsum["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_gp > 0:
                    helper_results_chemistry[element].append(w_gp*value[0])
            # Organic matter
            data_organic_matter = Organics(
                mineral="Organic Matter", data_type=True, traces_list=[]).generate_dataset(number=1)
            helper_rho_s.append(w_org*data_organic_matter["rho"][0])
            helper_vol.append(w_org*data_organic_matter["V"][0])
            helper_bulkmod.append(w_org*data_organic_matter["K"][0])
            helper_shearmod.append(w_org*data_organic_matter["G"][0])
            helper_poisson.append(w_org*data_organic_matter["nu"][0])
            helper_gr.append(w_org*data_organic_matter["GR"][0])
            helper_pe.append(w_org*data_organic_matter["PE"][0])
            helper_u.append(w_org*data_organic_matter["U"][0])
            if "Org" not in helper_results_minerals:
                helper_results_minerals["Org"] = []
            helper_results_minerals["Org"].append(100*w_org)
            self.helper_results_minerals["Org"].append(100*w_org)
            for element, value in data_organic_matter["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_org > 0:
                    helper_results_chemistry[element].append(w_org*value[0])
            # Sulfide minerals
            # Pyrite
            w_py = round(val_ratio_PySd*w_sulfid, 4)
            helper_rho_s.append(w_py*data_pyrite["rho"][0])
            helper_vol.append(w_py*data_pyrite["V"][0])
            helper_bulkmod.append(w_py*data_pyrite["K"][0])
            helper_shearmod.append(w_py*data_pyrite["G"][0])
            helper_poisson.append(w_py*data_pyrite["nu"][0])
            helper_gr.append(w_py*data_pyrite["GR"][0])
            helper_pe.append(w_py*data_pyrite["PE"][0])
            helper_u.append(w_py*data_pyrite["U"][0])
            if "Py" not in helper_results_minerals:
                helper_results_minerals["Py"] = []
            helper_results_minerals["Py"].append(100*w_py)
            self.helper_results_minerals["Py"].append(100*w_py)
            for element, value in data_pyrite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_py > 0:
                    helper_results_chemistry[element].append(w_py*value[0])
            # Siderite
            w_sd = round(w_sulfid - w_py, 4)
            helper_rho_s.append(w_sd*data_siderite["rho"][0])
            helper_vol.append(w_sd*data_siderite["V"][0])
            helper_bulkmod.append(w_sd*data_siderite["K"][0])
            helper_shearmod.append(w_sd*data_siderite["G"][0])
            helper_poisson.append(w_sd*data_siderite["nu"][0])
            helper_gr.append(w_sd*data_siderite["GR"][0])
            helper_pe.append(w_sd*data_siderite["PE"][0])
            helper_u.append(w_sd*data_siderite["U"][0])
            if "Sd" not in helper_results_minerals:
                helper_results_minerals["Sd"] = []
            helper_results_minerals["Sd"].append(100*w_sd)
            self.helper_results_minerals["Sd"].append(100*w_sd)
            for element, value in data_siderite["chemistry"].items():
                if element not in helper_results_chemistry:
                    helper_results_chemistry[element] = []
                if w_sd > 0:
                    helper_results_chemistry[element].append(w_sd*value[0])
            ## Geophysical results
            # Porosity
            helper_results["\u03C6"].append(phi)
            self.helper_results_geophysics["phi"].append(phi)
            # Bulk density (solid)
            rho_s = round(np.sum(helper_rho_s), 3)
            helper_results["rho(solid)"].append(rho_s)
            self.helper_results_geophysics["rho(solid)"].append(rho_s)
            # Bulk density (solid)
            rho = round((1 - phi)*rho_s + phi*data_water[2], 3)
            helper_results["rho"].append(rho)
            self.helper_results_geophysics["rho"].append(rho)
            # Molar volume
            molar_vol = round(np.sum(helper_vol), 3)
            helper_results["V(molar)"].append(molar_vol)
            self.helper_results_geophysics["V(molar)"].append(molar_vol)
            # Elastic properties
            bulkmod = round(np.sum(helper_bulkmod), 3)
            shearmod = round(np.sum(helper_shearmod), 3)
            youngs_mod = round((9*bulkmod*shearmod)/(3*bulkmod + shearmod), 3)
            lame = round(bulkmod - 2/3*shearmod, 3)
            vpmod = round(bulkmod + 4/3*shearmod, 3)
            poisson = round((3*bulkmod - 2*shearmod)/(2*(3*bulkmod + shearmod)), 3)
            helper_results["K"].append(bulkmod)
            helper_results["G"].append(shearmod)
            helper_results["E"].append(youngs_mod)
            helper_results["\u03BB"].append(lame)
            helper_results["\u03BD"].append(poisson)
            helper_results["M"].append(vpmod)
            self.helper_results_geophysics["K"].append(bulkmod)
            self.helper_results_geophysics["G"].append(shearmod)
            self.helper_results_geophysics["E"].append(youngs_mod)
            self.helper_results_geophysics["M"].append(vpmod)
            self.helper_results_geophysics["lame"].append(lame)
            self.helper_results_geophysics["poisson"].append(poisson)
            # Seismic velocities
            v_p = round(((bulkmod*10**9 + 4/3*shearmod*10**9)/(rho))**0.5, 3)
            v_s = round(((shearmod*10**9)/(rho))**0.5, 3)
            v_ps_ratio = round(v_p/v_s, 3)
            helper_results["vP"].append(v_p)
            helper_results["vS"].append(v_s)
            helper_results["vP/vS"].append(v_ps_ratio)
            self.helper_results_geophysics["vP"].append(v_p)
            self.helper_results_geophysics["vS"].append(v_s)
            self.helper_results_geophysics["vP/vS"].append(v_ps_ratio)
            # Gamma ray
            gr = round(np.sum(helper_gr), 3)
            helper_results["GR"].append(gr)
            self.helper_results_geophysics["GR"].append(gr)
            # Photoelectricity
            pe = round(np.sum(helper_pe), 3)
            helper_results["PE"].append(pe)
            self.helper_results_geophysics["PE"].append(pe)
            ## Geochemical results
            # Element concentrations
            # Oxygen
            val_o = round(np.sum(helper_results_chemistry["O"]), 3)
            helper_results_elements["O"].append(val_o)
            # Silicon
            val_si = np.sum(helper_results_chemistry["Si"])
            helper_results_elements["Si"].append(val_si)
            # Sodium
            val_na = round(np.sum(helper_results_chemistry["Na"]), 3)
            helper_results_elements["Na"].append(val_na)
            # Aluminium
            val_al = round(np.sum(helper_results_chemistry["Al"]), 3)
            helper_results_elements["Al"].append(val_al)
            # Potassium
            val_k = round(np.sum(helper_results_chemistry["K"]), 3)
            helper_results_elements["K"].append(val_k)
            # Calcium
            val_ca = round(np.sum(helper_results_chemistry["Ca"]), 3)
            helper_results_elements["Ca"].append(val_ca)
            # Carbon
            val_c = round(np.sum(helper_results_chemistry["C"]), 3)
            helper_results_elements["C"].append(val_c)
            # Magnesium
            val_mg = round(np.sum(helper_results_chemistry["Mg"]), 3)
            helper_results_elements["Mg"].append(val_mg)
            # Hydrogen
            val_h = round(np.sum(helper_results_chemistry["H"]), 3)
            helper_results_elements["H"].append(val_h)
            # Sulfur
            val_s = round(np.sum(helper_results_chemistry["S"]), 3)
            helper_results_elements["S"].append(val_s)
            # Nitrogen
            val_n = round(np.sum(helper_results_chemistry["N"]), 3)
            helper_results_elements["N"].append(val_n)
            # Iron
            val_fe = round(np.sum(helper_results_chemistry["Fe"]), 3)
            helper_results_elements["Fe"].append(val_fe)
            # Oxide concentrations
            # SiO2
            value = val_si*self.conversion_factors["SiO2"]
            if round(w_qz, 4) == 1.0:
                print(w_qz, val_si, self.conversion_factors["SiO2"], value)
            # if round(w_qz, 4) == 1.0:
            #     value = 1.0
            helper_results_oxides["SiO2"].append(value)
            # Na2O
            value = val_na*self.conversion_factors["Na2O"]
            helper_results_oxides["Na2O"].append(value)
            # Al2O3
            value = val_al*self.conversion_factors["Al2O3"]
            helper_results_oxides["Al2O3"].append(value)
            # K2O
            value = val_k*self.conversion_factors["K2O"]
            helper_results_oxides["K2O"].append(value)
            # CaO
            value = val_ca*self.conversion_factors["CaO"]
            helper_results_oxides["CaO"].append(value)
            # CO2
            value = val_c*self.conversion_factors["CO2"]
            helper_results_oxides["CO2"].append(value)
            # MgO
            value = val_mg*self.conversion_factors["MgO"]
            helper_results_oxides["MgO"].append(value)
            # H2O
            value = val_h*self.conversion_factors["H2O"]
            helper_results_oxides["H2O"].append(value)
            # SO2
            value = val_s*self.conversion_factors["SO2"]
            helper_results_oxides["SO2"].append(value)
            # Fe2O3
            value = val_fe*self.conversion_factors["Fe2O3"]
            helper_results_oxides["Fe2O3"].append(value)

        if len(self.tv_rockbuilder_01.get_children()) > 0:
            for item in self.tv_rockbuilder_01.get_children():
                self.tv_rockbuilder_01.delete(item)

        for category, dataset in helper_results.items():
            val_min = round(min(dataset), 3)
            val_max = round(max(dataset), 3)
            val_mean = round(np.mean(dataset), 3)
            val_std = round(np.std(dataset, ddof=1), 3)
            entries = [category, val_min, val_max, val_mean, val_std]
            self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        entries = ["---", "---", "---", "---", "---"]
        self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        for category, dataset in helper_results_minerals.items():
            val_min = round(min(dataset), 2)
            val_max = round(max(dataset), 2)
            val_mean = round(np.mean(dataset), 2)
            val_std = round(np.std(dataset, ddof=1), 2)
            entries = [category, val_min, val_max, val_mean, val_std]
            self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        entries = ["---", "---", "---", "---", "---"]
        self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        for category, dataset in helper_results_oxides.items():
            if len(dataset) > 0:
                val_min = round(min(dataset)*100, 2)
                val_max = round(max(dataset)*100, 2)
                val_mean = round(np.mean(dataset)*100, 2)
                val_std = round(np.std(dataset, ddof=1)*100, 2)
                entries = [category, val_min, val_max, val_mean, val_std]
                self.tv_rockbuilder_01.insert("", tk.END, values=entries)
            else:
                val_min = "-"
                val_max = "-"
                val_mean = "-"
                val_std = "-"
                entries = [category, val_min, val_max, val_mean, val_std]
                self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        entries = ["---", "---", "---", "---", "---"]
        self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        for category, dataset in helper_results_elements.items():
            if len(dataset) > 0:
                val_min = round(min(dataset)*100, 2)
                val_max = round(max(dataset)*100, 2)
                val_mean = round(np.mean(dataset)*100, 2)
                val_std = round(np.std(dataset*100, ddof=1), 2)
                entries = [category, val_min, val_max, val_mean, val_std]
                self.tv_rockbuilder_01.insert("", tk.END, values=entries)
            else:
                val_min = "-"
                val_max = "-"
                val_mean = "-"
                val_std = "-"
                entries = [category, val_min, val_max, val_mean, val_std]
                self.tv_rockbuilder_01.insert("", tk.END, values=entries)

        ## Initialization
        self.container_variables["Radiobuttons"]["Category Diagram"].set(0)
        self.container_variables["Radiobuttons"]["Category Data"].set(0)

        categories = [["rho", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "poisson"]]
        labels = [["kg/m^3", "API", "barns/e-"], ["m/s", "m/s", "1"], ["GPa", "GPa", "1"]]
        self.create_3x3_diagram(var_categories=categories, var_labels=labels)

    def change_diagram_category(self, var_rb, key):
        self.container_variables["Radiobuttons"][key].set(var_rb.get())
        if self.container_variables["Radiobuttons"]["Category Data"].get() == 0:    # Geophysical data
            categories = [["rho", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "poisson"]]
            labels = [["kg/m^3", "API", "barns/e-"], ["m/s", "m/s", "1"], ["GPa", "GPa", "1"]]

            self.create_3x3_diagram(var_categories=categories, var_labels=labels)
        else:                                                                       # Geochemical data
            categories = [["Qz", "Kfs", "Pl"], ["Cal", "Dol", "Ilt"], ["Mnt", "Gp", "Anh"], ["Org", "Py", "Sd"]]
            labels = [["wt.%", "wt.%", "wt.%"], ["wt.%", "wt.%", "wt.%"], ["wt.%", "wt.%", "wt.%"],
                      ["wt.%", "wt.%", "wt.%"]]

            max_avg = 0
            for mineral, values in sorted(self.helper_results_minerals.items()):
                avg = np.mean(values)
                if avg >= max_avg:
                    max_name = mineral
                    max_avg = avg

            self.create_4x3_diagram(var_categories=categories, var_labels=labels, max_key=max_name)

    def create_4x3_diagram(self, var_categories, var_labels, max_key="Qz"):
        self.fig_4x3 = Figure(figsize=(6, 6), tight_layout=True, facecolor=self.colors_gebpy["Background"])
        self.ax_4x3 = self.fig_4x3.subplots(nrows=4, ncols=3)

        min_max_key = np.min(self.helper_results_minerals[max_key])//5*5
        max_max_key = (np.max(self.helper_results_minerals[max_key]) + 5 - 1)//5*5
        if max_max_key > 100:
            max_max_key = 100

        delta_max_key = max_max_key - min_max_key

        if delta_max_key < 11:
            if delta_max_key < 6:
                stepsize_max_key = 1
            else:
                stepsize_max_key = 2
        elif delta_max_key < 21:
            if delta_max_key < 16:
                stepsize_max_key = 3
            else:
                stepsize_max_key = 4
        else:
            stepsize_max_key = 5

        lower_limit_x = min_max_key - stepsize_max_key
        if lower_limit_x < 0:
            lower_limit_x = 0

        upper_limit_x = max_max_key + stepsize_max_key
        if upper_limit_x > 100:
            upper_limit_x = 100

        for i, subcategories in enumerate(var_categories):
            for j, subcategory in enumerate(subcategories):
                if self.container_variables["Radiobuttons"]["Category Diagram"].get() == 0:    # Histogram
                    if subcategory in self.helper_results_minerals:
                        self.ax_4x3[i][j].hist(
                            self.helper_results_minerals[subcategory], bins=16, edgecolor="black",
                            facecolor=self.colors_gebpy["Option"])

                    self.ax_4x3[i][j].set_xlabel(subcategory + " (" + var_labels[i][j] + ")", fontsize=9)
                    self.ax_4x3[i][j].set_ylabel("Frequency (#)", labelpad=0.5, fontsize=9)
                elif self.container_variables["Radiobuttons"]["Category Diagram"].get() == 1:  # Scatter
                    if subcategory in self.helper_results_minerals:
                        self.ax_4x3[i][j].scatter(
                            self.helper_results_minerals[max_key], self.helper_results_minerals[subcategory],
                            edgecolor="black", color=self.colors_gebpy["Option"], alpha=0.75, s=50)

                    self.ax_4x3[i][j].set_xlabel(max_key + " (" + var_labels[i][j] + ")", fontsize=9)
                    self.ax_4x3[i][j].set_ylabel(subcategory + " (" + var_labels[i][j] + ")", labelpad=0.5, fontsize=9)

                    min_subcategory = np.min(self.helper_results_minerals[subcategory])//5*5
                    max_subcategory = (np.max(self.helper_results_minerals[subcategory]) + 5 - 1)//5*5

                    if min_subcategory < 0:
                        min_subcategory = 0
                    if max_subcategory > 100:
                        max_subcategory = 100

                    delta_subcategory = max_subcategory - min_subcategory

                    if delta_subcategory < 11:
                        if delta_subcategory < 6:
                            stepsize = 1
                        else:
                            stepsize = 2
                    elif delta_subcategory < 21:
                        if delta_subcategory < 16:
                            stepsize = 3
                        else:
                            stepsize = 4
                    else:
                        stepsize = 5

                    lower_limit_subcategory = min_subcategory - stepsize
                    if lower_limit_subcategory < 0:
                        lower_limit_subcategory = 0

                    upper_limit_subcategory = max_subcategory + stepsize
                    if upper_limit_subcategory > 100:
                        upper_limit_subcategory = 100

                    self.ax_4x3[i][j].set_xlim(lower_limit_x, upper_limit_x)
                    self.ax_4x3[i][j].set_xticks(
                        np.arange(lower_limit_x, upper_limit_x + stepsize_max_key, stepsize_max_key))
                    self.ax_4x3[i][j].set_ylim(lower_limit_subcategory, upper_limit_subcategory)
                    self.ax_4x3[i][j].set_yticks(
                        np.arange(lower_limit_subcategory, upper_limit_subcategory + stepsize, stepsize))

                self.ax_4x3[i][j].grid(True)
                self.ax_4x3[i][j].set_axisbelow(True)

        self.canvas_4x3 = FigureCanvasTkAgg(self.fig_4x3, master=self.subwindow_rockbuilder_sedimentary)
        self.canvas_4x3.get_tk_widget().grid(
            row=self.start_row + 2, column=3*self.n_columns_setup + 6, rowspan=self.n_rows_rbs - 2,
            columnspan=self.n_columns_rbs - (3*self.n_columns_setup + 6), sticky="nesw")

    def create_3x3_diagram(self, var_categories, var_labels):
        self.fig_3x3 = Figure(figsize=(6, 6), tight_layout=True, facecolor=self.colors_gebpy["Background"])
        self.ax_3x3 = self.fig_3x3.subplots(nrows=3, ncols=3)

        for i, subcategories in enumerate(var_categories):
            for j, subcategory in enumerate(subcategories):
                if self.container_variables["Radiobuttons"]["Category Diagram"].get() == 0:    # Histogram
                    if subcategory in self.helper_results_geophysics:
                        self.ax_3x3[i][j].hist(
                            self.helper_results_geophysics[subcategory], bins=16, edgecolor="black",
                            facecolor=self.colors_gebpy["Option"])

                    self.ax_3x3[i][j].set_xlabel(subcategory + " (" + var_labels[i][j] + ")", fontsize=9)
                    self.ax_3x3[i][j].set_ylabel("Frequency (#)", labelpad=0.5, fontsize=9)
                elif self.container_variables["Radiobuttons"]["Category Diagram"].get() == 1:  # Scatter
                    if subcategory in self.helper_results_geophysics:
                        if subcategory == "rho":
                            subcategory = "V(molar)"
                        self.ax_3x3[i][j].scatter(
                            self.helper_results_geophysics["rho"], self.helper_results_geophysics[subcategory],
                            edgecolor="black", color=self.colors_gebpy["Option"], alpha=0.75, s=50)

                    self.ax_3x3[i][j].set_xlabel("rho (kg/m^3)", fontsize=9)
                    self.ax_3x3[i][j].set_ylabel(subcategory + " (" + var_labels[i][j] + ")", labelpad=0.5, fontsize=9)

                self.ax_3x3[i][j].grid(True)
                self.ax_3x3[i][j].set_axisbelow(True)

        self.canvas_3x3 = FigureCanvasTkAgg(self.fig_3x3, master=self.subwindow_rockbuilder_sedimentary)
        self.canvas_3x3.get_tk_widget().grid(
            row=self.start_row + 2, column=3*self.n_columns_setup + 6, rowspan=self.n_rows_rbs - 2,
            columnspan=self.n_columns_rbs - (3*self.n_columns_setup + 6), sticky="nesw")

    def find_mineral_amounts(self):
        for key, item in self.helper_rockbuilder_sedimentary_variables.items():
            if key == "w(Qz,min)":
                val_wQz_min = float(item.get())
            elif key == "w(Qz,max)":
                val_wQz_max = float(item.get())
            elif key == "w(Kfs+Pl,min)":
                val_wKfsPl_min = float(item.get())
            elif key == "w(Kfs+Pl,max)":
                val_wKfsPl_max = float(item.get())
            elif key == "w(Cal+Dol,min)":
                val_wCalDol_min = float(item.get())
            elif key == "w(Cal+Dol,max)":
                val_wCalDol_max = float(item.get())
            elif key == "w(Ilt+Mnt,min)":
                val_wIltMnt_min = float(item.get())
            elif key == "w(Ilt+Mnt,max)":
                val_wIltMnt_max = float(item.get())
            elif key == "w(Anh+Gp,min)":
                val_wAnhGp_min = float(item.get())
            elif key == "w(Anh+Gp,max)":
                val_wAnhGp_max = float(item.get())
            elif key == "w(org,min)":
                val_wOrg_min = float(item.get())
            elif key == "w(org,max)":
                val_wOrg_max = float(item.get())
            elif key == "w(Py+Sd,min)":
                val_wPySd_min = float(item.get())
            elif key == "w(Py+Sd,max)":
                val_wPySd_max = float(item.get())

        w_qz_mean = (val_wQz_min + val_wQz_max)/2
        w_fsp_mean = (val_wKfsPl_min + val_wKfsPl_max)/2
        w_carb_mean = (val_wCalDol_min + val_wCalDol_max)/2
        w_clay_mean = (val_wIltMnt_min + val_wIltMnt_max)/2
        w_sulfat_mean = (val_wAnhGp_min + val_wAnhGp_max)/2
        w_org_mean = (val_wOrg_min + val_wOrg_max)/2
        w_sulfid_mean = (val_wPySd_min + val_wPySd_max)/2

        helper_data = {
            "Qz": [val_wQz_min, val_wQz_max], "Fsp": [val_wKfsPl_min, val_wKfsPl_max],
            "Carb": [val_wCalDol_min, val_wCalDol_max], "Clay": [val_wIltMnt_min, val_wIltMnt_max],
            "Sulfat": [val_wAnhGp_min, val_wAnhGp_max], "Org": [val_wOrg_min, val_wOrg_max],
            "Sulfid": [val_wPySd_min, val_wPySd_max]}
        dict_mean = {"Qz": w_qz_mean, "Fsp": w_fsp_mean, "Carb": w_carb_mean, "Clay": w_clay_mean,
                     "Sulfat": w_sulfat_mean, "Org": w_org_mean, "Sulfid": w_sulfid_mean}
        dict_mean_sorted = dict(sorted(dict_mean.items(), key=lambda item: item[1], reverse=True))

        n_minerals = len(dict_mean_sorted)
        condition = False
        while condition == False:
            helper_results = {}
            helper_total_amount = 100
            for index, (key, value) in enumerate(dict_mean_sorted.items()):
                values = helper_data[key]
                if index == 0:
                    val_amount = round(np.random.uniform(values[0], values[1]), 4)
                    helper_total_amount -= val_amount
                    helper_results[key] = val_amount
                else:
                    if index < n_minerals - 1:
                        if values[1] < helper_total_amount and values[1] != 0.0:
                            val_amount = round(np.random.uniform(values[0], values[1]), 4)
                        else:
                            val_amount = round(np.random.uniform(values[0], helper_total_amount), 4)

                        helper_total_amount -= val_amount
                        helper_results[key] = val_amount
                    else:
                        val_amount = round(helper_total_amount, 4)
                        helper_results[key] = val_amount

            helper_total_amount = 0
            for index, (key, value) in enumerate(helper_results.items()):
                helper_total_amount += value

            if round(helper_total_amount, 4) == 100.0:
                condition = True

        return helper_results

    def rock_builder(self):
        ## Initialization
        self.gui_variables["Entry"]["Number Datapoints"] = tk.IntVar()
        self.gui_variables["Entry"]["Number Datapoints"].set(100)
        self.gui_variables["Entry"]["Porosity Min"] = tk.IntVar()
        self.gui_variables["Entry"]["Porosity Min"].set(0)
        self.gui_variables["Entry"]["Porosity Max"] = tk.IntVar()
        self.gui_variables["Entry"]["Porosity Max"].set(10)
        self.gui_variables["Radiobutton"]["Mineralogy"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Mineralogy"].set(0)
        self.gui_variables["Radiobutton"]["Mineralogy Mode"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Mineralogy Mode"].set(0)
        self.last_mineralogy_mode = None
        self.gui_variables["Checkbox"]["Mineralogy"] = {}
        self.gui_variables["Entry"]["Mineralogy"] = {}
        self.gui_variables["Entry"]["Mineralogy"]["Minimum"] = {}
        self.gui_variables["Entry"]["Mineralogy"]["Maximum"] = {}
        self.gui_variables["Entry"]["Mineralogy"]["Clay Total"] = {}
        self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Minimum"] = tk.IntVar()
        self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Maximum"] = tk.IntVar()
        self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Minimum"].set(0)
        self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Maximum"].set(20)
        self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"] = {}
        self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Minimum"] = tk.IntVar()
        self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Maximum"] = tk.IntVar()
        self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Minimum"].set(80)
        self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Maximum"].set(100)
        self.gui_variables["Entry"]["Mineralogy"]["Ore Total"] = {}
        self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Minimum"] = tk.IntVar()
        self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Maximum"] = tk.IntVar()
        self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Minimum"].set(0)
        self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Maximum"].set(5)
        #
        self.gui_variables["Entry"]["Minimum"] = {}
        self.gui_variables["Entry"]["Maximum"] = {}
        self.gui_variables["Entry"]["Mean"] = {}
        self.gui_variables["Entry"]["Error"] = {}
        categories_short = ["rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi"]
        for category in categories_short:
            self.gui_variables["Entry"]["Minimum"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][category].set(0.0)
            self.gui_variables["Entry"]["Maximum"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][category].set(0.0)
            self.gui_variables["Entry"]["Mean"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][category].set(0.0)
            self.gui_variables["Entry"]["Error"][category] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][category].set(0.0)
        #
        self.rock_forming_minerals = [
            "Quartz", "Plagioclase", "Alkaline Feldspar", "Olivine", "Calcite", "Dolomite", "Muscovite", "Biotite",
            "Ca-Amphibole", "Na-Amphibole", "Mg-Fe-Pyroxene", "Ca-Pyroxene", "Al-Garnet", "Ca-Garnet",
            "Organic Matter"]
        self.ore_minerals = [
            "Acanthite", "Barite", "Boehmite", "Beryl", "Bornite", "Cassiterite", "Chalcocite", "Chalcopyrite",
            "Chromite", "Cinnabar", "Cobalite", "Coltan", "Galena", "Gold", "Hematite", "Ilmenite", "Magnetite",
            "Malachite", "Molybdenite", "Pentlandite", "Pyrolusite", "Scheelite", "Smithsonite", "Sperrylite",
            "Sphalerite", "Uraninite", "Wolframite", "Pollucite", "Pyrite", "Marmatite", "Au(III)-Oxide"]
        self.clay_minerals = [
            "Illite", "Kaolinite", "Montmorillonite", "Nontronite", "Saponite", "Chlorite", "Vermiculite"]
        #
        ## Labels
        lbl_title = SimpleElements(
            parent=self.parent, row_id=5, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Rock Builder", font_option="sans 14 bold", relief=tk.FLAT)
        lbl_mineralogy = SimpleElements(
            parent=self.parent, row_id=8, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Define Mineralogy", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_datapoints = SimpleElements(
            parent=self.parent, row_id=11, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Number of Datapoints", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_phi_min = SimpleElements(
            parent=self.parent, row_id=14, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Porosity Minimum", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_phi_max = SimpleElements(
            parent=self.parent, row_id=16, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Porosity Maximum", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Rockbuilder Static"]["Label"].extend(
            [lbl_title, lbl_mineralogy, lbl_datapoints, lbl_phi_min, lbl_phi_max])
        #
        ## Buttons
        btn_mineralogy = SimpleElements(
            parent=self.parent, row_id=8, column_id=16, n_rows=2, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Define Mineralogy", command=self.define_mineralogy)
        btn_simulation = SimpleElements(
            parent=self.parent, row_id=19, column_id=16, n_rows=2, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Run Simulation", command=self.run_simulation_rockbuilder)
        #
        self.gui_elements["Rockbuilder Static"]["Button"].extend([btn_mineralogy, btn_simulation])
        #
        ## Entries
        entr_datapoints = SimpleElements(
            parent=self.parent, row_id=11, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Number Datapoints"])
        entr_phi_min = SimpleElements(
            parent=self.parent, row_id=14, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Porosity Min"])
        entr_phi_max = SimpleElements(
            parent=self.parent, row_id=16, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Porosity Max"])
        #
        self.gui_elements["Rockbuilder Static"]["Entry"].extend([entr_datapoints, entr_phi_min, entr_phi_max])
        #
        ## TREEVIEWS
        list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
        list_width = list(75*np.ones(len(list_categories)))
        list_width = [int(item) for item in list_width]
        list_width[0] = 90
        list_width[-1] = 150
        #
        self.tv_petrology_results = SimpleElements(
            parent=self.parent, row_id=0, column_id=35, n_rows=30, n_columns=45,
            fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
            n_categories=len(list_categories), text_n=list_categories,
            width_n=list_width, individual=True)
        #
        scb_v = ttk.Scrollbar(self.parent, orient="vertical")
        scb_h = ttk.Scrollbar(self.parent, orient="horizontal")
        self.tv_petrology_results.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
        scb_v.config(command=self.tv_petrology_results.yview)
        scb_h.config(command=self.tv_petrology_results.xview)
        scb_v.grid(row=0, column=35 + 45, rowspan=30, columnspan=1, sticky="ns")
        scb_h.grid(row=30, column=35, rowspan=1, columnspan=45, sticky="ew")
    #
    def define_mineralogy(self):
        self.window_mineralogy = tk.Toplevel(self.parent)
        self.window_mineralogy.title("Mineralogy")
        self.window_mineralogy.geometry("1350x825")
        self.window_mineralogy.resizable(False, False)
        self.window_mineralogy["bg"] = self.colors_gebpy["Background"]
        #
        ## Cleaning
        categories = ["Frame", "Label", "Radiobutton", "Entry", "Checkbox"]
        priorities = ["Static", "Temporary"]
        for priority in priorities:
            for category in categories:
                if len(self.gui_elements_sub["Mineralogy"][priority][category]) > 0:
                    self.gui_elements_sub["Mineralogy"][priority][category].clear()
        #
        ## Geometry and Layout
        window_width = 1350
        window_heigth = 825
        row_min = 15
        n_rows = int(window_heigth / row_min)
        column_min = 15
        n_columns = int(window_width / column_min)
        #
        for x in range(n_columns):
            tk.Grid.columnconfigure(self.window_mineralogy, x, weight=1)
        for y in range(n_rows):
            tk.Grid.rowconfigure(self.window_mineralogy, y, weight=1)
        #
        # Rows
        for i in range(0, n_rows):
            self.window_mineralogy.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            self.window_mineralogy.grid_columnconfigure(i, minsize=column_min)
        #
        ## Frames
        frm_navigation = SimpleElements(
            parent=self.window_mineralogy, row_id=0, column_id=0, n_rows=n_rows, n_columns=19,
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Navigation"]).create_frame()
        frm_accent = SimpleElements(
            parent=self.window_mineralogy, row_id=0, column_id=20, n_rows=n_rows, n_columns=1,
            bg=self.colors_gebpy["Accent Blue"], fg=self.colors_gebpy["Accent Blue"]).create_frame()
        #
        self.gui_elements_sub["Mineralogy"]["Static"]["Frame"].extend([frm_navigation, frm_accent])
        #
        ## Labels
        lbl_title = SimpleElements(
            parent=self.window_mineralogy, row_id=0, column_id=0, n_rows=2, n_columns=19,
            bg=self.colors_gebpy["Accent Blue"], fg=self.colors_gebpy["Navigation"]).create_label(
            text="Mineral Class", font_option="sans 12 bold", relief=tk.FLAT)
        #
        self.gui_elements_sub["Mineralogy"]["Static"]["Label"].append(lbl_title)
        #
        ## Radiobuttons
        self.mineral_classes = [
            "Oxides", "Sulfides", "Sulfates", "Phosphates", "Phospides", "Carbonates", "Tectosilicates",
            "Nesosilicates", "Inosilicates", "Phyllosilicates", "Sorosilicates", "Cyclosilicates", "Miscellaneous"]
        self.mineral_classes.sort()
        # for index, mineral_class in enumerate(self.mineral_classes):
        #     rb_class = SimpleElements(
        #         parent=self.window_mineralogy, row_id=2*index + 2, column_id=0, n_rows=2, n_columns=12,
        #         bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
        #         text=mineral_class, var_rb=self.gui_variables["Radiobutton"]["Mineralogy"], value_rb=index,
        #         color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_mineral_class)
        #     #
        #     self.gui_elements_sub["Mineralogy"]["Static"]["Radiobutton"].append(rb_class)
        #
        self.mineral_classes_advanced = ["Rock-Forming Minerals", "Ore Minerals", "Clay Minerals"]
        self.mineral_classes_extended = [
            "Oxides", "Sulfides", "Sulfates", "Phosphates", "Phospides", "Carbonates", "Tectosilicates",
            "Nesosilicates", "Inosilicates", "Phyllosilicates", "Sorosilicates", "Cyclosilicates", "Miscellaneous"]
        self.mineral_classes_extended.sort()
        self.mineral_classes_extended.extend(self.mineral_classes_advanced)
        #
        # index_rockforming = 0
        # for index, mineral_class_adv in enumerate(self.mineral_classes_advanced, start=len(self.mineral_classes)):
        #     rb_class_adv = SimpleElements(
        #         parent=self.window_mineralogy, row_id=2*index + 3, column_id=0, n_rows=2,
        #         n_columns=12, bg=self.colors_gebpy["Navigation"],
        #         fg=self.colors_gebpy["Background"]).create_radiobutton(
        #         text=mineral_class_adv, var_rb=self.gui_variables["Radiobutton"]["Mineralogy"], value_rb=index,
        #         color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_mineral_class)
        #     #
        #     if mineral_class_adv == "Rock-Forming Minerals":
        #         index_rockforming = index
        #     #
        #     self.gui_elements_sub["Mineralogy"]["Static"]["Radiobutton"].append(rb_class_adv)
        # #
        # self.gui_variables["Radiobutton"]["Mineralogy"].set(index_rockforming)
        #
        modes = ["Simplified", "Complex"]
        for index, mode in enumerate(modes):
            rb_mode = SimpleElements(
                parent=self.window_mineralogy, row_id=2, column_id=9*index, n_rows=2, n_columns=9,
                bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                text=mode, var_rb=self.gui_variables["Radiobutton"]["Mineralogy Mode"], value_rb=index,
                color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_mineral_mode)
            #
            self.gui_elements_sub["Mineralogy"]["Static"]["Radiobutton"].append(rb_mode)
        #
        self.change_rb_mineral_mode()
        #self.change_rb_mineral_class()
    #
    def change_rb_mineral_mode(self):

        var_rb = self.gui_variables["Radiobutton"]["Mineralogy Mode"].get()
        #
        if var_rb != self.last_mineralogy_mode:
            ## Cleaning
            categories = ["Radiobutton"]
            for category in categories:
                for gui_item in self.gui_elements_sub["Mineralogy"]["Temporary"][category]:
                    gui_item.grid_remove()
                self.gui_elements_sub["Mineralogy"]["Temporary"][category].clear()
            #
            if var_rb == 0:
                index_rockforming = 0
                for index, mineral_class_adv in enumerate(self.mineral_classes_advanced, start=len(self.mineral_classes)):
                    rb_class_adv = SimpleElements(
                        parent=self.window_mineralogy, row_id=2*index + 5 - 2*len(self.mineral_classes), column_id=0,
                        n_rows=2, n_columns=12, bg=self.colors_gebpy["Navigation"],
                        fg=self.colors_gebpy["Background"]).create_radiobutton(
                        text=mineral_class_adv, var_rb=self.gui_variables["Radiobutton"]["Mineralogy"], value_rb=index,
                        color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_mineral_class)
                    #
                    if mineral_class_adv == "Rock-Forming Minerals":
                        index_rockforming = index
                    #
                    self.gui_elements_sub["Mineralogy"]["Temporary"]["Radiobutton"].append(rb_class_adv)
                #
                self.gui_variables["Radiobutton"]["Mineralogy"].set(index_rockforming)
                self.change_rb_mineral_class()
            else:
                for index, mineral_class in enumerate(self.mineral_classes):
                    rb_class = SimpleElements(
                        parent=self.window_mineralogy, row_id=2*index + 5, column_id=0, n_rows=2, n_columns=12,
                        bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                        text=mineral_class, var_rb=self.gui_variables["Radiobutton"]["Mineralogy"], value_rb=index,
                        color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_mineral_class)
                    #
                    self.gui_elements_sub["Mineralogy"]["Temporary"]["Radiobutton"].append(rb_class)
                #
                self.gui_variables["Radiobutton"]["Mineralogy"].set(0)
                self.change_rb_mineral_class()
        else:
            print("Nothing has changed!")
        #
        self.last_mineralogy_mode = var_rb
    #
    def change_rb_mineral_class(self):
        ## Cleaning
        categories = ["Label", "Entry", "Checkbox"]
        for category in categories:
            for gui_item in self.gui_elements_sub["Mineralogy"]["Temporary"][category]:
                gui_item.grid_remove()
            self.gui_elements_sub["Mineralogy"]["Temporary"][category].clear()
        #
        mineral_class_selected = self.mineral_classes_extended[self.gui_variables["Radiobutton"]["Mineralogy"].get()]
        #
        if mineral_class_selected == "Oxides":
            minerals_list = self.oxide_minerals
        elif mineral_class_selected == "Sulfides":
            minerals_list = self.sulfide_minerals
        elif mineral_class_selected == "Sulfates":
            minerals_list = self.sulfate_minerals
        elif mineral_class_selected == "Phosphates":
            minerals_list = self.phosphate_minerals
        elif mineral_class_selected == "Phospides":
            minerals_list = self.phospide_minerals
        elif mineral_class_selected == "Carbonates":
            minerals_list = self.carbonate_minerals
        elif mineral_class_selected == "Tectosilicates":
            minerals_list = self.tectosilicate_minerals
        elif mineral_class_selected == "Nesosilicates":
            minerals_list = self.nesosilicate_minerals
        elif mineral_class_selected == "Inosilicates":
            minerals_list = self.inosilicate_minerals
        elif mineral_class_selected == "Phyllosilicates":
            minerals_list = self.phyllosilicate_minerals
        elif mineral_class_selected == "Sorosilicates":
            minerals_list = self.sorosilicate_minerals
        elif mineral_class_selected == "Cyclosilicates":
            minerals_list = self.cyclosilicate_minerals
        elif mineral_class_selected == "Miscellaneous":
            minerals_list = self.miscellaneous_minerals
        elif mineral_class_selected == "Rock-Forming Minerals":
            minerals_list = [
                "Quartz", "Plagioclase", "Alkaline Feldspar", "Olivine", "Calcite", "Dolomite", "Muscovite", "Biotite",
                "Ca-Amphibole", "Na-Amphibole", "Mg-Fe-Pyroxene", "Ca-Pyroxene", "Al-Garnet", "Ca-Garnet",
                "Organic Matter"]
        elif mineral_class_selected == "Ore Minerals":
            minerals_list = ["Acanthite", "Barite", "Boehmite", "Beryl", "Bornite", "Cassiterite", "Chalcocite",
                             "Chalcopyrite", "Chromite", "Cinnabar", "Cobaltite", "Coltan", "Galena", "Gold", "Hematite",
                             "Ilmenite", "Magnetite", "Malachite", "Molybdenite", "Pentlandite", "Pyrolusite",
                             "Scheelite", "Smithsonite", "Sperrylite", "Sphalerite", "Uraninite", "Wolframite",
                             "Pollucite", "Pyrite", "Marmatite", "Au(III)-Oxide"]
        elif mineral_class_selected == "Clay Minerals":
            minerals_list = ["Illite", "Kaolinite", "Montmorillonite", "Nontronite", "Saponite", "Chlorite",
                             "Vermiculite"]
        #
        ## Labels
        minerals_list.sort()
        for index, mineral in enumerate(minerals_list):
            if index == 0 and len(self.gui_elements_sub["Mineralogy"]["Temporary"]["Label"]) == 0:
                lbl_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=22, n_rows=2, n_columns=12,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mineral", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=34, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Min", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=37, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Max", font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements_sub["Mineralogy"]["Temporary"]["Label"].extend(
                    [lbl_mineral, lbl_min, lbl_max])
            #
            elif index == 24:
                lbl_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=41, n_rows=2, n_columns=12,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mineral", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=53, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Min", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=56, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Max", font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements_sub["Mineralogy"]["Temporary"]["Label"].extend(
                    [lbl_mineral, lbl_min, lbl_max])
            #
            elif index == 48:
                lbl_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=60, n_rows=2, n_columns=12,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mineral", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=72, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Min", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=0, column_id=75, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Max", font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements_sub["Mineralogy"]["Temporary"]["Label"].extend(
                    [lbl_mineral, lbl_min, lbl_max])
            #
            if index < 24:
                ## Labels
                lbl_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2, column_id=22, n_rows=2, n_columns=10,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=mineral, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if mineral not in self.gui_variables["Checkbox"]["Mineralogy"]:
                    self.gui_variables["Checkbox"]["Mineralogy"][mineral] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Mineralogy"][mineral].set(0)
                #
                cb_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2, column_id=32, n_rows=2, n_columns=2,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Mineralogy"][mineral])
                #
                ## Entries
                if mineral not in self.gui_variables["Entry"]["Mineralogy"]["Minimum"]:
                    self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral] = tk.IntVar()
                    self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral].set(0)
                    self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral] = tk.IntVar()
                    self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral].set(100)
                #
                entr_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2, column_id=34, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral])
                entr_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2, column_id=37, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral])
            #
            elif 24 <= index < 48:
                ## Labels
                lbl_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 48, column_id=41, n_rows=2,
                    n_columns=10,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=mineral, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if mineral not in self.gui_variables["Checkbox"]["Mineralogy"]:
                    self.gui_variables["Checkbox"]["Mineralogy"][mineral] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Mineralogy"][mineral].set(0)
                cb_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 48, column_id=51, n_rows=2,
                    n_columns=2,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Mineralogy"][mineral])
                #
                ## Entries
                if mineral not in self.gui_variables["Entry"]["Mineralogy"]["Minimum"]:
                    self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral] = tk.IntVar()
                    self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral].set(0)
                    self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral] = tk.IntVar()
                    self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral].set(100)
                #
                entr_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 48, column_id=53, n_rows=2,
                    n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral])
                entr_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 48, column_id=56, n_rows=2,
                    n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral])
            #
            elif 48 <= index < 72:
                ## Labels
                lbl_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=2*index + 2 - 72, column_id=60, n_rows=2,
                    n_columns=10,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=mineral, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if mineral not in self.gui_variables["Checkbox"]["Mineralogy"]:
                    self.gui_variables["Checkbox"]["Mineralogy"][mineral] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Mineralogy"][mineral].set(0)
                cb_mineral = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 72, column_id=70, n_rows=2,
                    n_columns=2,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Mineralogy"][mineral])
                #
                ## Entries
                if mineral not in self.gui_variables["Entry"]["Mineralogy"]["Minimum"]:
                    self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral] = tk.IntVar()
                    self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral].set(0)
                    self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral] = tk.IntVar()
                    self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral].set(100)
                #
                entr_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 72, column_id=72, n_rows=2,
                    n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral])
                entr_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=2 * index + 2 - 72, column_id=75, n_rows=2,
                    n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral])
            #
            self.gui_elements_sub["Mineralogy"]["Temporary"]["Label"].append(lbl_mineral)
            self.gui_elements_sub["Mineralogy"]["Temporary"]["Checkbox"].append(cb_mineral)
            self.gui_elements_sub["Mineralogy"]["Temporary"]["Entry"].extend([entr_min, entr_max])
        #
        if mineral_class_selected in ["Rock-Forming Minerals", "Ore Minerals", "Clay Minerals"]:
            ## Entries
            if len(self.gui_elements_sub["Mineralogy"]["Static"]["Entry"]) == 0:
                entr_fraction_total_rf_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=5, column_id=12, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Minimum"])
                entr_fraction_total_rf_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=5, column_id=15, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Maximum"])
                #
                entr_fraction_total_ore_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=7, column_id=12, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Minimum"])
                entr_fraction_total_ore_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=7, column_id=15, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Maximum"])
                #
                entr_fraction_total_clay_min = SimpleElements(
                    parent=self.window_mineralogy, row_id=9, column_id=12, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Minimum"])
                entr_fraction_total_clay_max = SimpleElements(
                    parent=self.window_mineralogy, row_id=9, column_id=15, n_rows=2, n_columns=3,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Maximum"])
                #
                self.gui_elements_sub["Mineralogy"]["Static"]["Entry"].extend(
                    [entr_fraction_total_rf_min, entr_fraction_total_rf_max, entr_fraction_total_ore_min,
                     entr_fraction_total_ore_max, entr_fraction_total_clay_min, entr_fraction_total_clay_max])
            else:
                for category in ["Entry"]:
                    for gui_element in self.gui_elements_sub["Mineralogy"]["Static"][category]:
                        gui_element.grid()
        #
        else:
            for category in ["Entry"]:
                for gui_element in self.gui_elements_sub["Mineralogy"]["Static"][category]:
                    gui_element.grid_remove()
    #
    def run_simulation_rockbuilder(self):
        try:
            self.canvas_scatter.get_tk_widget().grid_forget()
        except:
            pass
        #
        ## Initialization
        self.gui_variables["Radiobutton"]["Analysis Mode"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Rock"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Rock"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Mineral"].set(0)
        self.gui_variables["Radiobutton"]["Diagram Type Element"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type Element"].set(0)
        self.gui_variables["Radiobutton"]["Concentration Type"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Concentration Type"].set(0)
        self.gui_variables["Radiobutton"]["Amount Mineral"] = tk.StringVar()
        self.gui_variables["Radiobutton"]["Amount Mineral"].set("Select Mineral")
        self.gui_variables["Radiobutton"]["Amount Element"] = tk.StringVar()
        self.gui_variables["Radiobutton"]["Amount Element"].set("Select Element")
        self.gui_variables["Radiobutton"]["Amount Oxide"] = tk.StringVar()
        self.gui_variables["Radiobutton"]["Amount Oxide"].set("Select Oxide")
        #
        self.last_rb_analysis_rock.set(42)
        #
        self.selected_minerals = {}
        n_digits = 8
        #
        ## Preparation Simulation
        var_name = "Custom Rock"
        data_water = Water.water("")
        #
        n_datapoints = self.gui_variables["Entry"]["Number Datapoints"].get()
        phi_min = self.gui_variables["Entry"]["Porosity Min"].get()
        phi_max = self.gui_variables["Entry"]["Porosity Max"].get()
        #
        rock_forming_min = self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Minimum"].get()
        rock_forming_max = self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Maximum"].get()
        n_rock_forming = 0
        mean_rock_forming = (rock_forming_min + rock_forming_max)/2
        ore_min = self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Minimum"].get()
        ore_max = self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Maximum"].get()
        n_ore = 0
        mean_ore = (ore_min + ore_max)/2
        clay_min = self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Minimum"].get()
        clay_max = self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Maximum"].get()
        n_clay = 0
        mean_clay = (clay_min + clay_max)/2
        #
        selected_minerals = {}
        #
        self.calculate_fractions()
        for mineral, value in self.gui_variables["Checkbox"]["Mineralogy"].items():
            if value.get() == 1:
                value_min = self.selected_minerals[mineral]["Min"]
                value_max = self.selected_minerals[mineral]["Max"]
                #
                selected_minerals[mineral] = {"Min": value_min, "Max": value_max}
                #
                if mineral in self.rock_forming_minerals:
                    n_rock_forming += 1
                elif mineral in self.ore_minerals:
                    n_ore += 1
                elif mineral in self.clay_minerals:
                    n_clay += 1
        #
        if n_rock_forming > 0:
            if n_ore > 0 and n_clay == 0:
                if mean_rock_forming >= mean_ore:
                    mean = (rock_forming_min + rock_forming_max)/2
                    sigma = (mean - rock_forming_min)/3
                    val_fraction = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)
                    #
                    fraction_rock_forming = val_fraction
                    fraction_ore = 100*np.ones(n_datapoints) - fraction_rock_forming
                    fraction_clay = np.zeros(n_datapoints)
                else:
                    mean = (ore_min + ore_max)/2
                    sigma = (mean - ore_min)/3
                    val_fraction = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)
                    #
                    fraction_ore = val_fraction
                    fraction_rock_forming = 100*np.ones(n_datapoints) - fraction_ore
                    fraction_clay = np.zeros(n_datapoints)
            elif n_ore == 0 and n_clay > 0:
                if mean_rock_forming >= mean_clay:
                    mean = (rock_forming_min + rock_forming_max)/2
                    sigma = (mean - rock_forming_min)/3
                    val_fraction = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)
                    #
                    fraction_rock_forming = val_fraction
                    fraction_clay = 100*np.ones(n_datapoints) - fraction_rock_forming
                    fraction_ore = np.zeros(n_datapoints)
                else:
                    mean = (clay_min + clay_max)/2
                    sigma = (mean - clay_min)/3
                    val_fraction = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)
                    #
                    fraction_clay = val_fraction
                    fraction_rock_forming = 100*np.ones(n_datapoints) - fraction_clay
                    fraction_ore = np.zeros(n_datapoints)
            elif n_ore > 0 and n_clay > 0:
                if mean_rock_forming >= mean_ore and mean_rock_forming >= mean_clay:
                    if mean_ore >= mean_clay:   # clay < ore < rock forming
                        mean_rf = (rock_forming_min + rock_forming_max)/2
                        sigma_rf = (mean_rf - rock_forming_min)/3
                        val_fraction_rf = np.random.normal(loc=mean_rf, scale=sigma_rf, size=n_datapoints)
                        #
                        mean_o = (ore_min + ore_max)/2
                        sigma_o = (mean_o - ore_min)/3
                        val_fraction_o = np.random.normal(loc=mean_o, scale=sigma_o, size=n_datapoints)
                        #
                        fraction_rock_forming = val_fraction_rf
                        fraction_ore = val_fraction_o
                        fraction_clay = 100*np.ones(n_datapoints) - fraction_rock_forming - fraction_ore
                    else:                       # ore < clay < rock forming
                        mean_rf = (rock_forming_min + rock_forming_max)/2
                        sigma_rf = (mean_rf - rock_forming_min)/3
                        val_fraction_rf = np.random.normal(loc=mean_rf, scale=sigma_rf, size=n_datapoints)
                        #
                        mean_c = (clay_min + clay_max)/2
                        sigma_c = (mean_c - clay_min)/3
                        val_fraction_c = np.random.normal(loc=mean_c, scale=sigma_c, size=n_datapoints)
                        #
                        fraction_rock_forming = val_fraction_rf
                        fraction_clay = val_fraction_c
                        fraction_ore = 100*np.ones(n_datapoints) - fraction_rock_forming - fraction_clay
                elif mean_rock_forming >= mean_ore and mean_rock_forming < mean_clay:   # ore < rock forming < clay
                    mean_c = (clay_min + clay_max)/2
                    sigma_c = (mean_c - clay_min)/3
                    val_fraction_c = np.random.normal(loc=mean_c, scale=sigma_c, size=n_datapoints)
                    #
                    mean_rf = (rock_forming_min + rock_forming_max)/2
                    sigma_rf = (mean_rf - rock_forming_min)/3
                    val_fraction_rf = np.random.normal(loc=mean_rf, scale=sigma_rf, size=n_datapoints)
                    #
                    fraction_clay = val_fraction_c
                    fraction_rock_forming = val_fraction_rf
                    fraction_ore = 100*np.ones(n_datapoints) - fraction_clay - fraction_rock_forming
                elif mean_rock_forming < mean_ore and mean_rock_forming >= mean_clay:   # clay < rock forming < ore
                    mean_o = (ore_min + ore_max)/2
                    sigma_o = (mean_o - ore_min)/3
                    val_fraction_o = np.random.normal(loc=mean_o, scale=sigma_o, size=n_datapoints)
                    #
                    mean_rf = (rock_forming_min + rock_forming_max)/2
                    sigma_rf = (mean_rf - rock_forming_min)/3
                    val_fraction_rf = np.random.normal(loc=mean_rf, scale=sigma_rf, size=n_datapoints)
                    #
                    fraction_ore = val_fraction_o
                    fraction_rock_forming = val_fraction_rf
                    fraction_clay = 100*np.ones(n_datapoints) - fraction_ore - fraction_rock_forming
                elif mean_rock_forming < mean_ore and mean_rock_forming < mean_clay:
                    if mean_ore >= mean_clay:   # rock forming < clay < ore
                        mean_o = (ore_min + ore_max)/2
                        sigma_o = (mean_o - ore_min)/3
                        val_fraction_o = np.random.normal(loc=mean_o, scale=sigma_o, size=n_datapoints)
                        #
                        mean_c = (clay_min + clay_max)/2
                        sigma_c = (mean_c - clay_min)/3
                        val_fraction_c = np.random.normal(loc=mean_c, scale=sigma_c, size=n_datapoints)
                        #
                        fraction_ore = val_fraction_o
                        fraction_clay = val_fraction_c
                        fraction_rock_forming = 100*np.ones(n_datapoints) - fraction_ore - fraction_clay
                    else:                       # rock forming < ore < clay
                        mean_c = (clay_min + clay_max)/2
                        sigma_c = (mean_c - clay_min)/3
                        val_fraction_c = np.random.normal(loc=mean_c, scale=sigma_c, size=n_datapoints)
                        #
                        mean_o = (ore_min + ore_max)/2
                        sigma_o = (mean_o - ore_min)/3
                        val_fraction_o = np.random.normal(loc=mean_o, scale=sigma_o, size=n_datapoints)
                        #
                        fraction_clay = val_fraction_c
                        fraction_ore = val_fraction_o
                        fraction_rock_forming = 100*np.ones(n_datapoints) - fraction_clay - fraction_ore
            else:
                fraction_rock_forming = 100*np.ones(n_datapoints)
                fraction_ore = np.zeros(n_datapoints)
                fraction_clay = np.zeros(n_datapoints)
        else:
            if n_ore > 0 and n_clay == 0:
                fraction_rock_forming = np.zeros(n_datapoints)
                fraction_ore = 100*np.ones(n_datapoints)
                fraction_clay = np.zeros(n_datapoints)
            elif n_ore == 0 and n_clay > 0:
                fraction_rock_forming = np.zeros(n_datapoints)
                fraction_ore = np.zeros(n_datapoints)
                fraction_clay = 100*np.ones(n_datapoints)
            else:
                if mean_ore >= mean_clay:
                    mean_o = (ore_min + ore_max)/2
                    sigma_o = (mean_o - ore_min)/3
                    val_fraction_o = np.random.normal(loc=mean_o, scale=sigma_o, size=n_datapoints)
                    #
                    val_fraction_c = 100*np.ones(n_datapoints) - val_fraction_o
                    #
                    fraction_ore = val_fraction_o
                    fraction_clay = val_fraction_c
                    fraction_rock_forming = np.zeros(n_datapoints)
                else:
                    mean_c = (clay_min + clay_max)/2
                    sigma_c = (mean_c - clay_min)/3
                    val_fraction_c = np.random.normal(loc=mean_c, scale=sigma_c, size=n_datapoints)
                    #
                    val_fraction_o = 100*np.ones(n_datapoints) - val_fraction_c
                    #
                    fraction_ore = val_fraction_o
                    fraction_clay = val_fraction_c
                    fraction_rock_forming = np.zeros(n_datapoints)

        fractions = {
            "Rock Forming": np.around(fraction_rock_forming/100, n_digits),
            "Ore": np.around(fraction_ore/100, n_digits),
            "Clay": np.around(fraction_clay/100, n_digits)}

        mineral_list = list(selected_minerals.keys())
        mineral_list.sort()
        pure_minerals = ["Qz"]
        pure_minerals.sort()
        #
        ## Simulation
        self.data_rock = {}
        categories = ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi", "fluid",
                      "mineralogy", "chemistry"]
        for category in categories:
            if category not in ["mineralogy", "chemistry"]:
                self.data_rock[category] = []
            else:
                self.data_rock[category] = {}
        #
        mineral_data = {}
        mineral_amounts = {}
        mineral_names = {}
        for mineral in mineral_list:
            if mineral in self.oxide_minerals:
                data_mineral = Oxides(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.sulfide_minerals:
                data_mineral = Sulfides(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.sulfate_minerals:
                data_mineral = Sulfates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.carbonate_minerals:
                data_mineral = Carbonates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.phospide_minerals:
                data_mineral = Phospides(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.phosphate_minerals:
                data_mineral = Phosphates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.miscellaneous_minerals:
                data_mineral = Organics(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.tectosilicate_minerals:
                data_mineral = Tectosilicates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.sorosilicate_minerals:
                data_mineral = Sorosilicates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.nesosilicate_minerals:
                data_mineral = Nesosilicates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.inosilicate_minerals:
                data_mineral = Inosilicates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.cyclosilicate_minerals:
                data_mineral = Cyclosilicates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)
            elif mineral in self.phyllosilicate_minerals:
                data_mineral = Phyllosilicates(
                    mineral=mineral, data_type=True, traces_list=[]).generate_dataset(number=n_datapoints)

            val_min = selected_minerals[mineral]["Min"]
            val_max = selected_minerals[mineral]["Max"]
            mean = (val_min + val_max)/2
            sigma = (mean - val_min)/3

            if mineral in self.rock_forming_minerals:
                fraction_factor = self.minerals_helper["Rock Forming"][mineral]["Fraction"]/100

                val_min = self.minerals_helper["Rock Forming"][mineral]["Min"]
                val_max = self.minerals_helper["Rock Forming"][mineral]["Max"]
                mean = (val_min + val_max)/2
                sigma = (mean - val_min)/3

            elif mineral in self.ore_minerals:
                fraction_factor = self.minerals_helper["Ore"][mineral]["Fraction"]/100

            elif mineral in self.clay_minerals:
                fraction_factor = self.minerals_helper["Clay"][mineral]["Fraction"]/100

            amounts_raw = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)

            amounts_ppm = np.array([int(fraction_factor*value) for index, value in enumerate(amounts_raw)])
            amounts_ppm[amounts_ppm < 0] = 0
            amounts_ppm[amounts_ppm == 0] = 1

            mineral_amounts[data_mineral["mineral"]] = list(amounts_ppm)
            mineral_data[data_mineral["mineral"]] = data_mineral
            mineral_names[mineral] = data_mineral["mineral"]

        data_amounts_minerals = self.calculate_mineral_fractions(
            var_minerals=list(mineral_data.keys()), var_data=self.minerals_helper, var_n=n_datapoints)

        for mineral, mineral_fractions in data_amounts_minerals.items():
            mineral_amounts[mineral_names[mineral]] = mineral_fractions

        n = 0
        while n < n_datapoints:
            porosity = round(rd.uniform(phi_min, phi_max)/100, 5)
            rho_solid = 0
            bulk_modulus = 0
            shear_modulus = 0
            gamma_ray = 0
            photoelectricity = 0
            w_minerals = {}
            w_elements = {}
            elements_list = []
            phi_list = []
            phi_minerals = {}
            velocity_solid = {"vP": 0, "vS": 0}
            amounts_minerals = {"mass": {}, "volume": {}}
            sum_amounts_w = 0
            for mineral, dataset in mineral_data.items():
                if mineral not in self.data_rock["mineralogy"]:
                    self.data_rock["mineralogy"][mineral] = []

                if mineral not in amounts_minerals["mass"]:
                    amounts_minerals["mass"][mineral] = []
                    amounts_minerals["volume"][mineral] = []
                mineral_amount = round(mineral_amounts[mineral][n]*10**(-2), 6)
                phi_minerals[mineral] = mineral_amount
                amounts_minerals["volume"][mineral].append(mineral_amount)
                phi_list.append(mineral_amount)
                if mineral in pure_minerals:
                    sum_amounts_w += mineral_amount*dataset["rho"][0]
                    velocity_solid["vP"] += mineral_amount*dataset["vP"][0]
                    velocity_solid["vS"] += mineral_amount*dataset["vS"][0]
                    rho_solid += mineral_amount*dataset["rho"][0]
                    bulk_modulus += mineral_amount*dataset["K"][0]
                    shear_modulus += mineral_amount*dataset["G"][0]
                    gamma_ray += mineral_amount*dataset["GR"][0]
                    photoelectricity += mineral_amount*dataset["PE"][0]
                else:
                    sum_amounts_w += mineral_amount*dataset["rho"][n]
                    velocity_solid["vP"] += mineral_amount*dataset["vP"][n]
                    velocity_solid["vS"] += mineral_amount*dataset["vS"][n]
                    rho_solid += mineral_amount*dataset["rho"][n]
                    bulk_modulus += mineral_amount*dataset["K"][n]
                    shear_modulus += mineral_amount*dataset["G"][n]
                    gamma_ray += mineral_amount*dataset["GR"][n]
                    photoelectricity += mineral_amount*dataset["PE"][n]
                for element, value in dataset["chemistry"].items():
                    if element not in elements_list:
                        elements_list.append(element)
                        w_elements[element] = 0.0

                    if element not in self.data_rock["chemistry"]:
                        self.data_rock["chemistry"][element] = []

            for mineral, dataset in mineral_data.items():
                if mineral in pure_minerals:
                    amount_w = round((mineral_amounts[mineral][n]*10**(-2)*dataset["rho"][0])/sum_amounts_w, 6)
                else:
                    amount_w = round((mineral_amounts[mineral][n]*10**(-2)*dataset["rho"][n])/sum_amounts_w, 6)
                amounts_minerals["mass"][mineral].append(amount_w)

            ## Density
            rho_solid = round(rho_solid, 3)
            rho = round((1 - porosity)*rho_solid + porosity*data_water[2], 3)
            ## Seismic Velocities
            vP_factor = 4/3
            vP = round(velocity_solid["vP"]*(1 - vP_factor*porosity), 3)
            vS_factor = 3/2
            vS = round(velocity_solid["vS"]*(1 - vS_factor*porosity), 3)
            vPvS = round(vP/vS, 6)
            ## Elastic Parameters
            bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio = SeismicVelocities(
                rho_solid=None, rho_fluid=None).calculate_elastic_properties(
                rho=rho, vP=vP, vS=vS)
            ## Gamma Ray
            gamma_ray = round(gamma_ray, 3)
            ## Photoelectricity
            photoelectricity = round(photoelectricity, 3)

            for mineral, dataset in mineral_amounts.items():
                if mineral in pure_minerals:
                    mineral_amount = round(
                        (mineral_amounts[mineral][n]*mineral_data[mineral]["rho"][0])/(100*rho_solid), n_digits)
                else:
                    mineral_amount = round(
                        (mineral_amounts[mineral][n]*mineral_data[mineral]["rho"][n])/(100*rho_solid), n_digits)
                w_minerals[mineral] = mineral_amount
                self.data_rock["mineralogy"][mineral].append(mineral_amount)

            old_index = elements_list.index("O")
            elements_list += [elements_list.pop(old_index)]
            w_elements_total = 0.0
            for element in elements_list:
                if element != "O":
                    for mineral, w_mineral in w_minerals.items():
                        if element in mineral_data[mineral]["chemistry"]:
                            if mineral in pure_minerals:
                                value = round(w_mineral*mineral_data[mineral]["chemistry"][element][0], n_digits)
                            else:
                                value = round(w_mineral*mineral_data[mineral]["chemistry"][element][n], n_digits)
                            w_elements[element] += value
                            w_elements_total += value

                            w_elements[element] = round(w_elements[element], n_digits)
                elif element == "O":
                    w_elements[element] += round(1 - w_elements_total, n_digits)

                    w_elements[element] = round(w_elements[element], n_digits)

            for element, value in w_elements.items():
                self.data_rock["chemistry"][element].append(value)

            ## Results
            self.data_rock["rho_s"].append(rho_solid)
            self.data_rock["rho"].append(rho)
            self.data_rock["vP"].append(vP)
            self.data_rock["vS"].append(vS)
            self.data_rock["vP/vS"].append(vPvS)
            self.data_rock["K"].append(bulk_modulus)
            self.data_rock["G"].append(shear_modulus)
            self.data_rock["E"].append(youngs_modulus)
            self.data_rock["nu"].append(poisson_ratio)
            self.data_rock["GR"].append(gamma_ray)
            self.data_rock["PE"].append(photoelectricity)
            self.data_rock["phi"].append(porosity*100)
            #
            self.list_elements_rock = list(w_elements.keys())
            self.list_minerals_rock = list(w_minerals.keys())
            #
            n += 1

        ## Labels
        lbl_analysis = SimpleElements(
            parent=self.parent, row_id=22, column_id=0, n_rows=6, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Analysis Mode", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Rockbuilder Static"]["Label"].append(lbl_analysis)
        #
        ## Button
        name_rock = var_name
        btn_export = SimpleElements(
            parent=self.parent, row_id=29, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Export Data", command=lambda var_dataset=self.data_rock, var_name=name_rock:
            self.export_rock_data(var_dataset, var_name))
        #
        self.gui_elements["Rockbuilder Static"]["Button"].append(btn_export)
        #
        ## Radiobuttons
        rb_geophysics = SimpleElements(
            parent=self.parent, row_id=22, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Physical Data", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=0,
            color_bg=self.colors_gebpy["Background"], command=self.change_rb_analysis_rocks)
        rb_geochemistry = SimpleElements(
            parent=self.parent, row_id=24, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Composition", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis_rocks)
        rb_geochemistry2 = SimpleElements(
            parent=self.parent, row_id=26, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Element Composition", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=2,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis_rocks)
        #
        self.gui_elements["Rockbuilder Static"]["Radiobutton"].extend(
            [rb_geophysics, rb_geochemistry, rb_geochemistry2])
        #
        for key, gui_element in self.gui_elements["Rockbuilder Temporary"].items():
            if key not in ["Canvas", "Button"]:
                if type(gui_element) == list:
                    for gui_item in gui_element:
                        gui_item.grid_remove()
            elif key == "Canvas":
                for key_2, gui_item in gui_element.items():
                    gui_item.get_tk_widget().grid_remove()
            gui_element.clear()
        #
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.change_rb_analysis_rocks()
        #
        for mineral in self.list_minerals_rock:
            self.gui_variables["Entry"]["Minimum"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][mineral].set(0.0)
            self.gui_variables["Entry"]["Maximum"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][mineral].set(0.0)
            self.gui_variables["Entry"]["Mean"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][mineral].set(0.0)
            self.gui_variables["Entry"]["Error"][mineral] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][mineral].set(0.0)
        #
        for element in self.list_elements_rock:
            self.gui_variables["Entry"]["Minimum"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Minimum"][element].set(0.0)
            self.gui_variables["Entry"]["Maximum"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Maximum"][element].set(0.0)
            self.gui_variables["Entry"]["Mean"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Mean"][element].set(0.0)
            self.gui_variables["Entry"]["Error"][element] = tk.StringVar()
            self.gui_variables["Entry"]["Error"][element].set(0.0)
        #
        if len(self.tv_petrology_results.get_children()) > 0:
            for item in self.tv_petrology_results.get_children():
                self.tv_petrology_results.delete(item)
        #
        for property, dataset in self.data_rock.items():
            if property not in ["fluid", "mineralogy", "chemistry"]:
                entries = [property]
                #
                n_digits = 2
                #
                var_entr_min = round(min(dataset), n_digits)
                var_entr_max = round(max(dataset), n_digits)
                var_entr_mean = round(np.mean(dataset), n_digits)
                var_entr_error = round(np.std(dataset, ddof=1), n_digits)
                #
                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
                #
                self.tv_petrology_results.insert("", tk.END, values=entries)
        #
        entries = ["-", "-", "-", "-", "-"]
        self.tv_petrology_results.insert("", tk.END, values=entries)
        #
        for mineral in np.sort(self.list_minerals_rock):
            dataset = self.data_rock["mineralogy"][mineral]
            entries = [str(mineral)+str(" (%)")]
            #
            n_digits = 2
            var_factor = 100
            #
            var_entr_min = round(var_factor*min(dataset), n_digits)
            var_entr_max = round(var_factor*max(dataset), n_digits)
            var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
            var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)
            #
            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            #
            self.tv_petrology_results.insert("", tk.END, values=entries)
        #
        entries = ["-", "-", "-", "-", "-"]
        self.tv_petrology_results.insert("", tk.END, values=entries)
        #
        for element in np.sort(self.list_elements_rock):
            dataset = self.data_rock["chemistry"][element]
            entries = [str(element)+str(" (%)")]
            #
            n_digits = 2
            var_factor = 100
            #
            var_entr_min = round(var_factor*min(dataset), n_digits)
            var_entr_max = round(var_factor*max(dataset), n_digits)
            var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
            var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)
            #
            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            #
            self.tv_petrology_results.insert("", tk.END, values=entries)

    def improve_xlimits_rock(self, delta_x, x_min, x_max):
        if delta_x < 1:
            n_digits = 3
        elif 1 <= delta_x < 5:
            n_digits = 2
        elif delta_x >= 5:
            n_digits = 0

        x_min = round(x_min - 0.1*delta_x, n_digits)
        x_max = round(x_max + 0.1*delta_x, n_digits)

        if delta_x > 100:
            delta = 10
            x_min = self.myround(x=x_min - delta, prec=0, base=3)
            x_max = self.myround(x=x_max + delta, prec=0, base=3)
        elif 100 >= delta_x > 25:
            delta = 3
            x_min = self.myround(x=x_min - delta, prec=0, base=3)
            x_max = self.myround(x=x_max + delta, prec=0, base=3)
        elif 25 >= delta_x > 10:
            delta = 2
            x_min = self.myround(x=x_min - delta, prec=0, base=3)
            x_max = self.myround(x=x_max + delta, prec=0, base=3)
        elif 10 >= delta_x > 3:
            delta = 1.5
            x_min = self.myround(x=x_min - delta, prec=0, base=3)
            x_max = self.myround(x=x_max + delta, prec=0, base=3)
        elif 3 >= delta_x > 1:
            delta = 0.9
            x_min = self.myround(x=x_min - delta, prec=1, base=1.5)
            x_max = self.myround(x=x_max + delta, prec=1, base=1.5)
        else:
            delta = 0.015
            x_min = self.myround(x=x_min - delta, prec=3, base=0.003)
            x_max = self.myround(x=x_max + delta, prec=3, base=0.003)

        return x_min, x_max

    def improve_ylimits_rock(self, y_min, y_max):
        if (y_max - y_min) > 50:
            delta = 9
            y_min = self.myround(x=y_min - delta, base=3)
            y_max = self.myround(x=y_max + delta, base=3)
        elif 50 >= (y_max - y_min) > 25:
            delta = 4.5
            y_min = self.myround(x=y_min - delta, base=3)
            y_max = self.myround(x=y_max + delta, base=3)
        elif 25 >= (y_max - y_min) > 10:
            delta = 3
            y_min = self.myround(x=y_min - delta, base=3)
            y_max = self.myround(x=y_max + delta, base=3)
        else:
            if y_min > 5:
                delta = 3
                y_min = self.myround(x=y_min - delta, base=1.5)
                y_max = self.myround(x=y_max + delta, base=1.5)
            elif 5 >= (y_max - y_min) > 1:
                delta = 1.5
                y_min = self.myround(x=y_min - delta, base=0.3)
                y_max = self.myround(x=y_max + delta, base=0.3)
            else:
                delta = 0.015
                y_min = self.myround(x=y_min - delta, base=0.003)
                y_max = self.myround(x=y_max + delta, base=0.003)

        if y_min < 0:
            y_min = 0

        return y_min, y_max

    def calculate_delta_xy(self, dataset_x, dataset_y):
        x_min = min(dataset_x)
        x_max = max(dataset_x)
        delta_x = round(x_max - x_min, 4)
        y_min = min(dataset_y)
        y_max = max(dataset_y)
        delta_y = round(y_max - y_min, 4)

        return x_min, x_max, delta_x, y_min, y_max, delta_y

    def change_rb_diagram_rocks(self):  # RB DIAGRAM ROCKS
        if self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 0:   # ROCK PHYSICS
            var_rb_diagram = self.gui_variables["Radiobutton"]["Diagram Type Rock"].get()
            if var_rb_diagram == 0:   # HISTOGRAM
                if self.last_rb_setting["Petrology"]["General Mode"].get() not in [1, 2]:
                    if self.last_rb_setting["Petrology"]["Rock Physics"].get() != 0:
                        if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                          "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                        ## Diagram
                        if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                            fig_petrology = Figure(
                                figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                            ax_rp_histo = fig_petrology.subplots(nrows=3, ncols=3)
                            #
                            categories = [["rho", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "nu"]]
                            labels = [["$\\varrho$ (kg/m\u00B3)", "GR (API)", "PE (barns/e\u207B)"],
                                      ["vP (m/s)", "vS (m/s)", "vP/vS (1)"], ["K (GPa)", "G (GPa)", "nu (1)"]]
                            #
                            for i, subcategories in enumerate(categories):
                                for j, key in enumerate(subcategories):
                                    y, x, _ = ax_rp_histo[i][j].hist(
                                        x=self.data_rock[key], color=self.colors_gebpy["Option"], edgecolor="black",
                                        bins=12)

                                    x_min, x_max, delta_x, y_min, y_max, delta_y = self.calculate_delta_xy(
                                        dataset_x=x, dataset_y=y)
                                    y_min = 0
                                    y_max = round(1.05*max(y), 2)

                                    if delta_x < 1:
                                        n_digits = 3
                                    elif 1 <= delta_x < 5:
                                        n_digits = 2
                                    elif delta_x >= 5:
                                        n_digits = 0

                                    x_min, x_max = self.improve_xlimits_rock(delta_x=delta_x, x_min=x_min, x_max=x_max)
                                    y_min, y_max = self.improve_ylimits_rock(y_min=y_min, y_max=y_max)

                                    if key != "nu":
                                        if x_min < 0:
                                            x_min = 0

                                    ax_rp_histo[i][j].set_xlim(left=x_min, right=x_max)
                                    ax_rp_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                                    ax_rp_histo[i][j].set_xticks(np.around(
                                        np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits))
                                    ax_rp_histo[i][j].set_yticks(np.around(
                                        np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 1))
                                    ax_rp_histo[i][j].xaxis.set_tick_params(labelsize=8)
                                    ax_rp_histo[i][j].yaxis.set_tick_params(labelsize=8)
                                    ax_rp_histo[i][j].set_xlabel(labels[i][j], fontsize=8)
                                    ax_rp_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=8)

                                    ax_rp_histo[i][j].grid(which="major", axis="both", linestyle="-")
                                    ax_rp_histo[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                    ax_rp_histo[i][j].minorticks_on()
                                    ax_rp_histo[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rp_histo[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rp_histo[i][j].set_axisbelow(True)

                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                                canvas_petrology.get_tk_widget().grid(
                                    row=0, column=81, rowspan=int(self.n_rows - 3), columnspan=int(self.n_columns - 81),
                                    sticky="nesw")
                                #
                                self.gui_elements["Temporary"]["Axis"]["Rock Physics Histogram"] = ax_rp_histo
                                self.gui_elements["Temporary"]["Figure"]["Petrology"] = fig_petrology
                                self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                            #
                        else:
                            ## Cleaning
                            for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Physics Histogram"]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("on")
                                    gui_axis.set_visible(True)
                            #
                            categories = ["Rock Physics Scatter", "Rock Chemistry Histogram", "Rock Chemistry Scatter",
                                          "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                    else:
                        pass
                else:
                    ## RECONSTRUCTION
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Physics Histogram"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)
                    #
                    categories = ["Rock Physics Scatter", "Rock Chemistry Histogram", "Rock Chemistry Scatter",
                                  "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                #
            elif var_rb_diagram == 1:  # SCATTER
                if self.last_rb_setting["Petrology"]["General Mode"].get() not in [1, 2]:
                    if self.last_rb_setting["Petrology"]["Rock Physics"].get() != 1:
                        if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                          "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                        if "Rock Physics Scatter" not in self.gui_elements["Temporary"]["Axis"]:
                            if "Petrology" not in self.gui_elements["Temporary"]["Figure"]:
                                fig_petrology = Figure(
                                    figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                            else:
                                fig_petrology = self.gui_elements["Temporary"]["Figure"]["Petrology"]
                            #
                            ax_rp_scatter = fig_petrology.subplots(nrows=3, ncols=3)
                            #
                            categories = [["phi", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "nu"]]
                            labels = [["$\\varphi$ (%)", "GR (API)", "PE (barns/e\u207B)"],
                                      ["vP (m/s)", "vS (m/s)", "vP/vS (1)"], ["K (GPa)", "G (GPa)", "nu (1)"]]
                            #
                            x_key = "rho"
                            dataset_x = self.data_rock[x_key]
                            #
                            for i, subcategories in enumerate(categories):
                                for j, key in enumerate(subcategories):
                                    dataset_y = self.data_rock[key]
                                    ax_rp_scatter[i][j].scatter(
                                        dataset_x, dataset_y, color=self.colors_gebpy["Option"],
                                        edgecolor="black", alpha=0.5)

                                    x_min, x_max, delta_x, y_min, y_max, delta_y = self.calculate_delta_xy(
                                        dataset_x=dataset_x, dataset_y=dataset_y)

                                    if delta_y < 1:
                                        n_digits_y = 3
                                    elif 1 <= delta_y < 5:
                                        n_digits_y = 2
                                    elif delta_y >= 5:
                                        n_digits_y = 0

                                    x_min = 0.995*x_min
                                    x_max = 1.005*x_max
                                    x_min, x_max = self.improve_xlimits_rock(delta_x=delta_x, x_min=x_min, x_max=x_max)
                                    y_min = 0.95*y_min
                                    y_max = 1.05*y_max
                                    y_min, y_max = self.improve_ylimits_rock(y_min=y_min, y_max=y_max)

                                    if key != "nu":
                                        if x_min < 0:
                                            x_min = 0
                                    if key != "nu":
                                        if y_min < 0:
                                            y_min = 0

                                    ax_rp_scatter[i][j].set_xlim(left=x_min, right=x_max)
                                    ax_rp_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                                    ax_rp_scatter[i][j].set_xticks(np.around(
                                        np.linspace(x_min, x_max, 4, dtype=int, endpoint=True), 0))
                                    ax_rp_scatter[i][j].set_yticks(np.around(
                                        np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), n_digits_y))
                                    ax_rp_scatter[i][j].xaxis.set_tick_params(labelsize=8)
                                    ax_rp_scatter[i][j].yaxis.set_tick_params(labelsize=8)
                                    ax_rp_scatter[i][j].set_xlabel("Density - kg/m$^3$", fontsize=8)
                                    ax_rp_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=8)

                                    ax_rp_scatter[i][j].grid(which="major", axis="both", linestyle="-")
                                    ax_rp_scatter[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                    ax_rp_scatter[i][j].minorticks_on()
                                    ax_rp_scatter[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rp_scatter[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rp_scatter[i][j].set_axisbelow(True)

                            if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                                #
                            else:
                                canvas_petrology = self.gui_elements["Temporary"]["Canvas"]["Petrology"]

                            canvas_petrology.draw()

                            self.gui_elements["Temporary"]["Axis"]["Rock Physics Scatter"] = ax_rp_scatter
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                        else:
                            ## Cleaning
                            for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Physics Scatter"]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("on")
                                    gui_axis.set_visible(True)
                            #
                            categories = ["Rock Physics Histogram", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                          "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                    else:
                        pass
                    #
                else:
                    ## RECONSTRUCTION
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Physics Scatter"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)
                    #
                    categories = ["Rock Physics Histogram", "Rock Chemistry Histogram", "Rock Chemistry Scatter",
                                  "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
            #
            self.last_rb_setting["Petrology"]["Rock Physics"].set(var_rb_diagram)
            #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 1:  # ROCK CHEMISTRY (MINERALS)
            var_rb_diagram = self.gui_variables["Radiobutton"]["Diagram Type Mineral"].get()
            if var_rb_diagram == 0:     # HISTOGRAM
                if self.last_rb_setting["Petrology"]["Rock Chemistry"].get() != 0:
                    if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                        categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                      "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                      "Rock Chemistry Scatter Element"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                    #
                    ## Diagram
                    if "Rock Chemistry Histogram" not in self.gui_elements["Temporary"]["Axis"]:
                        if "Petrology" not in self.gui_elements["Temporary"]["Figure"]:
                            fig_petrology = Figure(
                                figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                        else:
                            fig_petrology = self.gui_elements["Temporary"]["Figure"]["Petrology"]
                        #
                        ax_rc_histo = fig_petrology.subplots(nrows=3, ncols=3)
                        #
                        categories = []
                        labels = []
                        for index, mineral in enumerate(self.list_minerals_rock):
                            if index < 3:
                                if index == 0:
                                    categories.append([])
                                    labels.append([])
                                #
                                categories[0].append(mineral)
                                if mineral != "Urn":
                                    labels[0].append(str(mineral)+" (wt.%)")
                                else:
                                    labels[0].append(str(mineral) + " (ppm)")
                                #
                            elif 2 < index < 6:
                                if index == 3:
                                    categories.append([])
                                    labels.append([])
                                #
                                categories[1].append(mineral)
                                if mineral != "Urn":
                                    labels[1].append(str(mineral) + " (wt.%)")
                                else:
                                    labels[1].append(str(mineral) + " (ppm)")
                                #
                            elif 5 < index < 9:
                                if index == 6:
                                    categories.append([])
                                    labels.append([])
                                #
                                categories[2].append(mineral)
                                if mineral != "Urn":
                                    labels[2].append(str(mineral) + " (wt.%)")
                                else:
                                    labels[2].append(str(mineral) + " (ppm)")
                        #
                        for i, subcategories in enumerate(categories):
                            for j, key in enumerate(subcategories):
                                if key != "Urn":
                                    factor = 10**2
                                else:
                                    factor = 10**6
                                #
                                y, x, _ = ax_rc_histo[i][j].hist(
                                    x=np.array(self.data_rock["mineralogy"][key])*factor,
                                    color=self.colors_gebpy["Option"], edgecolor="black", bins=12)

                                x_min, x_max, delta_x, y_min, y_max, delta_y = self.calculate_delta_xy(
                                    dataset_x=x, dataset_y=y)
                                y_min = 0
                                y_max = round(1.05*max(y), 2)

                                if delta_x < 1:
                                    n_digits = 3
                                elif 1 <= delta_x < 5:
                                    n_digits = 2
                                elif delta_x >= 5:
                                    n_digits = 0

                                x_min, x_max = self.improve_xlimits_rock(delta_x=delta_x, x_min=x_min, x_max=x_max)
                                y_min, y_max = self.improve_ylimits_rock(y_min=y_min, y_max=y_max)

                                if x_min < 0:
                                    x_min = 0
                                if key != "Urn":
                                    if x_max > 100:
                                        x_max = 100

                                ax_rc_histo[i][j].set_xlim(left=x_min, right=x_max)
                                ax_rc_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                                ax_rc_histo[i][j].set_xticks(np.around(
                                    np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits))
                                ax_rc_histo[i][j].set_yticks(np.around(
                                    np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 1))
                                ax_rc_histo[i][j].xaxis.set_tick_params(labelsize=8)
                                ax_rc_histo[i][j].yaxis.set_tick_params(labelsize=8)
                                ax_rc_histo[i][j].set_xlabel(labels[i][j], fontsize=8)
                                ax_rc_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=8)

                                ax_rc_histo[i][j].grid(which="major", axis="both", linestyle="-")
                                ax_rc_histo[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                ax_rc_histo[i][j].minorticks_on()
                                ax_rc_histo[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                ax_rc_histo[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                ax_rc_histo[i][j].set_axisbelow(True)
                            if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                                #
                            else:
                                canvas_petrology = self.gui_elements["Temporary"]["Canvas"]["Petrology"]
                            #
                            canvas_petrology.draw()
                            #
                            self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Histogram"] = ax_rc_histo
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                        #
                    else:
                        ## Cleaning
                        for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Histogram"]:
                            for gui_axis in gui_axes:
                                gui_axis.axis("on")
                                gui_axis.set_visible(True)
                        #
                        categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rocks Chemistry Scatter",
                                      "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                    #
                elif self.last_rb_setting["Petrology"]["Rock Chemistry"].get() == 0 \
                        and self.last_rb_setting["Petrology"]["General Mode"].get() in [0, 2]:
                    ## RECONSTRUCTION
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Histogram"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)
                    #
                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rocks Chemistry Scatter",
                                  "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                    #
                else:
                    pass
                #
            elif var_rb_diagram == 1:   # SCATTER
                if self.last_rb_setting["Petrology"]["General Mode"].get() not in [0, 2]:
                    if self.last_rb_setting["Petrology"]["Rock Chemistry"].get() != 1:
                        if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                          "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                        if "Rock Chemistry Scatter" not in self.gui_elements["Temporary"]["Axis"]:
                            if "Petrology" not in self.gui_elements["Temporary"]["Figure"]:
                                fig_petrology = Figure(
                                    figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                            else:
                                fig_petrology = self.gui_elements["Temporary"]["Figure"]["Petrology"]
                            #
                            ax_rc_scatter = fig_petrology.subplots(nrows=3, ncols=3)
                            #
                            categories = []
                            labels = []
                            ref_mineral = None
                            ref_mean = 0
                            for index, mineral in enumerate(self.list_minerals_rock):
                                mineral_mean = np.mean(self.data_rock["mineralogy"][mineral])
                                if mineral_mean > ref_mean and mineral != "Urn":
                                    ref_mineral = mineral
                                    ref_mean = mineral_mean
                                #
                                if index < 3:
                                    if index == 0:
                                        categories.append([])
                                        labels.append([])
                                    #
                                    categories[0].append(mineral)
                                    if mineral != "Urn":
                                        labels[0].append(str(mineral) + " (wt.%)")
                                    else:
                                        labels[0].append(str(mineral) + " (ppm)")
                                    #
                                elif 2 < index < 6:
                                    if index == 3:
                                        categories.append([])
                                        labels.append([])
                                    #
                                    categories[1].append(mineral)
                                    if mineral != "Urn":
                                        labels[1].append(str(mineral) + " (wt.%)")
                                    else:
                                        labels[1].append(str(mineral) + " (ppm)")
                                    #
                                elif 5 < index < 9:
                                    if index == 6:
                                        categories.append([])
                                        labels.append([])
                                    #
                                    categories[2].append(mineral)
                                    if mineral != "Urn":
                                        labels[2].append(str(mineral) + " (wt.%)")
                                    else:
                                        labels[2].append(str(mineral) + " (ppm)")
                            #
                            dataset_x = np.array(self.data_rock["mineralogy"][ref_mineral])*100
                            #
                            for i, subcategories in enumerate(categories):
                                for j, key in enumerate(subcategories):
                                    if key != "Urn":
                                        factor = 10**2
                                        y_label = str(key)+" (wt.%)"
                                    else:
                                        factor = 10**6
                                        y_label = str(key) + " (ppm)"
                                    #
                                    dataset_y = np.array(self.data_rock["mineralogy"][key])*factor
                                    #
                                    ax_rc_scatter[i][j].scatter(
                                        dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black",
                                        alpha=0.5)

                                    x_min, x_max, delta_x, y_min, y_max, delta_y = self.calculate_delta_xy(
                                        dataset_x=dataset_x, dataset_y=dataset_y)

                                    x_min, x_max = self.improve_xlimits_rock(delta_x=delta_x, x_min=x_min, x_max=x_max)
                                    y_min, y_max = self.improve_ylimits_rock(y_min=y_min, y_max=y_max)

                                    if delta_x < 1:
                                        n_digits_x = 3
                                    elif 1 <= delta_x < 5:
                                        n_digits_x = 2
                                    elif delta_x >= 5:
                                        n_digits_x = 0
                                    #
                                    if delta_y < 1:
                                        n_digits_y = 3
                                    elif 1 <= delta_y < 5:
                                        n_digits_y = 2
                                    elif delta_y >= 5:
                                        n_digits_y = 0

                                    if x_min < 0:
                                        x_min = 0
                                    if y_min < 0:
                                        y_min = 0
                                    if key != "Urn":
                                        if x_max > 100:
                                            x_max = 100
                                        if y_max > 100:
                                            y_max = 100
                                    #
                                    ax_rc_scatter[i][j].set_xlim(left=x_min, right=x_max)
                                    ax_rc_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                                    ax_rc_scatter[i][j].set_xticks(np.around(
                                        np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits_x))
                                    ax_rc_scatter[i][j].set_yticks(np.around(
                                        np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), n_digits_y))
                                    ax_rc_scatter[i][j].xaxis.set_tick_params(labelsize=8)
                                    ax_rc_scatter[i][j].yaxis.set_tick_params(labelsize=8)
                                    ax_rc_scatter[i][j].set_xlabel(str(ref_mineral)+" (wt.%)", fontsize=8)
                                    ax_rc_scatter[i][j].set_ylabel(y_label, labelpad=0.5, fontsize=8)

                                    ax_rc_scatter[i][j].grid(which="major", axis="both", linestyle="-")
                                    ax_rc_scatter[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                    ax_rc_scatter[i][j].minorticks_on()
                                    ax_rc_scatter[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rc_scatter[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rc_scatter[i][j].set_axisbelow(True)
                            #
                            if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                                #
                            else:
                                canvas_petrology = self.gui_elements["Temporary"]["Canvas"]["Petrology"]
                            #
                            canvas_petrology.draw()
                            #
                            self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter"] = ax_rc_scatter
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                            #
                        else:
                            ## Cleaning
                            for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter"]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("on")
                                    gui_axis.set_visible(True)
                            #
                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                    else:
                        pass
                    #
                else:
                    ## RECONSTRUCTION
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)
                    #
                    categories = ["Rock Physics Histogram", "Rock Chemistry Histogram", "Rock Chemistry Histogram",
                                  "Rock Chemistry Histogram Element", "Rock Chemistry Scatter Element"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
            #
            self.last_rb_setting["Petrology"]["Rock Chemistry"].set(var_rb_diagram)
            #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 2:  # ROCK CHEMISTRY (ELEMENTS)
            var_rb_diagram = self.gui_variables["Radiobutton"]["Diagram Type Element"].get()
            if var_rb_diagram == 0:     # HISTOGRAM
                if self.last_rb_setting["Petrology"]["Rock Chemistry"].get() != 0 \
                        or self.last_rb_setting["Petrology"]["General Mode"].get() != 2:
                    if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                        categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                      "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                      "Rock Chemistry Scatter Element"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                    #
                    ## Diagram
                    if "Rock Chemistry Histogram Element" not in self.gui_elements["Temporary"]["Axis"]:
                        if "Petrology" not in self.gui_elements["Temporary"]["Figure"]:
                            fig_petrology = Figure(
                                figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                        else:
                            fig_petrology = self.gui_elements["Temporary"]["Figure"]["Petrology"]
                        #
                        ax_rce_histo = fig_petrology.subplots(nrows=3, ncols=3)
                        #
                        categories = []
                        labels = []
                        for index, element in enumerate(self.list_elements_rock):
                            if index < 3:
                                if index == 0:
                                    categories.append([])
                                    labels.append([])
                                #
                                categories[0].append(element)
                                if element != "U":
                                    labels[0].append(str(element)+" (wt.%)")
                                else:
                                    labels[0].append(str(element) + " (ppm)")
                                #
                            elif 2 < index < 6:
                                if index == 3:
                                    categories.append([])
                                    labels.append([])
                                #
                                categories[1].append(element)
                                if element != "U":
                                    labels[1].append(str(element) + " (wt.%)")
                                else:
                                    labels[1].append(str(element) + " (ppm)")
                                #
                            elif 5 < index < 9:
                                if index == 6:
                                    categories.append([])
                                    labels.append([])
                                #
                                categories[2].append(element)
                                if element != "U":
                                    labels[2].append(str(element) + " (wt.%)")
                                else:
                                    labels[2].append(str(element) + " (ppm)")
                        #
                        for i, subcategories in enumerate(categories):
                            for j, key in enumerate(subcategories):
                                if key != "U":
                                    factor = 10**2
                                else:
                                    factor = 10**6

                                y, x, _ = ax_rce_histo[i][j].hist(
                                    x=np.array(self.data_rock["chemistry"][key])*factor,
                                    color=self.colors_gebpy["Option"], edgecolor="black", bins=12)

                                x_min, x_max, delta_x, y_min, y_max, delta_y = self.calculate_delta_xy(
                                    dataset_x=x, dataset_y=y)
                                y_min = 0
                                y_max = round(1.05*max(y), 2)

                                if delta_x < 1:
                                    n_digits = 3
                                elif 1 <= delta_x < 5:
                                    n_digits = 2
                                elif delta_x >= 5:
                                    n_digits = 0

                                x_min, x_max = self.improve_xlimits_rock(delta_x=delta_x, x_min=x_min, x_max=x_max)
                                y_min, y_max = self.improve_ylimits_rock(y_min=y_min, y_max=y_max)

                                if x_min < 0:
                                    x_min = 0
                                if key != "U":
                                    if x_max > 100:
                                        x_max = 100

                                ax_rce_histo[i][j].set_xlim(left=x_min, right=x_max)
                                ax_rce_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                                ax_rce_histo[i][j].set_xticks(np.around(
                                    np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits))
                                ax_rce_histo[i][j].set_yticks(np.around(
                                    np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 1))
                                ax_rce_histo[i][j].xaxis.set_tick_params(labelsize=8)
                                ax_rce_histo[i][j].yaxis.set_tick_params(labelsize=8)
                                ax_rce_histo[i][j].set_xlabel(labels[i][j], fontsize=8)
                                ax_rce_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=8)

                                ax_rce_histo[i][j].grid(which="major", axis="both", linestyle="-")
                                ax_rce_histo[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                ax_rce_histo[i][j].minorticks_on()
                                ax_rce_histo[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                ax_rce_histo[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                ax_rce_histo[i][j].set_axisbelow(True)
                            if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                                #
                            else:
                                canvas_petrology = self.gui_elements["Temporary"]["Canvas"]["Petrology"]
                            #
                            canvas_petrology.draw()
                            #
                            self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Histogram Element"] = ax_rce_histo
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                        #
                    else:
                        ## Cleaning
                        for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Histogram Element"]:
                            for gui_axis in gui_axes:
                                gui_axis.axis("on")
                                gui_axis.set_visible(True)
                        #
                        categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rocks Chemistry Histogram",
                                      "Rocks Chemistry Scatter", "Rock Chemistry Scatter Element"]
                        for category in categories:
                            if category in self.gui_elements["Temporary"]["Axis"]:
                                for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                    for gui_axis in gui_axes:
                                        gui_axis.axis("off")
                                        gui_axis.set_visible(False)
                        #
                        self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                    #
                elif self.last_rb_setting["Petrology"]["Rock Chemistry"].get() == 0 \
                        and self.last_rb_setting["Petrology"]["General Mode"].get() in [0, 2]:
                    ## RECONSTRUCTION
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Histogram Element"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)
                    #
                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rocks Chemistry Histogram",
                                  "Rocks Chemistry Scatter", "Rock Chemistry Scatter Element"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                    #
                else:
                    pass
                #
            elif var_rb_diagram == 1:   # SCATTER
                if self.last_rb_setting["Petrology"]["General Mode"].get() not in [0, 1]:
                    if self.last_rb_setting["Petrology"]["Rock Chemistry"].get() != 1:
                        if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element",
                                          "Rock Chemistry Scatter Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)
                            #
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                        #
                        if "Rock Chemistry Scatter Element" not in self.gui_elements["Temporary"]["Axis"]:
                            if "Petrology" not in self.gui_elements["Temporary"]["Figure"]:
                                fig_petrology = Figure(
                                    figsize=(3, 3), dpi=150, tight_layout=True, facecolor=self.colors_gebpy["Background"])
                            else:
                                fig_petrology = self.gui_elements["Temporary"]["Figure"]["Petrology"]
                            #
                            ax_rce_scatter = fig_petrology.subplots(nrows=3, ncols=3)
                            #
                            categories = []
                            labels = []
                            ref_element = None
                            ref_mean = 0
                            for index, element in enumerate(self.list_elements_rock):
                                element_mean = np.mean(self.data_rock["chemistry"][element])
                                if element_mean > ref_mean and element != "U":
                                    ref_element = element
                                    ref_mean = element_mean
                                #
                                if index < 3:
                                    if index == 0:
                                        categories.append([])
                                        labels.append([])
                                    #
                                    categories[0].append(element)
                                    if element != "U":
                                        labels[0].append(str(element) + " (wt.%)")
                                    else:
                                        labels[0].append(str(element) + " (ppm)")
                                    #
                                elif 2 < index < 6:
                                    if index == 3:
                                        categories.append([])
                                        labels.append([])
                                    #
                                    categories[1].append(element)
                                    if element != "U":
                                        labels[1].append(str(element) + " (wt.%)")
                                    else:
                                        labels[1].append(str(element) + " (ppm)")
                                    #
                                elif 5 < index < 9:
                                    if index == 6:
                                        categories.append([])
                                        labels.append([])
                                    #
                                    categories[2].append(element)
                                    if element != "U":
                                        labels[2].append(str(element) + " (wt.%)")
                                    else:
                                        labels[2].append(str(element) + " (ppm)")
                            #
                            dataset_x = np.array(self.data_rock["chemistry"][ref_element])*100
                            #
                            for i, subcategories in enumerate(categories):
                                for j, key in enumerate(subcategories):
                                    if key != "U":
                                        factor = 10**2
                                        y_label = str(key)+" (wt.%)"
                                    else:
                                        factor = 10**6
                                        y_label = str(key) + " (ppm)"

                                    dataset_y = np.array(self.data_rock["chemistry"][key])*factor
                                    ax_rce_scatter[i][j].scatter(
                                        dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black",
                                        alpha=0.5)

                                    x_min, x_max, delta_x, y_min, y_max, delta_y = self.calculate_delta_xy(
                                        dataset_x=dataset_x, dataset_y=dataset_y)

                                    if delta_x < 1:
                                        n_digits_x = 3
                                    elif 1 <= delta_x < 5:
                                        n_digits_x = 2
                                    elif delta_x >= 5:
                                        n_digits_x = 0

                                    if delta_y < 1:
                                        n_digits_y = 3
                                    elif 1 <= delta_y < 5:
                                        n_digits_y = 2
                                    elif delta_y >= 5:
                                        n_digits_y = 0

                                    x_min, x_max = self.improve_xlimits_rock(delta_x=delta_x, x_min=x_min, x_max=x_max)
                                    y_min, y_max = self.improve_ylimits_rock(y_min=y_min, y_max=y_max)

                                    if x_min < 0:
                                        x_min = 0
                                    if y_min < 0:
                                        y_min = 0
                                    if key != "U":
                                        if x_max > 100:
                                            x_max = 100
                                        if y_max > 100:
                                            y_max = 100

                                    ax_rce_scatter[i][j].set_xlim(left=x_min, right=x_max)
                                    ax_rce_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                                    ax_rce_scatter[i][j].set_xticks(np.around(
                                        np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits_x))
                                    ax_rce_scatter[i][j].set_yticks(np.around(
                                        np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), n_digits_y))
                                    ax_rce_scatter[i][j].xaxis.set_tick_params(labelsize=8)
                                    ax_rce_scatter[i][j].yaxis.set_tick_params(labelsize=8)
                                    ax_rce_scatter[i][j].set_xlabel(str(ref_element)+" (wt.%)", fontsize=8)
                                    ax_rce_scatter[i][j].set_ylabel(y_label, labelpad=0.5, fontsize=8)

                                    ax_rce_scatter[i][j].grid(which="major", axis="both", linestyle="-")
                                    ax_rce_scatter[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                    ax_rce_scatter[i][j].minorticks_on()
                                    ax_rce_scatter[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rce_scatter[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4))
                                    ax_rce_scatter[i][j].set_axisbelow(True)

                            if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                            else:
                                canvas_petrology = self.gui_elements["Temporary"]["Canvas"]["Petrology"]

                            canvas_petrology.draw()

                            self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter Element"] = ax_rce_scatter
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology

                        else:
                            ## Cleaning
                            for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter Element"]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("on")
                                    gui_axis.set_visible(True)

                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element"]
                            for category in categories:
                                if category in self.gui_elements["Temporary"]["Axis"]:
                                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                        for gui_axis in gui_axes:
                                            gui_axis.axis("off")
                                            gui_axis.set_visible(False)

                            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()

                    else:
                        pass
                else:
                    ## RECONSTRUCTION
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter Element"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)

                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                  "Rock Chemistry Scatter", "Rock Chemistry Histogram Element"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)

                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()

            self.last_rb_setting["Petrology"]["Rock Chemistry"].set(var_rb_diagram)

    def change_rb_analysis_rocks(self):
        start_column = 35
        var_rb_mode = self.gui_variables["Radiobutton"]["Analysis Mode"].get()
        #
        if var_rb_mode == 0:   # ROCK PHYSICS
            if self.last_rb_setting["Petrology"]["General Mode"].get() != 0:
                ## Cleaning
                if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                  "Rock Chemistry Scatter"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                #
                for key, gui_items in self.gui_elements["Temporary"].items():
                    if len(gui_items) > 0:
                        if key not in ["Canvas"]:
                            if type(gui_items) == list:
                                for gui_item in gui_items:
                                    gui_item.grid_remove()
                                gui_items.clear()
                #
                for key, gui_items in self.gui_elements["Rockbuilder Temporary"].items():
                    if len(gui_items) > 0:
                        if key not in ["Canvas"]:
                            if type(gui_items) == list:
                                for gui_item in gui_items:
                                    gui_item.grid_remove()
                                gui_items.clear()
                #
                ## Labels
                lbl_addsetup = SimpleElements(
                    parent=self.parent, row_id=32, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
                    fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Additional Settings", font_option="sans 14 bold", relief=tk.FLAT)
                lbl_diagram_type = SimpleElements(
                    parent=self.parent, row_id=34, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
                #
                self.gui_elements["Rockbuilder Temporary"]["Label"].extend([lbl_addsetup, lbl_diagram_type])
                #
                ## Radiobuttons
                rb_diagram_type_01 = SimpleElements(
                    parent=self.parent, row_id=34, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Histogram", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Rock"], value_rb=0,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                rb_diagram_type_02 = SimpleElements(
                    parent=self.parent, row_id=36, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Scatter", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Rock"], value_rb=1,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                #
                self.gui_elements["Rockbuilder Temporary"]["Radiobutton"].extend(
                    [rb_diagram_type_01, rb_diagram_type_02])
                #
                self.change_rb_diagram_rocks()
                #
                self.last_rb_analysis_rock.set(var_rb_mode)
                #
            else:
                pass
        #
        elif var_rb_mode == 1: # ROCK CHEMISTRY
            if self.last_rb_setting["Petrology"]["General Mode"].get() != 1:
                ## Cleaning
                if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                  "Rock Chemistry Scatter"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                #
                for key, gui_items in self.gui_elements["Rockbuilder Temporary"].items():
                    if len(gui_items) > 0:
                        if key not in ["Canvas"]:
                            if type(gui_items) == list:
                                for gui_item in gui_items:
                                    gui_item.grid_remove()
                                gui_items.clear()
                #
                ## Labels
                lbl_addsetup = SimpleElements(
                    parent=self.parent, row_id=32, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
                    fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Additional Settings", font_option="sans 14 bold", relief=tk.FLAT)
                lbl_diagram_type = SimpleElements(
                    parent=self.parent, row_id=34, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
                lbl_mineral = SimpleElements(
                    parent=self.parent, row_id=38, column_id=0, n_rows=2, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Mineral Selection", font_option="sans 10 bold", relief=tk.FLAT)
                #
                self.gui_elements["Rockbuilder Temporary"]["Label"].extend(
                    [lbl_addsetup, lbl_diagram_type, lbl_mineral])
                #
                ## Radiobuttons
                rb_diagram_type_01 = SimpleElements(
                    parent=self.parent, row_id=34, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Histogram", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Mineral"], value_rb=0,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                rb_diagram_type_02 = SimpleElements(
                    parent=self.parent, row_id=36, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Scatter", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Mineral"], value_rb=1,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                #
                self.gui_elements["Rockbuilder Temporary"]["Radiobutton"].extend(
                    [rb_diagram_type_01, rb_diagram_type_02])
                #
                ## Option Menu
                self.list_minerals_rock.sort()
                list_opt = self.list_minerals_rock
                opt_mineral = SimpleElements(
                    parent=self.parent, row_id=38, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_option_menu(
                    var_opt=self.gui_variables["Radiobutton"]["Amount Mineral"],
                    var_opt_set=self.gui_variables["Radiobutton"]["Amount Mineral"].get(), opt_list=list_opt,
                    active_bg=self.colors_gebpy["Accent"])
                #
                self.gui_elements["Rockbuilder Temporary"]["Option Menu"].extend([opt_mineral])
                #
                self.change_rb_diagram_rocks()
                #
                self.last_rb_analysis_rock.set(var_rb_mode)
                #
            else:
                pass
        #
        elif var_rb_mode == 2: # Element Composition
            if self.last_rb_setting["Petrology"]["General Mode"].get() != 2:
                ## Cleaning
                if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                  "Rock Chemistry Scatter"]
                    for category in categories:
                        if category in self.gui_elements["Temporary"]["Axis"]:
                            for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("off")
                                    gui_axis.set_visible(False)
                    #
                    self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
                #
                for key, gui_items in self.gui_elements["Rockbuilder Temporary"].items():
                    if len(gui_items) > 0:
                        if key not in ["Canvas"]:
                            if type(gui_items) == list:
                                for gui_item in gui_items:
                                    gui_item.grid_remove()
                                gui_items.clear()
                #
                ## Labels
                lbl_addsetup = SimpleElements(
                    parent=self.parent, row_id=32, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
                    fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Additional Settings", font_option="sans 14 bold", relief=tk.FLAT)
                lbl_diagram_type = SimpleElements(
                    parent=self.parent, row_id=34, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
                lbl_concentration_setup = SimpleElements(
                    parent=self.parent, row_id=38, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Concentration Setup", font_option="sans 10 bold", relief=tk.FLAT)
                lbl_element = SimpleElements(
                    parent=self.parent, row_id=42, column_id=0, n_rows=2, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Element Selection", font_option="sans 10 bold", relief=tk.FLAT)
                lbl_oxide = SimpleElements(
                    parent=self.parent, row_id=44, column_id=0, n_rows=2, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Oxide Selection", font_option="sans 10 bold", relief=tk.FLAT)
                #
                self.gui_elements["Rockbuilder Temporary"]["Label"].extend(
                    [lbl_addsetup, lbl_diagram_type, lbl_concentration_setup, lbl_element, lbl_oxide])
                #
                ## Radiobuttons
                rb_diagram_type_01 = SimpleElements(
                    parent=self.parent, row_id=34, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Histogram", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Element"], value_rb=0,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                rb_diagram_type_02 = SimpleElements(
                    parent=self.parent, row_id=36, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Scatter", var_rb=self.gui_variables["Radiobutton"]["Diagram Type Element"], value_rb=1,
                    color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                rb_concentration_type_01 = SimpleElements(
                    parent=self.parent, row_id=38, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Element Concentrations", var_rb=self.gui_variables["Radiobutton"]["Concentration Type"],
                    value_rb=0, color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                rb_concentration_type_02 = SimpleElements(
                    parent=self.parent, row_id=40, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                    text="Oxide Concentrations", var_rb=self.gui_variables["Radiobutton"]["Concentration Type"],
                    value_rb=1, color_bg=self.colors_gebpy["Background"], command=self.change_rb_diagram_rocks)
                #
                self.gui_elements["Rockbuilder Temporary"]["Radiobutton"].extend(
                    [rb_diagram_type_01, rb_diagram_type_02, rb_concentration_type_01, rb_concentration_type_02])
                #
                ## Option Menu
                self.list_elements_rock.sort()
                list_opt_element = self.list_elements_rock
                list_opt_oxides = self.find_suitable_oxides(var_list_elements=list_opt_element)
                #
                opt_element = SimpleElements(
                    parent=self.parent, row_id=42, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_option_menu(
                    var_opt=self.gui_variables["Radiobutton"]["Amount Element"],
                    var_opt_set=self.gui_variables["Radiobutton"]["Amount Element"].get(), opt_list=list_opt_element,
                    active_bg=self.colors_gebpy["Accent"])
                opt_oxide = SimpleElements(
                    parent=self.parent, row_id=44, column_id=14, n_rows=2, n_columns=16,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_option_menu(
                    var_opt=self.gui_variables["Radiobutton"]["Amount Oxide"],
                    var_opt_set=self.gui_variables["Radiobutton"]["Amount Oxide"].get(), opt_list=list_opt_oxides,
                    active_bg=self.colors_gebpy["Accent"])
                #
                self.gui_elements["Rockbuilder Temporary"]["Option Menu"].extend([opt_element, opt_oxide])
                #
                self.change_rb_diagram_rocks()
                #
                self.last_rb_analysis_rock.set(var_rb_mode)
                #
            else:
                pass
        #
        self.last_rb_setting["Petrology"]["General Mode"].set(var_rb_mode)
    #
    ######################################
    ## G e n e r a l  F u n c t i o n s ##
    ######################################
    #
    def calculate_fractions(self, number=100, mode="simplified"):
        self.rock_fractions = {"Rock Forming": {}, "Ore": {}, "Clay": {}}
        #
        if mode == "simplified":
            rock_forming_min = self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Minimum"].get()
            rock_forming_max = self.gui_variables["Entry"]["Mineralogy"]["Rock Forming Total"]["Maximum"].get()
            n_rock_forming = 0
            mean_rock_forming = (rock_forming_min + rock_forming_max)/2
            ore_min = self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Minimum"].get()
            ore_max = self.gui_variables["Entry"]["Mineralogy"]["Ore Total"]["Maximum"].get()
            n_ore = 0
            mean_ore = (ore_min + ore_max)/2
            clay_min = self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Minimum"].get()
            clay_max = self.gui_variables["Entry"]["Mineralogy"]["Clay Total"]["Maximum"].get()
            n_clay = 0
            mean_clay = (clay_min + clay_max)/2
            #
            self.rock_fractions["Rock Forming"]["Min"] = float(rock_forming_min)
            self.rock_fractions["Rock Forming"]["Max"] = float(rock_forming_max)
            self.rock_fractions["Ore"]["Min"] = float(ore_min)
            self.rock_fractions["Ore"]["Max"] = float(ore_max)
            self.rock_fractions["Clay"]["Min"] = float(clay_min)
            self.rock_fractions["Clay"]["Max"] = float(clay_max)
            #
            self.fractions_sum = {"Rock Forming": 0, "Ore": 0, "Clay": 0}
            fractions_mean = {"Rock Forming": mean_rock_forming, "Ore": mean_ore, "Clay": mean_clay}
            fractions_mean = {k: v for k, v in sorted(fractions_mean.items(), key=lambda item: item[1], reverse=True)}
            order_fractions = {}
            #
            magicnumber_set = False
            for index, (key, value) in enumerate(fractions_mean.items()):
                if index == 0:
                    order_fractions[index] = key
                else:
                    if magicnumber_set == False:
                        magicnumber = rd.randint(0, 1)
                        if magicnumber == 1:
                            magicnumber_set = True
                    else:
                        magicnumber = 0
                    #
                    order_fractions[index + magicnumber] = key

            for mineral, value in self.gui_variables["Checkbox"]["Mineralogy"].items():
                if value.get() == 1:
                    var_min = self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral].get()
                    var_max = self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral].get()
                    var_mean = (var_min + var_max)/2
                    #
                    if mineral in self.rock_forming_minerals:
                        var_group = "Rock Forming"
                    elif mineral in self.ore_minerals:
                        var_group = "Ore"
                    elif mineral in self.clay_minerals:
                        var_group = "Clay"
                    #
                    self.selected_minerals[mineral] = {
                        "Min": var_min, "Max": var_max, "Mean": var_mean, "Group": var_group}
                    #
                    if mineral in self.rock_forming_minerals:
                        n_rock_forming += 1
                        self.fractions_sum["Rock Forming"] += 1
                    elif mineral in self.ore_minerals:
                        n_ore += 1
                        self.fractions_sum["Ore"] += 1
                    elif mineral in self.clay_minerals:
                        n_clay += 1
                        self.fractions_sum["Clay"] += 1
            #
            self.fractions_sum = dict(sorted(self.fractions_sum.items(), key=lambda item: item[1]))
            #
            w_total = 100
            w_actual = 0
            previous_keys = []
            for key, value in self.fractions_sum.items():
                if value > 0:
                    for previous_key in previous_keys:
                        w_actual += self.rock_fractions[previous_key]["Effective"]
                    #
                    value_calc = round(w_total - w_actual, 4)
                    #
                else:
                    value_calc = round(0.0, 4)
                #
                self.rock_fractions[key]["Effective"] = value_calc
                self.rock_fractions[key]["Size"] = value
                previous_keys.append(key)
            #
            previous_keys = []
            w_effective = 0
            temp_min = 0
            temp_max = 0
            for index, (key, value) in enumerate(self.fractions_sum.items()):
                if self.rock_fractions[key]["Size"] > 0:
                    if index < 2:
                        self.rock_fractions[key]["Effective"] = round(
                            (self.rock_fractions[key]["Min"] + self.rock_fractions[key]["Max"])/2, 4)
                        temp_min += self.rock_fractions[key]["Min"]
                        temp_max += self.rock_fractions[key]["Max"]
                    else:
                        var_delta = 100 - w_effective
                        value_min = 100 - temp_max
                        value_max = 100 - temp_min
                        if value_min < 0:
                            value_min = 0
                        if value_max < 0:
                            value_max = 0
                        value_eff = (value_min + value_max)/2
                        #
                        self.rock_fractions[key]["Min"] = round(value_min, 4)
                        self.rock_fractions[key]["Max"] = round(value_max, 4)
                        self.rock_fractions[key]["Effective"] = round(value_eff, 4)
                else:
                    self.rock_fractions[key]["Min"] = 0.0
                    self.rock_fractions[key]["Max"] = 0.0
                    self.rock_fractions[key]["Effective"] = round(
                        (self.rock_fractions[key]["Min"] + self.rock_fractions[key]["Max"])/2, 4)
                    temp_min += self.rock_fractions[key]["Min"]
                    temp_max += self.rock_fractions[key]["Max"]
                #
                w_effective += self.rock_fractions[key]["Effective"]
                previous_keys.append(key)
            amounts_fractions = {}
            last_amounts = 0
            for index, fraction in enumerate(fractions_mean):
                if self.fractions_sum[fraction] > 0:
                    var_key = fraction+" Total"
                    if index == 0:
                        var_min = self.gui_variables["Entry"]["Mineralogy"][var_key]["Minimum"].get()
                        var_max = self.gui_variables["Entry"]["Mineralogy"][var_key]["Maximum"].get()
                        var_amount = round(rd.uniform(var_min, var_max), 4)
                    elif index == len(fractions_mean) - 1:
                        var_amount = round((100 - last_amounts), 4)
                    else:
                        var_min = self.gui_variables["Entry"]["Mineralogy"][var_key]["Minimum"].get()
                        var_max = 100 - last_amounts
                        var_amount = round(rd.uniform(var_min, var_max), 4)
                    #
                    amounts_fractions[fraction] = round(var_amount, 4)
                    last_amounts += round(var_amount, 4)
                else:
                    var_amount = 0.0
                    amounts_fractions[fraction] = round(var_amount, 4)
                    last_amounts += round(var_amount, 4)
            #
            for fraction, amount in amounts_fractions.items():
                amounts_fractions[fraction] = round(amount/100, 4)
            #
            amount_fraction_rf = round(amounts_fractions["Rock Forming"], 4)
            amount_fraction_o = round(amounts_fractions["Ore"], 4)
            amount_fraction_c = round(amounts_fractions["Clay"], 4)
            #
            sorted_minerals = {}
            for fraction, size in self.fractions_sum.items():
                sorted_minerals[fraction] = {}
                if size > 0:
                    for mineral, values in self.selected_minerals.items():
                        if fraction == values["Group"]:
                            sorted_minerals[fraction][mineral] = values["Mean"]
            self.minerals_helper = {}
            for index_01, (key_01, value_01) in enumerate(self.fractions_sum.items()):
                if value_01 > 0:
                    self.minerals_helper[key_01] = {}
                    for index_02, (key_02, value_02) in enumerate(self.selected_minerals.items()):
                        if key_01 == value_02["Group"]:
                            main_mineral = self.find_maximum_in_dict(var_dict=sorted_minerals[key_01])
                            self.minerals_helper[key_01][key_02] = {}
                            if value_01 == 1:
                                value_min = self.rock_fractions[key_01]["Min"]*value_02["Min"]/100
                                value_max = self.rock_fractions[key_01]["Max"]*value_02["Max"]/100
                                value_eff = (value_min + value_max)/2
                                self.minerals_helper[key_01][key_02]["Min"] = round(value_min, 4)
                                self.minerals_helper[key_01][key_02]["Max"] = round(value_max, 4)
                                self.minerals_helper[key_01][key_02]["Effective"] = round(value_eff, 4)
                                self.minerals_helper[key_01][key_02]["Fraction"] = round(
                                    self.rock_fractions[key_01]["Effective"], 4)
                                #
                                self.selected_minerals[key_02]["Min"] = round(value_min, 4)
                                self.selected_minerals[key_02]["Max"] = round(value_max, 4)
                                self.selected_minerals[key_02]["Mean"] = round(value_eff, 4)
                                #
                            else:
                                value_min = self.rock_fractions[key_01]["Min"]*value_02["Min"]/100
                                value_max = self.rock_fractions[key_01]["Max"]*value_02["Max"]/100
                                value_eff = (value_min + value_max)/2
                                self.minerals_helper[key_01][key_02]["Min"] = round(value_min, 4)
                                self.minerals_helper[key_01][key_02]["Max"] = round(value_max, 4)
                                self.minerals_helper[key_01][key_02]["Effective"] = round(value_eff, 4)
                                self.minerals_helper[key_01][key_02]["Fraction"] = round(
                                    self.rock_fractions[key_01]["Effective"], 4)
                                #
                                self.selected_minerals[key_02]["Min"] = round(value_min, 4)
                                self.selected_minerals[key_02]["Max"] = round(value_max, 4)
                                self.selected_minerals[key_02]["Mean"] = round(value_eff, 4)
        #
        ## TESTING
        # print("")
        # print("Rock Fractions:")
        # for key, value in self.rock_fractions.items():
        #     print(key, value)
        # print("")
        # print("Minerals Helper:")
        # for key, value in self.minerals_helper.items():
        #     print(key, value)
        # print("")
        # print("Selected Minerals:")
        # for key, value in self.selected_minerals.items():
        #     print(key, value)
        # print("")
    #
    def find_suitable_oxides(self, var_list_elements):
        list_oxides = []
        #
        for element in var_list_elements:
            if element in ["H", "Li", "C", "Na", "Cl", "K", "Cu", "Br", "Rb", "Ag", "I", "Cs", "Au", "Hg", "Tl", "At",
                           "Fr"]:   # +1
                list_oxides.append(element+str("2O"))
                #
            if element in ["C", "Cu", "Hg", "Be", "Mg", "S", "Ca", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Zn",
                           "Ge", "Se", "Sr", "Pd", "Cd", "Sn", "Te", "Ba", "Pt", "Pb", "Eu", "Po", "Rn", "Ra", "No",
                           "Cn"]: # +2
                list_oxides.append(element + str("O"))
                #
            if element in ["C", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Eu", "Cl", "Br", "I", "Au", "Tl", "B", "N", "Al",
                           "P", "Sc", "Ga", "As", "Y", "Ru", "Rh", "In", "Sb", "La", "Pr", "Sm", "Gd", "Tb", "Tm", "Yb",
                           "Ir", "Bi", "Ce", "Nd", "Pm", "Dy", "Ho" "Er", "Lu", "Ac", "Am", "Bk", "Cf", "Es", "Fm",
                           "Md", "Cm", "Lr"]: # +3
                list_oxides.append(element + str("2O3"))
                #
            if element in ["C", "Ti", "V", "Mn", "Ru", "Ir", "Ce", "S", "Ge", "Se", "Pd", "Sn", "Te", "Pt", "Pb",
                           "Si", "Zr", "Mo", "Tc", "Hf", "W", "Re", "Os", "Th", "U", "Pu", "Po", "Rf"]: # +4
                list_oxides.append(element + str("O2"))
                #
            if element in ["V", "Cl", "Br", "I", "N", "P", "As", "Sb", "Nb", "Ta", "Pa", "Np", "Db"]: # +5
                list_oxides.append(element + str("2O5"))
                #
            if element in ["Mn", "S", "Se", "Te", "Mo", "W", "U", "Cr", "Sg"]: # +6
                list_oxides.append(element + str("O3"))
                #
            if element in ["Mn", "Cl", "I", "Tc", "Re", "Bh"]: # +7
                list_oxides.append(element + str("2O7"))
                #
            if element in ["Hs"]: # +8
                list_oxides.append(element + str("O4"))
        #
        return list_oxides
    #
    def find_maximum_in_dict(self, var_dict):
        v = list(var_dict.values())
        k = list(var_dict.keys())
        #
        return k[v.index(max(v))]
    #
    def select_all_checkboxes(self, var_cb):
        for key, cb_var in var_cb.items():
            cb_var.set(1)
    #
    def unselect_all_checkboxes(self, var_cb):
        for key, cb_var in var_cb.items():
            cb_var.set(0)
    #
    def restart_gebpy(self):
        self.parent.destroy()
        #
        root = tk.Tk()
        #
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        #
        GebPyGUI(parent=root, var_screen_width=screen_width, var_screen_height=screen_height)
        #
        root.mainloop()
    #
    def export_mineral_data(self, var_dataset, var_name):
        list_keys = list(var_dataset.keys())
        if "mineral" in list_keys:
            list_keys.remove("mineral")
        if "state" in list_keys:
            list_keys.remove("state")
        if "chemistry" in list_keys:
            list_keys.remove("chemistry")

        if "compounds" in list_keys:
            list_keys.remove("compounds")
        if "LA-ICP-MS" in list_keys:
            list_keys.remove("LA-ICP-MS")

        list_keys = ["POISSON" if item == "nu" else item for item in list_keys]
        list_elements = list(var_dataset["chemistry"].keys())
        if "compounds" in list_keys:
            list_compounds = list(var_dataset["compounds"].keys())

        report_file = filedialog.asksaveasfile(
            mode="w", initialfile="Report_Mineral_"+str(var_name), defaultextension=".csv")

        ## General Data
        report_file.write("REPORT (MINERALOGY)"+"\n")
        report_file.write("Mineral"+";"+str(var_name)+"\n")
        report_file.write("\n")

        ## Geophysical Data
        report_file.write("MINERAL DATA" + "\n")
        raw_line = "ID;"
        for key in list_keys:
            raw_line += str(key)
            raw_line += str(";")

        raw_line += str(";")
        for element in list_elements:
            raw_line += str(element)
            raw_line += str(";")

        if "compounds" in list_keys:
            raw_line += str(";")
            for compound in list_compounds:
                raw_line += str(compound)
                raw_line += str(";")

        raw_line += str("\n")
        report_file.write(raw_line)

        index = 0
        while index < len(var_dataset["rho"]):
            raw_line = str(index + 1) + ";"
            for key, values in var_dataset.items():
                if key in list_keys:
                    try:
                        raw_line += str(round(values[index], 3))
                    except:
                        raw_line += str(values[index])
                    raw_line += str(";")
                elif key == "nu":
                    raw_line += str(round(values[index], 3))
                    raw_line += str(";")

            raw_line += str(";")
            for element, values in var_dataset["chemistry"].items():
                if element in list_elements:
                    if element not in ["U"]:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")
                    else:
                        raw_line += str(round(values[index], 6))
                        raw_line += str(";")

            if "compounds" in list_keys:
                raw_line += str(";")
                for compound, values in var_dataset["compounds"].items():
                    if compound in list_compounds:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")

            report_file.write(raw_line+"\n")

            index += 1
        report_file.write("\n")

    def export_rock_data(self, var_dataset, var_name):
        list_keys = list(var_dataset.keys())
        list_keys.remove("mineralogy")
        list_keys.remove("chemistry")
        if "compounds" in list_keys:
            list_keys.remove("compounds")
        if "fluid" in list_keys:
            list_keys.remove("fluid")

        list_keys = ["POISSON" if item == "nu" else item for item in list_keys]
        list_minerals = list(var_dataset["mineralogy"].keys())
        list_elements = list(var_dataset["chemistry"].keys())
        if "compounds" in list_keys:
            list_compounds = list(var_dataset["compounds"].keys())

        report_file = filedialog.asksaveasfile(
            mode="w", initialfile="Report_Rock_"+str(var_name), defaultextension=".csv")

        ## General Data
        report_file.write("REPORT (PETROLOGY)" + "\n")
        report_file.write("Rock" + ";" + str(var_name) + "\n")
        report_file.write("\n")

        ## Geophysical Data
        report_file.write("ROCK DATA" + "\n")
        raw_line = "ID;"
        for key in list_keys:
            raw_line += str(key)
            raw_line += str(";")

        raw_line += str(";")
        for mineral in list_minerals:
            raw_line += str(mineral)
            raw_line += str(";")

        raw_line += str(";")
        for element in list_elements:
            raw_line += str(element)
            raw_line += str(";")

        if "compounds" in list_keys:
            raw_line += str(";")
            for compound in list_compounds:
                raw_line += str(compound)
                raw_line += str(";")

        raw_line += str("\n")

        report_file.write(raw_line)

        index = 0
        while index < len(var_dataset["rho"]):
            raw_line = str(index + 1) + ";"
            for key, values in var_dataset.items():
                if key in list_keys:
                    try:
                        raw_line += str(round(values[index], 3))
                    except:
                        raw_line += str(values[index])
                    raw_line += str(";")
                elif key == "nu":
                    raw_line += str(round(values[index], 3))
                    raw_line += str(";")

            raw_line += str(";")

            for mineral, values in var_dataset["mineralogy"].items():
                if mineral in list_minerals:
                    if mineral not in ["Urn"]:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")
                    else:
                        raw_line += str(round(values[index], 6))
                        raw_line += str(";")

            raw_line += str(";")
            for element, values in var_dataset["chemistry"].items():
                if element in list_elements:
                    if element not in ["U"]:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")
                    else:
                        raw_line += str(round(values[index], 6))
                        raw_line += str(";")

            if "compounds" in list_keys:
                raw_line += str(";")
                for compound, values in var_dataset["compounds"].items():
                    if compound in list_compounds:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")

            report_file.write(raw_line + "\n")

            index += 1

        report_file.write("\n")

    def calculate_mineral_fractions(self, var_minerals, var_data, var_n):
        ## Input
        # print("var_minerals:", var_minerals)
        # print("var_data:", var_data)
        # print("var_n:", var_n)
        # print("selected_minerals:", self.selected_minerals)

        data_amounts = {}
        n_minerals = len(var_minerals)
        mineral_information = {}
        for key_group, group_data in var_data.items():
            for mineral, values in group_data.items():
                if mineral not in data_amounts:
                    data_amounts[mineral] = []
                    mineral_information[mineral] = values
        n = 0
        while n < var_n:
            condition = False
            while condition == False:
                temp_values = {}
                phi_total = 0
                for index, (mineral, values) in enumerate(self.selected_minerals.items()):
                    limit_lower = values["Min"]
                    limit_upper = values["Max"]
                    if limit_lower < 0:
                        limit_lower = 0.0

                    if limit_upper > 100:
                        limit_upper = 100

                    if index < n_minerals - 1:
                        #value = round(0.01*values["Fraction"]*rd.randint(values["Min"], values["Max"]), 6)
                        value = round(rd.uniform(values["Min"], values["Max"]), 6)
                    else:
                        value = round(100 - phi_total, 6)

                    # limit_lower = 0.01*values["Fraction"]*values["Min"]
                    # limit_upper = 0.01*values["Fraction"]*values["Max"]
                    # limit_lower = values["Min"]
                    # limit_upper = values["Max"]
                    #
                    # if limit_lower < 0:
                    #     limit_lower = 0.0
                    #
                    # if limit_upper > 100:
                    #     limit_upper = 100

                    if limit_lower <= value <= limit_upper:
                        phi_total += value
                        temp_values[mineral] = value

                mineral_total = list(temp_values.values())
                if np.isclose(np.sum(mineral_total), 100.000000) == True:
                    for mineral, value in temp_values.items():
                        data_amounts[mineral].append(value)
                    n += 1
                    condition = True

        # print("RESULTS (data_amounts):")
        # for key, values in data_amounts.items():
        #     print(key, len(values))
        #
        return data_amounts
    #
    ###########################
    ## SEQUENCE STRATIGRAPHY ##
    ###########################
    #
    def real_sequences(self, var_unit):
        ################################################################################################################
        ## CLEANING
        for category in ["Label", "Button", "Entry", "Radiobutton", "Option Menu"]:
            for gui_element in self.gui_elements["Static"][category]:
                gui_element.grid_remove()
            for gui_element in self.gui_elements["Temporary"][category]:
                gui_element.grid_remove()
        #
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Rockbuilder Static"][category]:
                gui_element.grid_remove()
        #
        for category in ["Label", "Button", "Entry", "Radiobutton"]:
            for gui_element in self.gui_elements["Rockbuilder Temporary"][category]:
                gui_element.grid_remove()
        #
        for key, gui_items in self.gui_elements["Rockbuilder Temporary"].items():
            if len(gui_items) > 0:
                if key not in ["Canvas"]:
                    if type(gui_items) == list:
                        for gui_item in gui_items:
                            gui_item.grid_remove()
                        gui_items.clear()
                else:
                    for key_2, gui_item in gui_items.items():
                        gui_item.get_tk_widget().grid_remove()
        #
        if "Petrology" in self.gui_elements["Temporary"]["Canvas"]:
            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                          "Rock Chemistry Scatter"]
            for category in categories:
                if category in self.gui_elements["Temporary"]["Axis"]:
                    for gui_axes in self.gui_elements["Temporary"]["Axis"][category]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("off")
                            gui_axis.set_visible(False)
            #
            self.gui_elements["Temporary"]["Canvas"]["Petrology"].draw()
        #
        ################################################################################################################
        ## INITIALIZATION
        self.gui_variables["Radiobutton"]["Diagram Type"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Diagram Type"].set(0)
        self.gui_variables["Radiobutton"]["Subunit"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Subunit"].set(0)
        self.gui_variables["Option Menu"]["Lithological Focus"] = tk.StringVar()
        self.gui_variables["Option Menu"]["Lithological Focus"].set("Select Rock")
        #
        ## Labels
        lbl_01 = SimpleElements(
            parent=self.parent, row_id=5, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Subsurface Simulation", font_option="sans 14 bold", relief=tk.FLAT)
        lbl_02 = SimpleElements(
            parent=self.parent, row_id=7, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text=var_unit, font_option="sans 12 bold", relief=tk.FLAT)
        lbl_03 = SimpleElements(
            parent=self.parent, row_id=10, column_id=0, n_rows=6, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].extend([lbl_01, lbl_02, lbl_03])
        #
        ## Radiobuttons
        rb_03a = SimpleElements(
            parent=self.parent, row_id=10, column_id=14, n_rows=2, n_columns=16,
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Well-Log Analysis", var_rb=self.gui_variables["Radiobutton"]["Diagram Type"], value_rb=0,
            color_bg=self.colors_gebpy["Background"])
        rb_03b = SimpleElements(
            parent=self.parent, row_id=12, column_id=14, n_rows=2, n_columns=16,
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Histogram", var_rb=self.gui_variables["Radiobutton"]["Diagram Type"], value_rb=1,
            color_bg=self.colors_gebpy["Background"])
        rb_03c = SimpleElements(
            parent=self.parent, row_id=14, column_id=14, n_rows=2, n_columns=16,
            bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Scatter Plots", var_rb=self.gui_variables["Radiobutton"]["Diagram Type"], value_rb=2,
            color_bg=self.colors_gebpy["Background"])
        #
        self.gui_elements["Static"]["Radiobutton"].extend([rb_03a, rb_03b, rb_03c])
        #
        ## ADDITIONAL FEATURES
        self.build_unit(var_unit=var_unit)

    def build_unit(self, var_unit):
        if var_unit == "Muschelkalk":
            ## LABELS
            lbl_04 = SimpleElements(
                parent=self.parent, row_id=17, column_id=0, n_rows=8, n_columns=14, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Stratigraphic Focus", font_option="sans 10 bold", relief=tk.FLAT)
            lbl_05 = SimpleElements(
                parent=self.parent, row_id=26, column_id=0, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Lithological Focus", font_option="sans 10 bold", relief=tk.FLAT)
            #
            self.gui_elements["Temporary"]["Label"].extend([lbl_04, lbl_05])
            #
            ## RADIOBUTTONS
            rb_04a = SimpleElements(
                parent=self.parent, row_id=17, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Complete Muschelkalk", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=0,
                color_bg=self.colors_gebpy["Background"])
            rb_04b = SimpleElements(
                parent=self.parent, row_id=19, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Upper Muschelkalk", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=1,
                color_bg=self.colors_gebpy["Background"])
            rb_04c = SimpleElements(
                parent=self.parent, row_id=21, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Medium Muschelkalk", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=2,
                color_bg=self.colors_gebpy["Background"])
            rb_04d = SimpleElements(
                parent=self.parent, row_id=23, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Lower Muschelkalk", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=3,
                color_bg=self.colors_gebpy["Background"])
            #
            self.gui_elements["Temporary"]["Radiobutton"].extend([rb_04a, rb_04b, rb_04c, rb_04d])
            #
            ## OPTION MENUS
            list_rocks = Muschelkalk().export_lithological_keys()
            list_rocks.insert(0, "All Rocks")
            self.gui_variables["Option Menu"]["Lithological Focus"].set(list_rocks[0])
            #
            self.rock_data = {}
            for rock in list_rocks:
                self.rock_data[rock] = {
                    "Physics": {"rho": [], "rho_s": [], "phi": [], "vP": [], "vS": [], "vPvS": [], "K": [], "G": [],
                                "E": [], "v": [], "GR": [], "PE": []},
                    "Mineralogy": {}, "Chemistry": {}}
            #
            self.stratigraphy_data = {"Lithology": [], "Bottom": [], "Top": [], "Thickness": []}
            #
            opt_05a = SimpleElements(
                parent=self.parent, row_id=26, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Option"],
                fg=self.colors_gebpy["Navigation"]).create_option_menu(
                var_opt=self.gui_variables["Option Menu"]["Lithological Focus"],
                var_opt_set=self.gui_variables["Option Menu"]["Lithological Focus"].get(), opt_list=list_rocks,
                active_bg=self.colors_gebpy["Accent"],
                command=lambda var_opt=self.gui_variables["Option Menu"]["Lithological Focus"]:
                self.stratigraphy_change_lithology(var_opt))
            #
            self.gui_elements["Temporary"]["Option Menu"].extend([opt_05a])
            #
            ## BUTTONS
            name_unit = var_unit
            btn_06 = SimpleElements(
                parent=self.parent, row_id=29, column_id=14, n_rows=2, n_columns=16,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                text="Run Simulation", command=lambda var_unit=name_unit: self.generate_stratigraphic_data(var_unit))
            #
            self.gui_elements["Temporary"]["Button"].extend([btn_06])
            #
        elif var_unit == "Zechstein":
            ## LABELS
            lbl_04 = SimpleElements(
                parent=self.parent, row_id=17, column_id=0, n_rows=8, n_columns=14, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Stratigraphic Focus", font_option="sans 10 bold", relief=tk.FLAT)
            lbl_05 = SimpleElements(
                parent=self.parent, row_id=30, column_id=0, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Lithological Focus", font_option="sans 10 bold", relief=tk.FLAT)
            #
            self.gui_elements["Temporary"]["Label"].extend([lbl_04, lbl_05])
            #
            ## RADIOBUTTONS
            rb_04a = SimpleElements(
                parent=self.parent, row_id=17, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Complete Zechstein", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=0,
                color_bg=self.colors_gebpy["Background"])
            rb_04b = SimpleElements(
                parent=self.parent, row_id=19, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Zechstein Z1 (Werra)", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=1,
                color_bg=self.colors_gebpy["Background"])
            rb_04c = SimpleElements(
                parent=self.parent, row_id=21, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Zechstein Z2 (Straßfurt)", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=2,
                color_bg=self.colors_gebpy["Background"])
            rb_04d = SimpleElements(
                parent=self.parent, row_id=23, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Zechstein Z3 (Leine)", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=3,
                color_bg=self.colors_gebpy["Background"])
            rb_04e = SimpleElements(
                parent=self.parent, row_id=25, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Zechstein Z4 (Aller)", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=4,
                color_bg=self.colors_gebpy["Background"])
            rb_04f = SimpleElements(
                parent=self.parent, row_id=27, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Zechstein Z5 (Ohre)", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=5,
                color_bg=self.colors_gebpy["Background"])
            #
            self.gui_elements["Temporary"]["Radiobutton"].extend([rb_04a, rb_04b, rb_04c, rb_04d, rb_04e, rb_04f])
            #
            ## OPTION MENUS
            list_rocks = Zechstein().export_lithological_keys()
            list_rocks.insert(0, "All Rocks")
            self.gui_variables["Option Menu"]["Lithological Focus"].set(list_rocks[0])
            #
            self.rock_data = {}
            for rock in list_rocks:
                self.rock_data[rock] = {
                    "Physics": {"rho": [], "rho_s": [], "phi": [], "vP": [], "vS": [], "vPvS": [], "K": [], "G": [],
                                "E": [], "v": [], "GR": [], "PE": []},
                    "Mineralogy": {}, "Chemistry": {}}
            #
            self.stratigraphy_data = {"Lithology": [], "Bottom": [], "Top": [], "Thickness": []}
            #
            opt_05a = SimpleElements(
                parent=self.parent, row_id=30, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Option"],
                fg=self.colors_gebpy["Navigation"]).create_option_menu(
                var_opt=self.gui_variables["Option Menu"]["Lithological Focus"],
                var_opt_set=self.gui_variables["Option Menu"]["Lithological Focus"].get(), opt_list=list_rocks,
                active_bg=self.colors_gebpy["Accent"],
                command=lambda var_opt=self.gui_variables["Option Menu"]["Lithological Focus"]:
                self.stratigraphy_change_lithology(var_opt))
            #
            self.gui_elements["Temporary"]["Option Menu"].extend([opt_05a])
            #
            ## BUTTONS
            name_unit = var_unit
            btn_06 = SimpleElements(
                parent=self.parent, row_id=33, column_id=14, n_rows=2, n_columns=16,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                text="Run Simulation", command=lambda var_unit=name_unit: self.generate_stratigraphic_data(var_unit))
            #
            self.gui_elements["Temporary"]["Button"].extend([btn_06])
            #
        elif var_unit in ["Buntsandstein", "Upper Buntsandstein", "Medium Buntsandstein", "Lower Buntsandstein"]:
            ## LABELS
            lbl_04 = SimpleElements(
                parent=self.parent, row_id=17, column_id=0, n_rows=8, n_columns=14, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Stratigraphic Focus", font_option="sans 10 bold", relief=tk.FLAT)
            lbl_05 = SimpleElements(
                parent=self.parent, row_id=26, column_id=0, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_label(
                text="Lithological Focus", font_option="sans 10 bold", relief=tk.FLAT)
            #
            self.gui_elements["Temporary"]["Label"].extend([lbl_04, lbl_05])
            #
            ## RADIOBUTTONS
            rb_04a = SimpleElements(
                parent=self.parent, row_id=17, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Complete Buntsandstein", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=0,
                color_bg=self.colors_gebpy["Background"])
            rb_04b = SimpleElements(
                parent=self.parent, row_id=19, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Upper Buntsandstein", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=1,
                color_bg=self.colors_gebpy["Background"])
            rb_04c = SimpleElements(
                parent=self.parent, row_id=21, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Medium Buntsandstein", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=2,
                color_bg=self.colors_gebpy["Background"])
            rb_04d = SimpleElements(
                parent=self.parent, row_id=23, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
                fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Lower Buntsandstein", var_rb=self.gui_variables["Radiobutton"]["Subunit"], value_rb=3,
                color_bg=self.colors_gebpy["Background"])
            #
            self.gui_elements["Temporary"]["Radiobutton"].extend([rb_04a, rb_04b, rb_04c, rb_04d])
            #
            ## OPTION MENUS
            list_rocks = Buntsandstein().export_lithological_keys()
            list_rocks.insert(0, "All Rocks")
            self.gui_variables["Option Menu"]["Lithological Focus"].set(list_rocks[0])
            #
            self.rock_data = {}
            for rock in list_rocks:
                self.rock_data[rock] = {
                    "Physics": {"rho": [], "rho_s": [], "phi": [], "vP": [], "vS": [], "vPvS": [], "K": [], "G": [],
                                "E": [], "v": [], "GR": [], "PE": []},
                    "Mineralogy": {}, "Chemistry": {}}
            #
            self.stratigraphy_data = {"Lithology": [], "Bottom": [], "Top": [], "Thickness": []}
            #
            opt_05a = SimpleElements(
                parent=self.parent, row_id=26, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Option"],
                fg=self.colors_gebpy["Navigation"]).create_option_menu(
                var_opt=self.gui_variables["Option Menu"]["Lithological Focus"],
                var_opt_set=self.gui_variables["Option Menu"]["Lithological Focus"].get(), opt_list=list_rocks,
                active_bg=self.colors_gebpy["Accent"],
                command=lambda var_opt=self.gui_variables["Option Menu"]["Lithological Focus"]:
                self.stratigraphy_change_lithology(var_opt))
            #
            self.gui_elements["Temporary"]["Option Menu"].extend([opt_05a])
            #
            ## BUTTONS
            name_unit = var_unit
            btn_06 = SimpleElements(
                parent=self.parent, row_id=29, column_id=14, n_rows=2, n_columns=16,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                text="Run Simulation", command=lambda var_unit=name_unit: self.generate_stratigraphic_data(var_unit))

            self.gui_elements["Temporary"]["Button"].extend([btn_06])

        ## INITIALIZATION
        self.generate_stratigraphic_data(var_unit=var_unit)
        self.stratigraphy_change_lithology(var_opt=self.gui_variables["Option Menu"]["Lithological Focus"].get())
        self.show_well_log_diagram()

    def generate_stratigraphic_data(self, var_unit):
        thickness_complete = rd.randrange(900, 1500, 100)
        #
        if var_unit == "Muschelkalk":
            thickness_muschelkalk_oberer_random = int(rd.uniform(0.3, 0.4) * thickness_complete)
            thickness_muschelkalk_mittlerer_random = int(rd.uniform(0.3, 0.4) * thickness_complete)
            thickness_muschelkalk_unterer_random = int(thickness_complete - thickness_muschelkalk_oberer_random
                                                       - thickness_muschelkalk_mittlerer_random)
            #
            data_muschelkalk_oberer = Muschelkalk(
                actual_thickness=0, resolution=10).create_muschelkalk_oberer(
                top_unit=0, thickness_unit=thickness_muschelkalk_oberer_random)
            data_muschelkalk_mittlerer = Muschelkalk(
                actual_thickness=0, resolution=10).create_muschelkalk_mittlerer(
                top_unit=thickness_muschelkalk_oberer_random, thickness_unit=thickness_muschelkalk_mittlerer_random)
            data_muschelkalk_unterer = Muschelkalk(
                actual_thickness=0).create_muschelkalk_unterer(
                top_unit=thickness_muschelkalk_oberer_random+thickness_muschelkalk_mittlerer_random,
                thickness_unit=thickness_muschelkalk_unterer_random)
            #
            data_units = data_muschelkalk_oberer + data_muschelkalk_mittlerer + data_muschelkalk_unterer
            #
        elif var_unit == "Zechstein":
            thickness_z5_random = int(rd.uniform(0.02, 0.06)*thickness_complete)
            thickness_z4_random = int(rd.uniform(0.18, 0.22)*thickness_complete)
            thickness_z3_random = int(rd.uniform(0.31, 0.35)*thickness_complete)
            thickness_z2_random = int(rd.uniform(0.30, 0.34)*thickness_complete)
            thickness_z1_random = int(thickness_complete - thickness_z5_random - thickness_z4_random -
                                      thickness_z3_random - thickness_z2_random)
            #
            data_z5 = Zechstein(actual_thickness=0).create_zechstein_z5(
                top_z=0, thickness_z5=thickness_z5_random)   # Ohre
            data_z4 = Zechstein(actual_thickness=0).create_zechstein_z4(
                top_z=thickness_z5_random, thickness_z4=thickness_z4_random)   # Aller
            data_z3 = Zechstein(actual_thickness=0).create_zechstein_z3(
                top_z=thickness_z5_random+thickness_z4_random, thickness_z3=thickness_z3_random)   # Leine
            data_z2 = Zechstein(actual_thickness=0).create_zechstein_z2(
                top_z=thickness_z5_random+thickness_z4_random+thickness_z3_random,
                thickness_z2=thickness_z2_random)   # Straßfurt
            data_z1 = Zechstein(actual_thickness=300).create_zechstein_z1(
                top_z=thickness_z5_random+thickness_z4_random+thickness_z3_random+thickness_z2_random,
                thickness_z1=thickness_z1_random)   # Werra
            #
            data_units = data_z5 + data_z4 + data_z3 + data_z2 + data_z1
            #
        elif var_unit == "Buntsandstein":
            thickness_buntsandstein_upper = int(rd.uniform(0.24, 0.28)*thickness_complete)
            thickness_buntsandstein_medium = int(rd.uniform(0.28, 0.32)*thickness_complete)
            thickness_buntsandstein_lower = int(thickness_complete - thickness_buntsandstein_upper -
                                                thickness_buntsandstein_medium)
            #
            data_buntsandstein_upper = Buntsandstein(actual_thickness=0).create_buntsandstein_upper(
                top_unit=0, thickness_unit=thickness_buntsandstein_upper)
            data_buntsandstein_medium = Buntsandstein(actual_thickness=0).create_buntsandstein_medium(
                top_unit=thickness_buntsandstein_upper, thickness_unit=thickness_buntsandstein_medium)
            data_buntsandstein_lower = Buntsandstein(
                actual_thickness=0).create_buntsandstein_lower(
                top_unit=thickness_buntsandstein_upper + thickness_buntsandstein_medium,
                thickness_unit=thickness_buntsandstein_lower)
            #
            data_units = data_buntsandstein_upper + data_buntsandstein_medium + data_buntsandstein_lower
        #
        elif var_unit in ["Upper Buntsandstein", "Medium Buntsandstein", "Lower Buntsandstein"]:
            thickness_buntsandstein = int(thickness_complete)
            #
            if var_unit == "Upper Buntsandstein":
                data_buntsandstein = Buntsandstein(actual_thickness=0).create_buntsandstein_upper(
                    top_unit=0, thickness_unit=thickness_buntsandstein)
            elif var_unit == "Medium Buntsandstein":
                data_buntsandstein = Buntsandstein(actual_thickness=0).create_buntsandstein_medium(
                    top_unit=0, thickness_unit=thickness_buntsandstein)
            elif var_unit == "Lower Buntsandstein":
                data_buntsandstein = Buntsandstein(actual_thickness=0).create_buntsandstein_lower(
                    top_unit=0, thickness_unit=thickness_buntsandstein)
            #
            data_units = data_buntsandstein
        #
        self.unit_sections = {}
        #
        n = 0
        for index, item in enumerate(data_units):
            for key, subitem in item.items():
                var_rock = subitem["rock"]
                var_minerals_list = list(subitem["mineralogy"].keys())
                var_elements_list = list(subitem["chemistry"].keys())
                #
                if var_rock != "All Rocks":
                    if var_rock == "Sandstone":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "tan"}
                    elif var_rock in ["Shale", "Mudstone"]:
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "olivedrab"}
                    elif var_rock in ["Granite", "Gabbro", "Diorite"]:
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "darkorange"}
                    elif var_rock == "Kupferschiefer":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "gray"}
                    elif var_rock in ["limestone", "Limestone"]:
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "skyblue"}
                    elif var_rock == "Anhydrite":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "orchid"}
                    elif var_rock == "Dolostone":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "lightcyan"}
                    elif var_rock == "Rock Salt":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "lavender"}
                    elif var_rock == "Potash":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "yellowgreen"}
                    elif var_rock == "Marl":
                        if var_rock not in self.unit_sections:
                            self.unit_sections[var_rock] = {"Intervals": [], "Color": "moccasin"}
                #
                var_bottom = float(key)
                self.stratigraphy_data["Lithology"].append(var_rock)
                if n > 0:
                    var_top = self.stratigraphy_data["Bottom"][-1]
                    self.stratigraphy_data["Top"].append(var_top)
                else:
                    var_top = 0
                    self.stratigraphy_data["Top"].append(var_top)
                self.stratigraphy_data["Bottom"].append(var_bottom)
                self.stratigraphy_data["Thickness"].append(var_bottom - var_top)
                #
                self.unit_sections[var_rock]["Intervals"].append([var_top, var_bottom])
                #
                n += 1
                #
                ## Physics
                try:
                    self.rock_data[var_rock]["Physics"]["rho"].extend(subitem["rho"])
                    self.rock_data[var_rock]["Physics"]["rho_s"].extend(subitem["rho_s"])
                    self.rock_data[var_rock]["Physics"]["phi"].extend(subitem["phi"])
                    self.rock_data[var_rock]["Physics"]["vP"].extend(subitem["vP"])
                    self.rock_data[var_rock]["Physics"]["vS"].extend(subitem["vS"])
                    self.rock_data[var_rock]["Physics"]["vPvS"].extend(subitem["vP/vS"])
                    self.rock_data[var_rock]["Physics"]["K"].extend(subitem["K"])
                    self.rock_data[var_rock]["Physics"]["G"].extend(subitem["G"])
                    self.rock_data[var_rock]["Physics"]["E"].extend(subitem["E"])
                    self.rock_data[var_rock]["Physics"]["v"].extend(subitem["nu"])
                    self.rock_data[var_rock]["Physics"]["GR"].extend(subitem["GR"])
                    self.rock_data[var_rock]["Physics"]["PE"].extend(subitem["PE"])
                    #
                    self.rock_data["All Rocks"]["Physics"]["rho"].extend(subitem["rho"])
                    self.rock_data["All Rocks"]["Physics"]["rho_s"].extend(subitem["rho_s"])
                    self.rock_data["All Rocks"]["Physics"]["phi"].extend(subitem["phi"])
                    self.rock_data["All Rocks"]["Physics"]["vP"].extend(subitem["vP"])
                    self.rock_data["All Rocks"]["Physics"]["vS"].extend(subitem["vS"])
                    self.rock_data["All Rocks"]["Physics"]["vPvS"].extend(subitem["vP/vS"])
                    self.rock_data["All Rocks"]["Physics"]["K"].extend(subitem["K"])
                    self.rock_data["All Rocks"]["Physics"]["G"].extend(subitem["G"])
                    self.rock_data["All Rocks"]["Physics"]["E"].extend(subitem["E"])
                    self.rock_data["All Rocks"]["Physics"]["v"].extend(subitem["nu"])
                    self.rock_data["All Rocks"]["Physics"]["GR"].extend(subitem["GR"])
                    self.rock_data["All Rocks"]["Physics"]["PE"].extend(subitem["PE"])
                except:
                    self.rock_data[var_rock]["Physics"]["rho"].append(subitem["rho"])
                    self.rock_data[var_rock]["Physics"]["rho_s"].append(subitem["rho_s"])
                    self.rock_data[var_rock]["Physics"]["phi"].append(subitem["phi"])
                    self.rock_data[var_rock]["Physics"]["vP"].append(subitem["vP"])
                    self.rock_data[var_rock]["Physics"]["vS"].append(subitem["vS"])
                    self.rock_data[var_rock]["Physics"]["vPvS"].append(subitem["vP/vS"])
                    self.rock_data[var_rock]["Physics"]["K"].append(subitem["K"])
                    self.rock_data[var_rock]["Physics"]["G"].append(subitem["G"])
                    self.rock_data[var_rock]["Physics"]["E"].append(subitem["E"])
                    self.rock_data[var_rock]["Physics"]["v"].append(subitem["nu"])
                    self.rock_data[var_rock]["Physics"]["GR"].append(subitem["GR"])
                    self.rock_data[var_rock]["Physics"]["PE"].append(subitem["PE"])
                    #
                    self.rock_data["All Rocks"]["Physics"]["rho"].append(subitem["rho"])
                    self.rock_data["All Rocks"]["Physics"]["rho_s"].append(subitem["rho_s"])
                    self.rock_data["All Rocks"]["Physics"]["phi"].append(subitem["phi"])
                    self.rock_data["All Rocks"]["Physics"]["vP"].append(subitem["vP"])
                    self.rock_data["All Rocks"]["Physics"]["vS"].append(subitem["vS"])
                    self.rock_data["All Rocks"]["Physics"]["vPvS"].append(subitem["vP/vS"])
                    self.rock_data["All Rocks"]["Physics"]["K"].append(subitem["K"])
                    self.rock_data["All Rocks"]["Physics"]["G"].append(subitem["G"])
                    self.rock_data["All Rocks"]["Physics"]["E"].append(subitem["E"])
                    self.rock_data["All Rocks"]["Physics"]["v"].append(subitem["nu"])
                    self.rock_data["All Rocks"]["Physics"]["GR"].append(subitem["GR"])
                    self.rock_data["All Rocks"]["Physics"]["PE"].append(subitem["PE"])
                    #
                ## Mineralogy
                for mineral in var_minerals_list:
                    if mineral not in self.rock_data[var_rock]["Mineralogy"]:
                        self.rock_data[var_rock]["Mineralogy"][mineral] = []
                        self.rock_data["All Rocks"]["Mineralogy"][mineral] = []
                        #
                        try:
                            self.rock_data[var_rock]["Mineralogy"][mineral].extend(subitem["mineralogy"][mineral])
                            self.rock_data["All Rocks"]["Mineralogy"][mineral].extend(subitem["mineralogy"][mineral])
                        except:
                            self.rock_data[var_rock]["Mineralogy"][mineral].append(subitem["mineralogy"][mineral])
                            self.rock_data["All Rocks"]["Mineralogy"][mineral].append(subitem["mineralogy"][mineral])
                    else:
                        try:
                            self.rock_data[var_rock]["Mineralogy"][mineral].extend(subitem["mineralogy"][mineral])
                            self.rock_data["All Rocks"]["Mineralogy"][mineral].extend(subitem["mineralogy"][mineral])
                        except:
                            self.rock_data[var_rock]["Mineralogy"][mineral].append(subitem["mineralogy"][mineral])
                            self.rock_data["All Rocks"]["Mineralogy"][mineral].append(subitem["mineralogy"][mineral])
                ## Chemistry
                for element in var_elements_list:
                    if element not in self.rock_data[var_rock]["Chemistry"]:
                        self.rock_data[var_rock]["Chemistry"][element] = []
                        self.rock_data["All Rocks"]["Chemistry"][element] = []
                        #
                        try:
                            self.rock_data[var_rock]["Chemistry"][element].extend(subitem["chemistry"][element])
                            self.rock_data["All Rocks"]["Chemistry"][element].extend(subitem["chemistry"][element])
                        except:
                            self.rock_data[var_rock]["Chemistry"][element].append(subitem["chemistry"][element])
                            self.rock_data["All Rocks"]["Chemistry"][element].append(subitem["chemistry"][element])
                    else:
                        try:
                            self.rock_data[var_rock]["Chemistry"][element].extend(subitem["chemistry"][element])
                            self.rock_data["All Rocks"]["Chemistry"][element].extend(subitem["chemistry"][element])
                        except:
                            self.rock_data[var_rock]["Chemistry"][element].append(subitem["chemistry"][element])
                            self.rock_data["All Rocks"]["Chemistry"][element].append(subitem["chemistry"][element])
        #
        ## TREE VIEW
        categories = ["rho (kg/m\u00B3)", "phi (%)", "vP (m/s)", "vS (m/s)", "vP/vS (1)", "K (GPa)", "G (GPa)",
                      "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
        categories_short = ["rho", "phi", "vP", "vS", "vPvS", "K", "G", "E", "v", "GR", "PE"]
        #
        list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
        list_width = list(75*np.ones(len(list_categories)))
        list_width = [int(item) for item in list_width]
        list_width[0] = 90
        list_width[-1] = 150
        #
        self.tv_strat_results = SimpleElements(
            parent=self.parent, row_id=0, column_id=35, n_rows=30, n_columns=45,
            fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
            n_categories=len(list_categories), text_n=list_categories,
            width_n=list_width, individual=True)
        #
        scb_v = ttk.Scrollbar(self.parent, orient="vertical")
        scb_h = ttk.Scrollbar(self.parent, orient="horizontal")
        self.tv_strat_results.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
        scb_v.config(command=self.tv_strat_results.yview)
        scb_h.config(command=self.tv_strat_results.xview)
        scb_v.grid(row=0, column=35 + 45, rowspan=30, columnspan=1, sticky="ns")
        scb_h.grid(row=30, column=35, rowspan=1, columnspan=45, sticky="ew")
        #
    def stratigraphy_change_lithology(self, var_opt):
        if len(self.tv_strat_results.get_children()) > 0:
            for item in self.tv_strat_results.get_children():
                self.tv_strat_results.delete(item)
        #
        categories = ["rho (kg/m\u00B3)", "phi (%)", "vP (m/s)", "vS (m/s)", "vP/vS (1)", "K (GPa)", "G (GPa)",
                      "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
        categories_short = ["rho", "phi", "vP", "vS", "vPvS", "K", "G", "E", "v", "GR", "PE"]
        #
        for index, category in enumerate(categories):
            entries = [category]
            #
            n_digits = 2
            if categories_short[index] == "phi":
                var_factor = 100
            else:
                var_factor = 1
                #
            var_entr_min = round(
                var_factor*min(self.rock_data[var_opt]["Physics"][categories_short[index]]), n_digits)
            var_entr_max = round(
                var_factor*max(self.rock_data[var_opt]["Physics"][categories_short[index]]), n_digits)
            var_entr_mean = round(
                var_factor*np.mean(self.rock_data[var_opt]["Physics"][categories_short[index]]), n_digits)
            var_entr_error = round(
                var_factor*np.std(self.rock_data[var_opt]["Physics"][categories_short[index]], ddof=1), n_digits)
            #
            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            #
            self.tv_strat_results.insert("", tk.END, values=entries)
        #
        entries = ["-", "-", "-", "-", "-"]
        self.tv_strat_results.insert("", tk.END, values=entries)
        #
        for mineral, dataset in self.rock_data[var_opt]["Mineralogy"].items():
            entries = [str(mineral)+str(" (%)")]
            #
            n_digits = 2
            var_factor = 100
            #
            var_entr_min = round(var_factor*min(dataset), n_digits)
            var_entr_max = round(var_factor*max(dataset), n_digits)
            var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
            var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)
            #
            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            #
            self.tv_strat_results.insert("", tk.END, values=entries)
        #
        entries = ["-", "-", "-", "-", "-"]
        self.tv_strat_results.insert("", tk.END, values=entries)
        #
        for element, dataset in self.rock_data[var_opt]["Chemistry"].items():
            entries = [str(element)+str(" (%)")]
            #
            n_digits = 2
            var_factor = 100
            #
            var_entr_min = round(var_factor*min(dataset), n_digits)
            var_entr_max = round(var_factor*max(dataset), n_digits)
            var_entr_mean = round(var_factor*np.mean(dataset), n_digits)
            var_entr_error = round(var_factor*np.std(dataset, ddof=1), n_digits)
            #
            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            #
            self.tv_strat_results.insert("", tk.END, values=entries)
    #
    def show_well_log_diagram(self):
        self.canvas = None
        max_thickness = max(self.stratigraphy_data["Bottom"])
        if max_thickness <= 100:
            step_depth = 10
        elif 100 < max_thickness <= 500:
            step_depth = 50
        elif 500 < max_thickness <= 1500:
            step_depth = 100
        elif max_thickness > 1500:
            step_depth = 200
        #
        self.fig, (self.ax1, self.ax2, self.ax3, self.ax4, self.ax5) = plt.subplots(
            1, 5, sharey="row", gridspec_kw={"wspace": 0.25}, figsize=(12, 24),
            facecolor=self.colors_gebpy["Background"])
        self.fig.subplots_adjust(wspace=0.25)
        # 1
        self.ax1.plot(self.rock_data["All Rocks"]["Physics"]["GR"], self.stratigraphy_data["Top"], color="#00549F",
                      linewidth=2)
        self.ax1.set_xlabel("GR [API]")
        self.ax1.set_ylabel("Depth [m]")
        if max(self.rock_data["All Rocks"]["Physics"]["GR"]) > 250:
            self.ax1.set_xlim(-0.5, max(self.rock_data["All Rocks"]["Physics"]["GR"]))
            self.ax1.set_xticks(np.arange(-0.5, max(self.rock_data["All Rocks"]["Physics"]["GR"]) + 200, 200))
            self.ax1.set_xscale("symlog")
        elif 50 < max(self.rock_data["All Rocks"]["Physics"]["GR"]) < 100:
            self.ax1.set_xlim(-0.5, max(self.rock_data["All Rocks"]["Physics"]["GR"]))
            self.ax1.set_xticks(np.arange(-0.5, max(self.rock_data["All Rocks"]["Physics"]["GR"]) + 20, 20))
            self.ax1.set_xscale("symlog")
        else:
            self.ax1.set_xlim(-1, max(self.rock_data["All Rocks"]["Physics"]["GR"]))
            self.ax1.set_xticks(np.arange(0, max(self.rock_data["All Rocks"]["Physics"]["GR"]) + 10, 10))
        #
        self.ax1.set_ylim(0, max_thickness)
        self.ax1.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax1.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 2
        vP_edit = np.array(self.rock_data["All Rocks"]["Physics"]["vP"])/1000
        vS_edit = np.array(self.rock_data["All Rocks"]["Physics"]["vS"])/1000
        self.ax2.plot(vP_edit, self.stratigraphy_data["Top"], color="#00549F", linewidth=2)
        self.ax2.set_xlabel("$v_P$ [km/s]")
        self.ax2.set_xlim(0, max(vP_edit))
        self.ax2.set_xticks(np.arange(0, max(vP_edit)+2.0, 2.0))
        self.ax2.xaxis.label.set_color("#00549F")
        self.ax2.set_ylim(0, max_thickness)
        self.ax2.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax2.grid(color="grey", linestyle="dashed")
        self.ax2_2 = self.ax2.twiny()
        self.ax2_2.plot(vS_edit, self.stratigraphy_data["Top"], color="#CC071E", linewidth=2)
        self.ax2_2.set_xlabel("$v_S$ [km/s]")
        self.ax2_2.set_xlim(0, max(vP_edit))
        self.ax2_2.set_xticks(np.arange(0, max(vP_edit)+2.0, 2.0))
        self.ax2_2.minorticks_on()
        self.ax2_2.xaxis.label.set_color("#CC071E")
        self.ax2_2.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # # 3
        phi_edit = np.array(self.rock_data["All Rocks"]["Physics"]["phi"])*100
        self.ax3.plot(np.array(
            self.rock_data["All Rocks"]["Physics"]["rho"])/1000, self.stratigraphy_data["Top"], color="#57AB27",
                      linewidth=2)
        self.ax3.set_xlabel("$\\varrho$ [g/cm$^3$]")
        self.ax3.set_xlim(1.7, 3.2)
        self.ax3.set_xticks(np.around(np.linspace(1.7, 3.2, 4, endpoint=True), decimals=1))
        self.ax3.xaxis.label.set_color("#57AB27")
        self.ax3.set_ylim(0, max_thickness)
        self.ax3.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax3.grid(color="grey", linestyle="dashed")
        self.ax3_2 = self.ax3.twiny()
        self.ax3_2.plot(phi_edit, self.stratigraphy_data["Top"], color="#00549F", linewidth=2)
        self.ax3_2.set_xlabel("$\\varphi$ [1]")
        self.ax3_2.set_xlim(60, -30)
        self.ax3_2.set_xticks(np.around(np.linspace(60, -30, 6, endpoint=True), decimals=0))
        self.ax3_2.minorticks_on()
        self.ax3_2.xaxis.label.set_color("#00549F")
        self.ax3_2.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # # 4
        self.ax4.plot(
            self.rock_data["All Rocks"]["Physics"]["PE"], self.stratigraphy_data["Top"], color="#00549F", linewidth=2)
        self.ax4.set_xlabel("PE [barns/electron]")
        #
        if max(self.rock_data["All Rocks"]["Physics"]["PE"]) > 10:
            self.ax4.set_xlim(-0.5, max(self.rock_data["All Rocks"]["Physics"]["PE"]))
            self.ax4.set_xticks(np.arange(-0.5, len(str(int(max(self.rock_data["All Rocks"]["Physics"]["PE"])))), 3))
            self.ax4.set_xscale("symlog")
        else:
            self.ax4.set_xlim(0, max(self.rock_data["All Rocks"]["Physics"]["PE"]))
            self.ax4.set_xticks(np.arange(0, max(self.rock_data["All Rocks"]["Physics"]["PE"]) + 1, 1))
        #
        self.ax4.set_ylim(0, max_thickness)
        self.ax4.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax4.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # # 5
        if self.unit_sections == None:
            n_units = []
            units_sorted = []
            for rock in self.list_rocks_short:
                units_sorted.append([rock])
                n_units.append(sum(self.results_plot[rock]["thickness"]))
                for index, value in enumerate(self.results_plot[rock]["thickness"], start=0):
                    if self.results_plot[rock]["Top"][index] != self.results_plot[rock]["Bottom"][index]:
                        units_sorted[-1].append(
                            [self.results_plot[rock]["Top"][index], self.results_plot[rock]["Bottom"][index]])
                    else:
                        units_sorted[-1].append(
                            [self.results_plot[rock]["Bottom"][index-1], self.results_plot[rock]["Bottom"][index]])
                if rock == "Sandstone":
                    units_sorted[-1].append("tan")
                elif rock in ["Shale", "Mudstone"]:
                    units_sorted[-1].append("olivedrab")
                elif rock in ["Granite", "Gabbro", "Diorite"]:
                    units_sorted[-1].append("darkorange")
                elif rock == "Kupferschiefer":
                    units_sorted[-1].append("gray")
                elif rock in ["limestone", "Limestone"]:
                    units_sorted[-1].append("skyblue")
                elif rock == "Anhydrite":
                    units_sorted[-1].append("orchid")
                elif rock == "Dolomite":
                    units_sorted[-1].append("lightcyan")
                elif rock == "Rock Salt":
                    units_sorted[-1].append("lavender")
                elif rock == "Potash":
                    units_sorted[-1].append("yellowgreen")
                elif rock == "Marl":
                    units_sorted[-1].append("moccasin")
            legend_lithology = []
            for i in range(len(units_sorted)):
                legend_lithology.append(mpatches.Patch(facecolor=units_sorted[i][-1], edgecolor="black", hatch="",
                                                       label=units_sorted[i][0]))
            for i in range(len(n_units)):
                for j in range(1, len(units_sorted[i])-1):
                    self.ax5.hist(x=np.linspace(units_sorted[i][j][0], units_sorted[i][j][1]), bins=len(n_units),
                                  color=units_sorted[i][-1], orientation="horizontal")
        else:
            n_units = len(self.unit_sections)
            legend_lithology = []
            for key, value in self.unit_sections.items():
                legend_lithology.append(
                    mpatches.Patch(facecolor=value["Color"], edgecolor="black", hatch="", label=key))
            for key, value in self.unit_sections.items():
                for interval in value["Intervals"]:
                    self.ax5.hist(x=np.linspace(interval[0], interval[1]), bins=n_units,
                                  color=value["Color"], orientation="horizontal")
        self.ax5.set_xlabel("Lithology")
        self.ax5.set_xlim(0, 5)
        self.ax5.set_xticks([])
        self.ax5.set_ylim(0, max_thickness)
        self.ax5.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax5.margins(0.3, 0.0)
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        self.ax5.legend(handles=legend_lithology, loc="lower left", bbox_to_anchor=(0, -0.125), shadow=False, ncol=2,
                        prop={'size': 7}, frameon=False)
        # #plt.tight_layout()
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.parent)
        self.canvas.get_tk_widget().grid(row=0, column=81, rowspan=self.n_rows, columnspan=self.n_columns - 80,
                                         sticky="nesw")
        self.canvas.draw()
#
if __name__ == "__main__":
    root = tk.Tk()
    #
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    #
    GebPyGUI(parent=root, var_screen_width=screen_width, var_screen_height=screen_height)
    #
    root.mainloop()