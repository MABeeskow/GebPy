#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gebpy_app.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.05.2023

#-----------------------------------------------

## MODULES
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
# Sequence Stratigraphy
from modules.series import Muschelkalk, Zechstein, Buntsandstein
from modules.petrophysics import SeismicVelocities

## GUI
class GebPyGUI(tk.Frame):
    #
    def __init__(self, parent, var_screen_width, var_screen_height):
        tk.Frame.__init__(self, parent)
        #
        var_screen_width = var_screen_width
        var_screen_height = var_screen_height
        #
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
            gebpy_logo = tk.PhotoImage(file="documents/readme_images/GebPy_Logo_new.png")
            gebpy_logo = gebpy_logo.subsample(5, 5)
            img = tk.Label(self.parent, image=gebpy_logo, bg=self.colors_gebpy["Navigation"])
            img.image = gebpy_logo
            img.grid(row=0, column=0, rowspan=5, columnspan=32, sticky="nesw")
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
                    "Gibbsite", "Au(III)-Oxide", "Brucite"]
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
                    "Alkaline Feldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline"]
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
                    "Siliciclastic": ["Sandstone", "Shale", "Mudstone", "Conglomerate", "Greywacke (Huckenholz)"],
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
                    "Fe-Ore": ["Itabirite", "Compact Hematite", "Friable Hematite", "Goethite Hematite",
                               "Al-rich Itabirite", "Compact Quartz Itabirite", "Friable Quartz Itabirite",
                               "Goethite Itabirite"]}
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
            label="Rock Comparison", command=self.rock_comparison)
        petrology_menu.add_command(
            label="Rock Builder", command=self.rock_builder)
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
        triassic_units = ["Keuper", "Muschelkalk", "Buntsandstein"]
        for unit in triassic_units:
            sub_triassic.add_command(label=unit, command=lambda var_unit=unit: self.real_sequences(var_unit))
        #
        # Real Sequences
        stratigraphy_menu.add_cascade(
            label="Triassic Units",
            menu=sub_triassic)
        stratigraphy_menu.add_cascade(
            label="Permian Units",
            menu=sub_permian)
        #
        stratigraphy_menu.add_separator()
        stratigraphy_menu.add_command(
            label="Create Sequences")
        #
        menubar.add_cascade(
            label="Stratigraphy",
            menu=stratigraphy_menu)
        #
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
        #
        ## Buttons
        btn_quit = SimpleElements(
            parent=self.parent, row_id=n_rows - 3, column_id=16, n_rows=3, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Quit GebPy", command=self.parent.quit)
        btn_restart = SimpleElements(
            parent=self.parent, row_id=n_rows - 3, column_id=1, n_rows=3, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Restart GebPy", command=self.restart_gebpy)
        #
        #self.gui_elements["Static"]["Button"].extend([btn_quit, btn_restart])
    #
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
        categories_short = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]
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
        rb_trace_with = SimpleElements(
            parent=self.parent, row_id=13, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="With Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
        #
        self.gui_elements["Static"]["Radiobutton"].extend([rb_trace_without, rb_trace_with])
        #
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
    #
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
                        labels = [["M (kg/mol)", "V (A$^3$/mol", "rho (kg/m$^3$"], ["vP (m/s)", "vS (m/s)", "vP/vS (1)"],
                                  ["GR (API)", "PE (barns/e$^-$)", "nu (1)"]]
                        #
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
                                #
                                x_min = round(x_min - 0.1*delta_x, n_digits)
                                x_max = round(x_max + 0.1*delta_x, n_digits)
                                #
                                if key != "nu":
                                    if x_min < 0:
                                        x_min = 0
                                #
                                ax_mp_histo[i][j].set_xlim(left=x_min, right=x_max)
                                ax_mp_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                                ax_mp_histo[i][j].set_xticks(np.around(
                                    np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits))
                                ax_mp_histo[i][j].set_yticks(np.around(
                                    np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 1))
                                ax_mp_histo[i][j].xaxis.set_tick_params(labelsize=8)
                                ax_mp_histo[i][j].yaxis.set_tick_params(labelsize=8)
                                ax_mp_histo[i][j].set_xlabel(labels[i][j], fontsize=8)
                                ax_mp_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=8)
                                ax_mp_histo[i][j].grid(True)
                                ax_mp_histo[i][j].set_axisbelow(True)
                            #
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
                        labels = [["M (kg/mol)", "V (A$^3$/mol", "rho (kg/m$^3$"], ["vP (m/s)", "vS (m/s)", "vP/vS (1)"],
                                  ["GR (API)", "PE (barns/e$^-$)", "nu (1)"]]
                        #
                        dataset_x = self.data_mineral["rho"]
                        #
                        for i, subcategories in enumerate(categories):
                            for j, key in enumerate(subcategories):
                                dataset_y = self.data_mineral[key]
                                ax_mp_scatter[i][j].scatter(
                                    dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black",
                                    alpha=0.5)
                                #
                                x_min = min(dataset_x)
                                x_max = max(dataset_x)
                                delta_x = round(x_max - x_min, 4)
                                y_min = min(dataset_y)
                                y_max = max(dataset_y)
                                delta_y = round(y_max - y_min, 4)
                                #
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
                                #
                                x_min = round(x_min - 0.1*delta_x, n_digits_x)
                                x_max = round(x_max + 0.1*delta_x, n_digits_x)
                                y_min = round(y_min - 0.1*delta_y, n_digits_y)
                                y_max = round(y_max + 0.1*delta_y, n_digits_y)
                                #
                                if key != "nu":
                                    if y_min < 0:
                                        y_min = 0
                                #
                                ax_mp_scatter[i][j].set_xlim(left=x_min, right=x_max)
                                ax_mp_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                                ax_mp_scatter[i][j].set_xticks(np.around(
                                    np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits_x))
                                ax_mp_scatter[i][j].set_yticks(np.around(
                                    np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), n_digits_y))
                                ax_mp_scatter[i][j].xaxis.set_tick_params(labelsize=8)
                                ax_mp_scatter[i][j].yaxis.set_tick_params(labelsize=8)
                                ax_mp_scatter[i][j].set_xlabel("Density - kg/m$^3$", fontsize=8)
                                ax_mp_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=8)
                                ax_mp_scatter[i][j].grid(True)
                                ax_mp_scatter[i][j].set_axisbelow(True)
                        #
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
                        ax_mc_histo[i][j].grid(True)
                        ax_mc_histo[i][j].set_axisbelow(True)
                #
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
                            dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black")
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
                        #
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
                        ax_mc_scatter[i][j].grid(True)
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
                lbl_title = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=2, n_columns=45,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mineral Physics", font_option="sans 12 bold", relief=tk.GROOVE)
                lbl_results = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Results", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 9, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Minimum", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 18, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Maximum", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_mean = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 27, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mean", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_error = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 36, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Error", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_diagram_type = SimpleElements(
                    parent=self.parent, row_id=33, column_id=0, n_rows=4, n_columns=14,
                    bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_label(
                    text="Diagram Type", font_option="sans 10 bold", relief=tk.FLAT)
                #
                self.gui_elements["Temporary"]["Label"].extend(
                    [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error, lbl_diagram_type])
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
                categories = [
                    "M\n (kg/mol)", "V\n (\u00C5\u00B3/mol)", "rho\n (kg/m\u00B3)", "vP\n (m/s)", "vS\n (m/s)",
                    "vP/vS\n (1)", "K\n (GPa)", "G\n (GPa)", "E\n (GPa)", "nu\n (1)", "GR\n (API)", "PE\n (barns/e\u207B)"]
                categories_short = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]
                for index, category in enumerate(categories):
                    lbl_category = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=category, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Temporary"]["Label"].append(lbl_category)
                    #
                    ## Entries
                    n_digits = 2
                    var_entr_min = round(min(self.data_mineral[categories_short[index]]), n_digits)
                    var_entr_max = round(max(self.data_mineral[categories_short[index]]), n_digits)
                    var_entr_mean = round(np.mean(self.data_mineral[categories_short[index]]), n_digits)
                    var_entr_error = round(np.std(self.data_mineral[categories_short[index]], ddof=1), n_digits)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 9, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][categories_short[index]],
                        var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 18, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][categories_short[index]],
                        var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 27, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][categories_short[index]],
                        var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 36, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][categories_short[index]],
                        var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                    #
                self.change_rb_diagram()
                #
                # ## TREE VIEW
                # categories = [
                #     "M (kg/mol)", "V (\u00C5\u00B3/mol)", "rho (kg/m\u00B3)", "vP (m/s)", "vS (m/s)", "vP/vS (1)",
                #     "K (GPa)", "G (GPa)", "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
                # list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
                # list_width = list(75*np.ones(len(list_categories)))
                # list_width = [int(item) for item in list_width]
                # list_width[0] = 90
                # list_width[-1] = 150
                # #
                # tv_ma_results = SimpleElements(
                #     parent=self.parent, row_id=0, column_id=start_column, n_rows=20, n_columns=45,
                #     fg=self.colors_gebpy["Black"], bg=self.colors_gebpy["White"]).create_treeview(
                #     n_categories=len(list_categories), text_n=list_categories,
                #     width_n=list_width, individual=True)
                # #
                # scb_v = ttk.Scrollbar(self.parent, orient="vertical")
                # scb_h = ttk.Scrollbar(self.parent, orient="horizontal")
                # tv_ma_results.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
                # scb_v.config(command=tv_ma_results.yview)
                # scb_h.config(command=tv_ma_results.xview)
                # scb_v.grid(row=0, column=start_column + 45, rowspan=20, columnspan=1, sticky="ns")
                # scb_h.grid(row=20, column=start_column, rowspan=1, columnspan=45, sticky="ew")
                # #
                # for index, category in enumerate(categories):
                #     entries = [category]
                #     #
                #     n_digits = 2
                #     var_entr_min = round(min(self.data_mineral[categories_short[index]]), n_digits)
                #     var_entr_max = round(max(self.data_mineral[categories_short[index]]), n_digits)
                #     var_entr_mean = round(np.mean(self.data_mineral[categories_short[index]]), n_digits)
                #     var_entr_error = round(np.std(self.data_mineral[categories_short[index]], ddof=1), n_digits)
                #     #
                #     entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
                #     #
                #     tv_ma_results.insert("", tk.END, values=entries)
                #
            else:
                pass
            #
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
                lbl_title = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=2, n_columns=45,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mineral Chemistry", font_option="sans 12 bold", relief=tk.GROOVE)
                lbl_results = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Results", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 9, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Minimum", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 18, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Maximum", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_mean = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 27, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mean", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_error = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 36, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Error", font_option="sans 10 bold", relief=tk.GROOVE)
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
                    [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error, lbl_diagram_type, lbl_element,
                     lbl_compound, lbl_concentration_setup, lbl_element])
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
                for index, element in enumerate(self.list_elements):
                    lbl_element = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column, n_rows=2, n_columns=9,
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
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 9, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][element], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 18, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][element], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 27, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][element], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 36, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][element], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                    self.gui_elements_sub["Mineral Chemistry"]["Element Concentration"]["Entry"].extend(
                        [entr_min, entr_max, entr_mean, entr_error])
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
                ## Labels
                lbl_title = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=2, n_columns=45,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Synthetic LA-ICP-MS", font_option="sans 12 bold", relief=tk.GROOVE)
                lbl_results = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Results", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_min = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 9, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Minimum", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_max = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 18, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Maximum", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_mean = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 27, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mean", font_option="sans 10 bold", relief=tk.GROOVE)
                lbl_error = SimpleElements(
                    parent=self.parent, row_id=2, column_id=start_column + 36, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Error", font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements["Temporary"]["Label"].extend(
                    [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
                #
                self.list_elements.sort()
                for index, element in enumerate(self.list_elements):
                    lbl_element = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=element, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Temporary"]["Label"].append(lbl_element)
                    #
                    ## Entries
                    var_entr_min = int(min(self.data_mineral["chemistry"][element])*10**6)
                    var_entr_max = int(max(self.data_mineral["chemistry"][element])*10**6)
                    var_entr_mean = int(np.mean(self.data_mineral["chemistry"][element])*10**6)
                    var_entr_error = int(np.std(self.data_mineral["chemistry"][element], ddof=1)*10**6)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 9, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][element], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 18, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][element], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 27, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][element], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 36, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][element], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
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
                    ax_laicpms[0][0].set_ylim(1, 10**9)
                    ax_laicpms[0][0].set_xticks(np.arange(0, 60 + 5, 5))
                    ax_laicpms[0][0].set_yscale("log")
                    ax_laicpms[0][0].grid(True)
                    #
                    ax_laicpms[0][0].grid(which="major", axis="both", linestyle="-")
                    ax_laicpms[0][0].minorticks_on()
                    ax_laicpms[0][0].grid(which="minor", axis="both", linestyle="-", alpha=0.25)
                    #
                    ax_laicpms[0][0].set_axisbelow(True)
                    #
                    ax_laicpms[0][0].legend()
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
    def simulate_srm(self, var_srm="NIST 610 (GeoRem)"):
        data_srm = {}
        if var_srm == "NIST 610 (GeoRem)":
            concentration_data = {
                "Li": 485, "Be": 466, "B": 356, "Na": 102970, "Mg": 465, "Al": 10791, "Si": 327091, "P": 409, "S": 570,
                "Cl": 470, "K": 465, "Ca": 81833, "Sc": 452, "Ti": 460, "V": 442, "Cr": 405, "Mn": 433, "Fe": 457,
                "Co": 405, "Ni": 459, "Cu": 430, "Zn": 456, "Ga": 438, "Ge": 426, "As": 317, "Se": 112, "Br": 33,
                "Rb": 426, "Sr": 516, "Y": 458, "Zr": 437, "Nb": 485, "Mo": 410, "Rh": 1, "Pd": 1, "Ag": 239, "Cd": 259,
                "In": 441, "Sn": 427, "Sb": 405, "Te": 327, "Cs": 357, "Ba": 454, "La": 440, "Ce": 458, "Pr": 443,
                "Nd": 437, "Sm": 453, "Eu": 444, "Gd": 456, "Tb": 440, "Dy": 436, "Ho": 440, "Er": 456, "Tm": 423,
                "Yb": 455, "Lu": 440, "Hf": 421, "Ta": 482, "W": 447, "Re": 50, "Pt": 3, "Au": 23, "Tl": 61, "Pb": 426,
                "Bi": 358, "Th": 457, "U": 462}
            #
            normalized_sensitivity_data = {
                "Li": 2555, "Be": 348, "B": 408, "Na": 5400, "Al": 3825, "Si": 78, "Ti": 462, "Mn": 8282, "Fe": 340,
                "Sn": 3439}
            #
            analytical_sensitivity_data = {
                "Li": 32.96, "Be": 4.49, "B": 5.27, "Na": 69.66, "Al": 49.34, "Si": 1.00, "Ti": 5.95, "Mn": 106.83,
                "Fe": 4.38, "Sn": 44.35}
            #
        data_srm["Concentration"] = concentration_data
        data_srm["Normalized Sensitivity"] = normalized_sensitivity_data
        data_srm["Analytical Sensitivity"] = analytical_sensitivity_data
        #
        return data_srm
        #
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
        if self.data_mineral["mineral"] == "Qz":
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
                #
                mean_bg = np.random.randint(1, 100)
                error_bg = np.random.uniform(0.1, 0.5)*mean_bg
                mean_sig = np.mean(self.data_mineral["chemistry"][element])*total_ppm
                if self.data_mineral["mineral"] == "Qz":
                    if element == "Si":
                        mean_sig_is = np.mean(self.data_mineral["chemistry"][element])*10**6*self.data_mineral[
                            "LA-ICP-MS"][element]
                    else:
                        mean_sig = data_srm["Analytical Sensitivity"][element]*mean_sig_is*np.mean(
                            self.data_mineral["chemistry"][element])/np.mean(self.data_mineral["chemistry"]["Si"])
                # error_sig = np.random.uniform(0, 0.025)*mean_sig
                #
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
                    elif 50 < time_value <= 52:
                        if last_value_sig > 10**4:
                            value = int(amount_end[index_end]*(last_value_sig + np.random.randint(10**2, 10**3) + value_bg))
                        else:
                            value = int(amount_end[index_end]*(last_value_sig + np.random.randint(10**1, 10**2) + value_bg))
                        last_value = value
                        index_end += 1
                    elif 52 < time_value <= 54:
                        if last_value_sig > 10**4:
                            value = int(amount_end2[index_end2]*(last_value + np.random.randint(10**2, 10**5)) + value_bg)
                        else:
                            value = int(amount_end2[index_end2]*(last_value + np.random.randint(10**2, 10**3)) + value_bg)
                        index_end2 += 1
                    elif time_value > 54:
                        value = int(np.random.normal(loc=mean_bg, scale=error_bg, size=1)[0])
                        if value < 0:
                            value = 0
                    #
                    intensity_data[element].append(value)
            #
        return time_data, intensity_data
    #
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
            #
            if self.btn_traces == None:
                self.btn_traces = SimpleElements(
                    parent=self.parent, row_id=15, column_id=16, n_rows=2, n_columns=15,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                    text="Trace Elements",
                    command=lambda var_traces=self.trace_elements_all:
                    self.select_trace_elements(var_traces))
            #
            self.gui_elements["Temporary"]["Button"].append(self.btn_traces)
    #
    def select_trace_elements(self, var_traces):
        #
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
        #
        self.oxidation_states = list(var_traces.keys())
        self.oxidation_states.remove("All")
        self.oxidation_states.sort(reverse=True)
        for index, oxidation_state in enumerate(self.oxidation_states):
            rb_oxidation_state = SimpleElements(
                parent=self.window_trace_elements, row_id=2*index + 2, column_id=0, n_rows=2, n_columns=12,
                bg=self.colors_gebpy["Navigation"], fg=self.colors_gebpy["Background"]).create_radiobutton(
                text="Mainly "+str(oxidation_state), var_rb=self.gui_variables["Radiobutton"]["Oxidation State"],
                value_rb=index, color_bg=self.colors_gebpy["Navigation"],
                command=lambda var_traces=var_traces: self.change_oxidation_state(var_traces))
            #
            self.gui_elements_sub["Trace Elements"]["Static"]["Radiobutton"].append(rb_oxidation_state)
        #
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
        #
        oxidation_state = self.oxidation_states[self.gui_variables["Radiobutton"]["Oxidation State"].get()]
        #
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
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
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
        rb_laicpms = SimpleElements(
            parent=self.parent, row_id=25, column_id=14, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Synthetic LA-ICP-MS", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=2,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis)
        #
        self.gui_elements["Static"]["Radiobutton"].extend([rb_geophysics, rb_geochemistry, rb_laicpms])
        #
        ## Button
        btn_export = SimpleElements(
            parent=self.parent, row_id=28, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Export Data", command=lambda var_dataset=self.data_mineral, var_name=var_name:
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
        btn_simulation = SimpleElements(
            parent=self.parent, row_id=17, column_id=16, n_rows=2, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Run Simulation", command=lambda var_name=var_name: self.run_simulation_petrology(var_name))
        #
        self.gui_elements["Static"]["Button"].extend([btn_simulation])
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
        #
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
            data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone_alternative(
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
                rock="Granite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Granodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Granodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Tonalite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Tonalite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Gabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Gabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Ca")
        elif var_name == "Norite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Norite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl=[0.1, 0.5])
        elif var_name == "Diorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Diorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Na")
        elif var_name == "Monzodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Monzodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Na")
        elif var_name == "Monzogabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Monzogabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(), enrichment_pl="Ca")
        elif var_name == "Monzonite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Monzonite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Syenite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Syenite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Granitoid (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Granitoid", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Quarzolite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Quarzolite", number=self.gui_variables["Entry"]["Number Datapoints"].get())
        elif var_name == "Foid-bearing Syenite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Syenite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
        elif var_name == "Foid-bearing Monzonite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Monzonite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
        elif var_name == "Foid-bearing Monzodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Monzodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Na")
        elif var_name == "Foid-bearing Monzogabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid-bearing Monzogabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Ca")
        elif var_name == "Foid Monzosyenite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid Monzosyenite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
        elif var_name == "Foid Monzodiorite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid Monzodiorite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Na")
        elif var_name == "Foid Monzogabbro (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foid Monzogabbro", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False, enrichment_pl="Ca")
        elif var_name == "Foidolite (Streckeisen)":
            data = Plutonic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[
                    self.gui_variables["Entry"]["Porosity Min"].get()/100,
                    self.gui_variables["Entry"]["Porosity Max"].get()/100]).create_plutonic_rock_streckeisen(
                rock="Foidolite", number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                upper_streckeisen=False)
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
        #
        self.data_rock = {}
        self.rock_data = {}
        categories = ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi", "phi_true", "fluid",
                      "mineralogy", "chemistry"]
        for category in categories:
            if category in ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]:
                self.data_rock[category] = data[category]
                self.rock_data[category] = data[category]
            elif category in ["mineralogy", "chemistry"]:
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
        #
        self.list_elements_rock = list(self.data_rock["chemistry"].keys())
        self.list_minerals_rock = list(self.data_rock["mineralogy"].keys())
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
        self.gui_elements["Static"]["Radiobutton"].extend([rb_geophysics, rb_geochemistry, rb_geochemistry2])
        #
        ## Button
        btn_export = SimpleElements(
            parent=self.parent, row_id=29, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Export Data", command=lambda var_dataset=self.data_rock, var_name=var_name:
            self.export_rock_data(var_dataset, var_name))
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
        self.last_rb_analysis_rock.set(42)
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
        for property, dataset in self.rock_data.items():
            entries = [property]
            #
            if property not in ["mineralogy", "chemistry"]:
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
            dataset = self.rock_data["mineralogy"][mineral]
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
            dataset = self.rock_data["chemistry"][element]
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
        #
    #
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
        #
    #
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
        #
        fractions = {
            "Rock Forming": np.around(fraction_rock_forming/100, n_digits),
            "Ore": np.around(fraction_ore/100, n_digits),
            "Clay": np.around(fraction_clay/100, n_digits)}
        #
        mineral_list = list(selected_minerals.keys())
        mineral_list.sort()
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
            #
            val_min = selected_minerals[mineral]["Min"]
            val_max = selected_minerals[mineral]["Max"]
            mean = (val_min + val_max)/2
            sigma = (mean - val_min)/3
            #
            if mineral in self.rock_forming_minerals:
                fraction_factor = self.minerals_helper["Rock Forming"][mineral]["Fraction"]/100
                #
                val_min = self.minerals_helper["Rock Forming"][mineral]["Min"]
                val_max = self.minerals_helper["Rock Forming"][mineral]["Max"]
                mean = (val_min + val_max)/2
                sigma = (mean - val_min)/3
                #
            elif mineral in self.ore_minerals:
                fraction_factor = self.minerals_helper["Ore"][mineral]["Fraction"]/100
                #
            elif mineral in self.clay_minerals:
                fraction_factor = self.minerals_helper["Clay"][mineral]["Fraction"]/100
            #
            amounts_raw = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)
            #
            amounts_ppm = np.array([int(fraction_factor*value) for index, value in enumerate(amounts_raw)])
            amounts_ppm[amounts_ppm < 0] = 0
            amounts_ppm[amounts_ppm == 0] = 1
            #
            mineral_amounts[data_mineral["mineral"]] = list(amounts_ppm)
            mineral_data[data_mineral["mineral"]] = data_mineral
            mineral_names[mineral] = data_mineral["mineral"]
            #
        data_amounts_minerals = self.calculate_mineral_fractions(
            var_minerals=list(mineral_data.keys()), var_data=self.minerals_helper, var_n=n_datapoints)
        #
        for mineral, mineral_fractions in data_amounts_minerals.items():
            mineral_amounts[mineral_names[mineral]] = mineral_fractions
        #
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
            velocities_minerals = {}
            velocity_solid = {"vP": 0, "vS": 0}
            amounts_minerals = {"mass": {}, "volume": {}}
            sum_amounts_w = 0
            #
            for mineral, dataset in mineral_data.items():
                if mineral not in self.data_rock["mineralogy"]:
                    self.data_rock["mineralogy"][mineral] = []
                if mineral not in amounts_minerals["mass"]:
                    amounts_minerals["mass"][mineral] = []
                    amounts_minerals["volume"][mineral] = []
                #
                mineral_amount = round(mineral_amounts[mineral][n]*10**(-2), 6)
                phi_minerals[mineral] = mineral_amount
                sum_amounts_w += mineral_amount*dataset["rho"][n]
                amounts_minerals["volume"][mineral].append(mineral_amount)
                #
                velocities_minerals[mineral] = {}
                velocities_minerals[mineral]["vP"] = dataset["vP"][n]
                velocities_minerals[mineral]["vS"] = dataset["vS"][n]
                velocity_solid["vP"] += mineral_amount*dataset["vP"][n]
                velocity_solid["vS"] += mineral_amount*dataset["vS"][n]
                #
                rho_solid += mineral_amount*dataset["rho"][n]
                bulk_modulus += mineral_amount*dataset["K"][n]
                shear_modulus += mineral_amount*dataset["G"][n]
                gamma_ray += mineral_amount*dataset["GR"][n]
                photoelectricity += mineral_amount*dataset["PE"][n]
                phi_list.append(mineral_amount)
                #
                for element, value in dataset["chemistry"].items():
                    if element not in elements_list:
                        elements_list.append(element)
                        w_elements[element] = 0.0
                    #
                    if element not in self.data_rock["chemistry"]:
                        self.data_rock["chemistry"][element] = []
            #
            for mineral, dataset in mineral_data.items():
                amount_w = round((mineral_amounts[mineral][n]*10**(-2)*dataset["rho"][n])/sum_amounts_w, 6)
                amounts_minerals["mass"][mineral].append(amount_w)
            #
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
            #
            for mineral, dataset in mineral_amounts.items():
                mineral_amount = round(
                    (mineral_amounts[mineral][n]*mineral_data[mineral]["rho"][n])/(100*rho_solid), n_digits)
                w_minerals[mineral] = mineral_amount
                self.data_rock["mineralogy"][mineral].append(mineral_amount)
            #
            old_index = elements_list.index("O")
            elements_list += [elements_list.pop(old_index)]
            w_elements_total = 0.0
            for element in elements_list:
                if element != "O":
                    for mineral, w_mineral in w_minerals.items():
                        if element in mineral_data[mineral]["chemistry"]:
                            value = round(w_mineral*mineral_data[mineral]["chemistry"][element][n], n_digits)
                            w_elements[element] += value
                            w_elements_total += value
                            #
                            w_elements[element] = round(w_elements[element], n_digits)
                elif element == "O":
                    w_elements[element] += round(1 - w_elements_total, n_digits)
                    #
                    w_elements[element] = round(w_elements[element], n_digits)
            #
            for element, value in w_elements.items():
                self.data_rock["chemistry"][element].append(value)
            #
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
        #
        ## Labels
        lbl_analysis = SimpleElements(
            parent=self.parent, row_id=22, column_id=0, n_rows=6, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Analysis Mode", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Rockbuilder Static"]["Label"].append(lbl_analysis)
        #
        ## Button
        btn_export = SimpleElements(
            parent=self.parent, row_id=29, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Option"],
            fg=self.colors_gebpy["Navigation"]).create_button(
            text="Export Data", command=lambda var_dataset=self.data_rock, var_name=var_name:
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
        #
    #
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
                                    #
                                    x_min = round(x_min - 0.1*delta_x, n_digits)
                                    x_max = round(x_max + 0.1*delta_x, n_digits)
                                    #
                                    if key != "nu":
                                        if x_min < 0:
                                            x_min = 0
                                    #
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
                                    ax_rp_histo[i][j].grid(True)
                                    ax_rp_histo[i][j].set_axisbelow(True)
                                #
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
                                    #
                                    x_min = min(dataset_x)
                                    x_max = max(dataset_x)
                                    delta_x = round(x_max - x_min, 4)
                                    y_min = min(dataset_y)
                                    y_max = max(dataset_y)
                                    delta_y = round(y_max - y_min, 4)
                                    #
                                    if delta_x < 1:
                                        n_digits_x = 3
                                        x_factor = 0.9
                                    elif 1 <= delta_x < 5:
                                        n_digits_x = 2
                                        x_factor = 0.5
                                    elif delta_x >= 5:
                                        n_digits_x = 0
                                        x_factor = 0.1
                                    #
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
                                    #
                                    x_min = round(x_min - x_factor*delta_x, n_digits_x)
                                    x_max = round(x_max + x_factor*delta_x, n_digits_x)
                                    y_min = round(y_min - y_factor*delta_y, n_digits_y)
                                    y_max = round(y_max + y_factor*delta_y, n_digits_y)
                                    #
                                    if x_key in ["rho"]:
                                        x_min = round(round(min(dataset_x) - 150, -2), n_digits_x)
                                        x_max = round(round(max(dataset_x) + 150, -2), n_digits_x)
                                    if key in ["vP", "vS", "rho"]:
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
                                    #
                                    if key != "nu":
                                        if x_min < 0:
                                            x_min = 0
                                    if key != "nu":
                                        if y_min < 0:
                                            y_min = 0
                                    #
                                    step_x = (x_max - x_min)/4
                                    step_x = int(round(step_x, n_digits_x))
                                    x_ticks = np.arange(x_min, x_max + step_x, step_x)
                                    step_y = (y_max - y_min)/4
                                    if var_dtype == float:
                                        step_y = round(step_y, n_digits_y)
                                    else:
                                        step_y = int(round(step_y, n_digits_y))
                                    y_ticks = np.arange(y_min, y_max + step_y, step_y)
                                    #
                                    ax_rp_scatter[i][j].set_xlim(left=x_min, right=x_max)
                                    ax_rp_scatter[i][j].set_ylim(bottom=y_min, top=y_max)
                                    ax_rp_scatter[i][j].set_xticks(x_ticks)
                                    ax_rp_scatter[i][j].set_yticks(y_ticks)
                                    ax_rp_scatter[i][j].xaxis.set_tick_params(labelsize=8)
                                    ax_rp_scatter[i][j].yaxis.set_tick_params(labelsize=8)
                                    ax_rp_scatter[i][j].set_xlabel("Density - kg/m$^3$", fontsize=8)
                                    ax_rp_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=8)
                                    #
                                    ax_rp_scatter[i][j].grid(which="major", axis="both", linestyle="-")
                                    ax_rp_scatter[i][j].grid(which="minor", axis="both", linestyle=":", alpha=0.5)
                                    ax_rp_scatter[i][j].minorticks_on()
                                    ax_rp_scatter[i][j].xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(3))
                                    ax_rp_scatter[i][j].yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(3))
                                    ax_rp_scatter[i][j].set_axisbelow(True)
                            #
                            if "Petrology" not in self.gui_elements["Temporary"]["Canvas"]:
                                canvas_petrology = FigureCanvasTkAgg(fig_petrology, master=self.parent)
                                #
                            else:
                                canvas_petrology = self.gui_elements["Temporary"]["Canvas"]["Petrology"]
                            #
                            canvas_petrology.draw()
                            #
                            self.gui_elements["Temporary"]["Axis"]["Rock Physics Scatter"] = ax_rp_scatter
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                            #
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
                                #
                                x_min = round(x_min - 0.1*delta_x, n_digits)
                                x_max = round(x_max + 0.1*delta_x, n_digits)
                                #
                                if x_min < 0:
                                    x_min = 0
                                if key != "Urn":
                                    if x_max > 100:
                                        x_max = 100
                                #
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
                                ax_rc_histo[i][j].grid(True)
                                ax_rc_histo[i][j].set_axisbelow(True)
                            #
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
                                    #
                                    x_min = min(dataset_x)
                                    x_max = max(dataset_x)
                                    delta_x = round(x_max - x_min, 4)
                                    y_min = min(dataset_y)
                                    y_max = max(dataset_y)
                                    delta_y = round(y_max - y_min, 4)
                                    #
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
                                    #
                                    x_min = round(x_min - 0.1*delta_x, n_digits_x)
                                    x_max = round(x_max + 0.1*delta_x, n_digits_x)
                                    y_min = round(y_min - 0.1*delta_y, n_digits_y)
                                    y_max = round(y_max + 0.1*delta_y, n_digits_y)
                                    #
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
                                    ax_rc_scatter[i][j].grid(True)
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
                                #
                                y, x, _ = ax_rce_histo[i][j].hist(
                                    x=np.array(self.data_rock["chemistry"][key])*factor,
                                    color=self.colors_gebpy["Option"], edgecolor="black", bins=12)
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
                                #
                                x_min = round(x_min - 0.1*delta_x, n_digits)
                                x_max = round(x_max + 0.1*delta_x, n_digits)
                                #
                                if x_min < 0:
                                    x_min = 0
                                if key != "U":
                                    if x_max > 100:
                                        x_max = 100
                                #
                                ax_rce_histo[i][j].set_xlim(left=x_min, right=x_max)
                                ax_rce_histo[i][j].set_ylim(bottom=y_min, top=y_max)
                                ax_rce_histo[i][j].set_xticks(np.around(
                                    np.linspace(x_min, x_max, 4, dtype=float, endpoint=True), n_digits))
                                ax_rce_histo[i][j].set_yticks(np.around(
                                    np.linspace(y_min, y_max, 4, dtype=float, endpoint=True), 1))
                                ax_rce_histo[i][j].xaxis.set_tick_params(labelsize=8)
                                ax_rce_histo[i][j].yaxis.set_tick_params(labelsize=8)
                                #ax_rce_histo[i][j].set_xlabel(str(key) + " (wt.%)", fontsize=8)
                                ax_rce_histo[i][j].set_xlabel(labels[i][j], fontsize=8)
                                ax_rce_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=8)
                                ax_rce_histo[i][j].grid(True)
                                ax_rce_histo[i][j].set_axisbelow(True)
                            #
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
                            ax_rc_scatter = fig_petrology.subplots(nrows=3, ncols=3)
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
                                    #
                                    dataset_y = np.array(self.data_rock["chemistry"][key])*factor
                                    #
                                    ax_rc_scatter[i][j].scatter(
                                        dataset_x, dataset_y, color=self.colors_gebpy["Option"], edgecolor="black",
                                        alpha=0.5)
                                    #
                                    x_min = min(dataset_x)
                                    x_max = max(dataset_x)
                                    delta_x = round(x_max - x_min, 4)
                                    y_min = min(dataset_y)
                                    y_max = max(dataset_y)
                                    delta_y = round(y_max - y_min, 4)
                                    #
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
                                    #
                                    x_min = round(x_min - 0.1*delta_x, n_digits_x)
                                    x_max = round(x_max + 0.1*delta_x, n_digits_x)
                                    y_min = round(y_min - 0.1*delta_y, n_digits_y)
                                    y_max = round(y_max + 0.1*delta_y, n_digits_y)
                                    #
                                    if x_min < 0:
                                        x_min = 0
                                    if y_min < 0:
                                        y_min = 0
                                    if key != "U":
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
                                    ax_rc_scatter[i][j].set_xlabel(str(ref_element)+" (wt.%)", fontsize=8)
                                    ax_rc_scatter[i][j].set_ylabel(y_label, labelpad=0.5, fontsize=8)
                                    ax_rc_scatter[i][j].grid(True)
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
                            self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter Element"] = ax_rc_scatter
                            self.gui_elements["Temporary"]["Canvas"]["Petrology"] = canvas_petrology
                            #
                        else:
                            ## Cleaning
                            for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter Element"]:
                                for gui_axis in gui_axes:
                                    gui_axis.axis("on")
                                    gui_axis.set_visible(True)
                            #
                            categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                          "Rock Chemistry Scatter", "Rock Chemistry Histogram Element"]
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
                    for gui_axes in self.gui_elements["Temporary"]["Axis"]["Rock Chemistry Scatter Element"]:
                        for gui_axis in gui_axes:
                            gui_axis.axis("on")
                            gui_axis.set_visible(True)
                    #
                    categories = ["Rock Physics Histogram", "Rock Physics Scatter", "Rock Chemistry Histogram",
                                  "Rock Chemistry Scatter", "Rock Chemistry Histogram Element"]
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
        #
    #
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
            #
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
            #
            for index, (key, value) in enumerate(self.fractions_sum.items()):
                if self.rock_fractions[key]["Size"] > 0:
                    if index < 2:
                        self.rock_fractions[key]["Effective"] = round(
                            (self.rock_fractions[key]["Min"] + self.rock_fractions[key]["Max"])/2, 4)
                        temp_min += self.rock_fractions[key]["Min"]
                        temp_max += self.rock_fractions[key]["Max"]
                        #
                    else:
                        var_delta = 100 - w_effective
                        #
                        value_min = 100 - temp_max
                        value_max = 100 - temp_min
                        value_eff = (value_min + value_max)/2
                        #
                        self.rock_fractions[key]["Min"] = round(value_min, 4)
                        self.rock_fractions[key]["Max"] = round(value_max, 4)
                        self.rock_fractions[key]["Effective"] = round(value_eff, 4)
                        #
                        # if float(self.rock_fractions[key]["Min"]) < var_delta < float(self.rock_fractions[key]["Max"]):
                        #     delta_min = var_delta - self.rock_fractions[key]["Min"]
                        #     delta_max = self.rock_fractions[key]["Max"] - var_delta
                        #     #
                        #     if delta_max < delta_min:
                        #         self.rock_fractions[key]["Effective"] = round(100 - w_effective, 4)
                        #         self.rock_fractions[key]["Min"] = round(
                        #             2*self.rock_fractions[key]["Effective"] - self.rock_fractions[key]["Max"], 4)
                        #     elif delta_min < delta_max:
                        #         self.rock_fractions[key]["Effective"] = round(100 - w_effective, 4)
                        #         self.rock_fractions[key]["Max"] = round(
                        #             2*self.rock_fractions[key]["Effective"] - self.rock_fractions[key]["Min"], 4)
                        #     else:
                        #         self.rock_fractions[key]["Effective"] = round(var_delta, 4)
                        #         self.rock_fractions[key]["Min"] = round(var_delta - delta_min, 4)
                        #         self.rock_fractions[key]["Max"] = round(var_delta + delta_max, 4)
                        #     #
                        # elif var_delta < float(self.rock_fractions[key]["Min"]):
                        #     self.rock_fractions[key]["Effective"] = round(100 - w_effective, 4)
                        #     self.rock_fractions[key]["Min"] = round(
                        #         2*self.rock_fractions[key]["Effective"] - self.rock_fractions[key]["Max"], 4)
                        #     #
                        # elif var_delta == float(self.rock_fractions[key]["Min"]) \
                        #         or var_delta == float(self.rock_fractions[key]["Max"]):
                        #     self.rock_fractions[key]["Effective"] = round(var_delta, 4)
                        #     self.rock_fractions[key]["Min"] = round(var_delta, 4)
                        #     self.rock_fractions[key]["Max"] = round(var_delta, 4)
                    #
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
            #
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
            #
            self.minerals_helper = {}
            for index_01, (key_01, value_01) in enumerate(self.fractions_sum.items()):
                if value_01 > 0:
                    self.minerals_helper[key_01] = {}
                    for index_02, (key_02, value_02) in enumerate(self.selected_minerals.items()):
                        if key_01 == value_02["Group"]:
                            main_mineral = self.find_maximum_in_dict(var_dict=sorted_minerals[key_01])
                            self.minerals_helper[key_01][key_02] = {}
                            #
                            if value_01 == 1:
                                value_min = self.rock_fractions[key_01]["Min"]*value_02["Min"]/100
                                value_max = self.rock_fractions[key_01]["Max"]*value_02["Max"]/100
                                value_eff = (value_min + value_max)/2
                                #
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
                                #
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
            # amounts_helper = {}
            # for group, value in self.fractions_sum.items():
            #     index = 1
            #     amount_now = 0
            #     for mineral, dataset in self.selected_minerals.items():
            #         if group == dataset["Group"]:
            #             if group == "Rock Forming":
            #                 if index == n_rock_forming:
            #                     var_amount = amount_fraction_rf - amount_now
            #                 elif index == 1:
            #                     var_amount = amount_fraction_rf*(rd.uniform(
            #                         self.selected_minerals[mineral]["Min"], self.selected_minerals[mineral]["Max"])/100)
            #                 else:
            #                     var_amount = amount_fraction_rf*(rd.uniform(
            #                         self.selected_minerals[mineral]["Min"], (1 - amount_now)*100)/100)
            #                 #
            #                 amounts_helper[mineral] = round(var_amount, 4)
            #                 amount_now += round(var_amount, 4)
            #                 index += 1
            #             elif group == "Ore":
            #                 if index == n_ore:
            #                     var_amount = amount_fraction_o - amount_now
            #                 else:
            #                     var_amount = amount_fraction_o*(rd.uniform(
            #                         self.selected_minerals[mineral]["Min"], self.selected_minerals[mineral]["Max"])/100)
            #                 #
            #                 amounts_helper[mineral] = round(var_amount, 4)
            #                 amount_now += round(var_amount, 4)
            #                 index += 1
            #             elif group == "Clay":
            #                 if index == n_clay:
            #                     var_amount = amount_fraction_c - amount_now
            #                 else:
            #                     var_amount = amount_fraction_c*(rd.uniform(
            #                         self.selected_minerals[mineral]["Min"], self.selected_minerals[mineral]["Max"])/100)
            #                 #
            #                 amounts_helper[mineral] = round(var_amount, 4)
            #                 amount_now += round(var_amount, 4)
            #                 index += 1

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
        # for key, values in var_dataset.items():
        #     print(key)
        #     print(values)
        list_keys = list(var_dataset.keys())
        list_keys.remove("mineral")
        list_keys.remove("state")
        list_keys.remove("chemistry")
        list_keys = ["POISSON" if item == "nu" else item for item in list_keys]
        list_elements = list(var_dataset["chemistry"].keys())
        #
        report_file = filedialog.asksaveasfile(mode="w", initialfile="Report_Mineral", defaultextension=".csv")
        #
        ## General Data
        report_file.write("REPORT (MINERALOGY)"+"\n")
        report_file.write("Mineral"+";"+str(var_name)+"\n")
        report_file.write("\n")
        #
        ## Geophysical Data
        report_file.write("MINERAL DATA" + "\n")
        raw_line = "ID;"
        for key in list_keys:
            raw_line += str(key)
            raw_line += str(";")
        #
        raw_line += str(";")
        #
        for element in list_elements:
            raw_line += str(element)
            raw_line += str(";")
        raw_line += str("\n")
        report_file.write(raw_line)
        #
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
            #
            raw_line += str(";")
            #
            for element, values in var_dataset["chemistry"].items():
                if element in list_elements:
                    if element not in ["U"]:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")
                    else:
                        raw_line += str(round(values[index], 6))
                        raw_line += str(";")
                #
            report_file.write(raw_line+"\n")
            #
            index += 1
        #
        report_file.write("\n")
    #
    def export_rock_data(self, var_dataset, var_name):
        # for key, values in var_dataset.items():
        #     print(key)
        #     print(values)
        #
        list_keys = list(var_dataset.keys())
        list_keys.remove("mineralogy")
        list_keys.remove("chemistry")
        if "fluid" in list_keys:
            list_keys.remove("fluid")
        #
        list_keys = ["POISSON" if item == "nu" else item for item in list_keys]
        list_minerals = list(var_dataset["mineralogy"].keys())
        list_elements = list(var_dataset["chemistry"].keys())
        #
        report_file = filedialog.asksaveasfile(mode="w", initialfile="Report_Rock", defaultextension=".csv")
        #
        ## General Data
        report_file.write("REPORT (PETROLOGY)" + "\n")
        report_file.write("Rock" + ";" + str(var_name) + "\n")
        report_file.write("\n")
        #
        ## Geophysical Data
        report_file.write("ROCK DATA" + "\n")
        raw_line = "ID;"
        for key in list_keys:
            raw_line += str(key)
            raw_line += str(";")
        #
        raw_line += str(";")
        #
        for mineral in list_minerals:
            raw_line += str(mineral)
            raw_line += str(";")
        #
        raw_line += str(";")
        #
        for element in list_elements:
            raw_line += str(element)
            raw_line += str(";")
        #
        raw_line += str("\n")
        report_file.write(raw_line)
        #
        index = 0
        while index < len(var_dataset["rho"]):
            raw_line = str(index + 1) + ";"
            #
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
            #
            raw_line += str(";")
            #
            for mineral, values in var_dataset["mineralogy"].items():
                if mineral in list_minerals:
                    if mineral not in ["Urn"]:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")
                    else:
                        raw_line += str(round(values[index], 6))
                        raw_line += str(";")
            #
            raw_line += str(";")
            #
            for element, values in var_dataset["chemistry"].items():
                if element in list_elements:
                    if element not in ["U"]:
                        raw_line += str(round(values[index], 4))
                        raw_line += str(";")
                    else:
                        raw_line += str(round(values[index], 6))
                        raw_line += str(";")
            #
            report_file.write(raw_line + "\n")
            #
            index += 1
        #
        report_file.write("\n")
        #
    #
    def calculate_mineral_fractions(self, var_minerals, var_data, var_n):
        ## Input
        # print("var_minerals:", var_minerals)
        # print("var_data:", var_data)
        # print("var_n:", var_n)
        # print("selected_minerals:", self.selected_minerals)
        #
        data_amounts = {}
        #
        n_minerals = len(var_minerals)
        mineral_information = {}
        #
        for key_group, group_data in var_data.items():
            for mineral, values in group_data.items():
                if mineral not in data_amounts:
                    data_amounts[mineral] = []
                    mineral_information[mineral] = values
        #
        n = 0
        while n < var_n:
            condition = False
            while condition == False:
                temp_values = {}
                phi_total = 0
                for index, (mineral, values) in enumerate(self.selected_minerals.items()):
                    if index < n_minerals - 1:
                        #value = round(0.01*values["Fraction"]*rd.randint(values["Min"], values["Max"]), 6)
                        value = round(rd.uniform(values["Min"], values["Max"]), 6)
                    else:
                        value = round(100 - phi_total, 6)
                    #
                    # limit_lower = 0.01*values["Fraction"]*values["Min"]
                    # limit_upper = 0.01*values["Fraction"]*values["Max"]
                    limit_lower = values["Min"]
                    limit_upper = values["Max"]
                    #
                    if limit_lower <= value <= limit_upper:
                        phi_total += value
                        temp_values[mineral] = value
                #
                mineral_total = list(temp_values.values())
                if np.isclose(np.sum(mineral_total), 100.000000) == True:
                    for mineral, value in temp_values.items():
                        data_amounts[mineral].append(value)
                    n += 1
                    condition = True
        #
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
        #
    #
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
            btn_06 = SimpleElements(
                parent=self.parent, row_id=29, column_id=14, n_rows=2, n_columns=16,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                text="Run Simulation", command=lambda var_unit=var_unit: self.generate_stratigraphic_data(var_unit))
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
            btn_06 = SimpleElements(
                parent=self.parent, row_id=33, column_id=14, n_rows=2, n_columns=16,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                text="Run Simulation", command=lambda var_unit=var_unit: self.generate_stratigraphic_data(var_unit))
            #
            self.gui_elements["Temporary"]["Button"].extend([btn_06])
            #
        elif var_unit == "Buntsandstein":
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
            btn_06 = SimpleElements(
                parent=self.parent, row_id=29, column_id=14, n_rows=2, n_columns=16,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                text="Run Simulation", command=lambda var_unit=var_unit: self.generate_stratigraphic_data(var_unit))
            #
            self.gui_elements["Temporary"]["Button"].extend([btn_06])
        #
        ## INITIALIZATION
        self.generate_stratigraphic_data(var_unit=var_unit)
        self.stratigraphy_change_lithology(var_opt=self.gui_variables["Option Menu"]["Lithological Focus"].get())
        self.show_well_log_diagram()
    #
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
                    elif var_rock == "Dolomite":
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