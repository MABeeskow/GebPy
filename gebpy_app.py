#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gebpy_app.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		23.11.2022

#-----------------------------------------------

## MODULES
import tkinter as tk
import collections
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from modules.gui_elements import SimpleElements
from modules.oxides import Oxides
from modules.carbonates import Carbonates, CarbonateRocks
from modules.sulfides import Sulfides
from modules.sulfates import Sulfates
from modules.halides import Halides
from modules.phospides import Phospides
from modules.phosphates import Phosphates
from modules.silicates import Phyllosilicates, Tectosilicates, Inosilicates, Nesosilicates, Sorosilicates, \
    Cyclosilicates
from modules.organics import Organics
from modules.fluids import Water
from modules.siliciclastics import Sandstone

## GUI
class GebPyGUI(tk.Frame):
    #
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        #
        ### Container
        self.gui_elements = {}
        self.gui_elements_sub = {}
        gui_elements = ["Frame", "Label", "Button", "Radiobutton", "Checkbox", "Entry", "Option Menu", "Canvas",
                        "Figure", "Axis"]
        gui_priority = ["Static", "Temporary", "Rockbuilder Static", "Rockbuilder Temporary"]
        gui_subwindows = ["Trace Elements", "Mineralogy"]
        for subwindow in gui_subwindows:
            self.gui_elements_sub[subwindow] = {}
            for priority in gui_priority:
                self.gui_elements_sub[subwindow][priority] = {}
                for gui_element in gui_elements:
                    if gui_element not in ["Canvas", "Figure", "Axis"]:
                        self.gui_elements_sub[subwindow][priority][gui_element] = []
                    else:
                        self.gui_elements_sub[subwindow][priority][gui_element] = {}
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
        self.last_rb_analysis_rock = tk.IntVar()
        self.last_rb_analysis_rock.set(42)
        #
        ### General Settings
        self.parent = parent
        self.parent.title("GebPy")
        self.parent.geometry("1800x975+0+0")
        self.parent.resizable(False, False)
        self.parent["bg"] = self.colors_gebpy["Background"]
        #
        ## Geometry and Layout
        window_width = 1800
        window_heigth = 975
        row_min = 15
        n_rows = int(window_heigth/row_min)
        column_min = 10
        n_columns = int(window_width/column_min)
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
            label="New")
        project_menu.add_command(
            label="Open")
        project_menu.add_command(
            label="Save")
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
                    "Coltan", "Crocoite", "Wulfenite", "Goethite", "Wolframite", "Huebnerite", "Ferberite", "Boehmite"]
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
                    "Mg-Fe-Pyroxene", "Ca-Pyroxene", "Donpeacorite", "Orthopyroxene"]
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
            "Spinel Group", "Hematite Group", "Rutile Group", "Periclase Group", "Wulfenite Group"]
        special_groups.sort()
        for special_group in special_groups:
            sub_special_groups.add_command(
                label=special_group)
        #
        sub_ore_minerals = tk.Menu(mineralogy_menu, tearoff=0)
        ore_minerals = [
            "Fe Ores", "Pb Ores", "Zn Ores"]
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
            "Sedimentary Rocks", "Igneous Rocks", "Metamorphic Rocks", "Evaporite Rocks"]
        rock_groups.sort()
        for rock_group in rock_groups:
            if rock_group == "Sedimentary Rocks":
                sub_sedimentary = tk.Menu(petrology_menu, tearoff=0)
                sedimentary_rocks = {
                    "Siliciclastic": ["Sandstone", "Shale", "Mudstone", "Conglomerate"],
                    "Carbonate": ["Limestone", "Dolomite Rock"]}
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
                igneous_rocks = {"Plutonic": ["Granite", "Diorite", "Gabbro"], "Volcanic": ["Rhyolite", "Basalt"]} #
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
            label="Rock Builder", command=self.rock_builder)
        #
        menubar.add_cascade(
            label="Petrology",
            menu=petrology_menu)
        #
        ## Stratigraphy
        stratigraphy_menu = tk.Menu(menubar, tearoff=0)
        stratigraphy_menu.add_command(
            label="Real Sequences")
        stratigraphy_menu.add_command(
            label="Create Sequences")

        menubar.add_cascade(
            label="Stratigraphy",
            menu=stratigraphy_menu)
        #
        ## Database
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
            parent=self.parent, row_id=62, column_id=16, n_rows=3, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Quit GebPy", command=self.parent.quit)
        btn_restart = SimpleElements(
            parent=self.parent, row_id=62, column_id=1, n_rows=3, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Restart GebPy", command=self.restart_gebpy)
        #
        self.gui_elements["Static"]["Button"].extend([btn_quit, btn_restart])
    #
    #########################
    ## M i n e r a l o g y ##
    #########################
    #
    def select_mineral(self, name):
        ## Initialization
        self.gui_variables["Entry"]["Number Samples"] = tk.IntVar()
        self.gui_variables["Entry"]["Number Samples"].set(100)
        self.gui_variables["Radiobutton"]["Analysis Mode"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Analysis Mode"].set(0)
        self.gui_variables["Radiobutton"]["Trace Elements"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Trace Elements"].set(0)
        self.gui_variables["Radiobutton"]["Oxidation State"] = tk.IntVar()
        self.gui_variables["Radiobutton"]["Oxidation State"].set(0)
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
        #
        self.gui_elements["Static"]["Label"].extend(
            [lbl_mineralogy, lbl_name, lbl_samples, lbl_traces])
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
    def change_rb_analysis(self):
        #
        start_row = 0
        start_column = 35
        #
        if self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 0:   # Mineral Physics
            ## Cleaning
            for key, gui_items in self.gui_elements["Temporary"].items():
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
            #
            self.gui_elements["Temporary"]["Label"].extend(
                [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
            #
            categories = [
                "M\n (kg/mol)", "V\n (\u00C5\u00B3/mol)", "rho\n (kg/m\u00B3)", "vP\n (m/s)", "vS\n (m/s)",
                "vP/vS\n (1)", "K\n (GPa)", "G\n (GPa)", "E\n (GPa)", "nu\n (1)", "GR\n (API)", "PE\n (barns/e\u207B)"]
            categories_short = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]
            for index, category in enumerate(categories):
                lbl_category = SimpleElements(
                    parent=self.parent, row_id=(2*index + 4), column_id=start_column, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=category, font_option="sans 10 bold", relief=tk.GROOVE)
                #
                self.gui_elements["Temporary"]["Label"].append(lbl_category)
                #
                ## Entries
                #
                var_entr_min = round(min(self.data_mineral[categories_short[index]]), 6)
                var_entr_max = round(max(self.data_mineral[categories_short[index]]), 6)
                var_entr_mean = round(np.mean(self.data_mineral[categories_short[index]]), 6)
                var_entr_error = round(np.std(self.data_mineral[categories_short[index]], ddof=1), 6)
                #
                entr_min = SimpleElements(
                    parent=self.parent, row_id=(2*index + 4), column_id=start_column + 9, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Minimum"][categories_short[index]], var_entr_set=var_entr_min)
                entr_max = SimpleElements(
                    parent=self.parent, row_id=(2*index + 4), column_id=start_column + 18, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Maximum"][categories_short[index]], var_entr_set=var_entr_max)
                entr_mean = SimpleElements(
                    parent=self.parent, row_id=(2*index + 4), column_id=start_column + 27, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mean"][categories_short[index]], var_entr_set=var_entr_mean)
                entr_error = SimpleElements(
                    parent=self.parent, row_id=(2*index + 4), column_id=start_column + 36, n_rows=2, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Error"][categories_short[index]], var_entr_set=var_entr_error)
                #
                self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                #
            ## Diagram
            if "Mineral Physics Scatter" not in self.gui_elements["Temporary"]["Canvas"]:
                fig_scatter, ax_scatter = plt.subplots(
                    ncols=3, nrows=3, figsize=(9, 9), facecolor=self.colors_gebpy["Background"])
                #
                categories = [["M", "V", "rho"], ["vP", "vS", "vP/vS"], ["GR", "PE", "nu"]]
                labels = [["M (kg/mol)", "V (A$^3$/mol", "rho (kg/m$^3$"], ["vP (m/s)", "vS (m/s)", "vP/vS (1)"],
                          ["GR (API)", "PE (barns/e$^-$)", "nu (1)"]]
                for i, subcategories in enumerate(categories):
                    for j, key in enumerate(subcategories):
                        ax_scatter[i][j].scatter(
                            self.data_mineral["rho"], self.data_mineral[key], color=self.colors_gebpy["Accent"],
                            edgecolor="black", alpha=0.5)
                        #
                        ax_scatter[i][j].set_xlabel("Density - kg/m$^3$", fontsize=9)
                        ax_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=9)
                        ax_scatter[i][j].grid(True)
                        ax_scatter[i][j].set_axisbelow(True)
                #
                fig_scatter.tight_layout()
                #
                canvas_scatter = FigureCanvasTkAgg(fig_scatter, master=self.parent)
                canvas_scatter.get_tk_widget().grid(
                    row=0, column=90, rowspan=65, columnspan=90, sticky="nesw")
                #
                self.gui_elements["Temporary"]["Canvas"]["Mineral Physics Scatter"] = canvas_scatter
                self.gui_elements["Temporary"]["Figure"]["Mineral Physics Scatter"] = fig_scatter
                self.gui_elements["Temporary"]["Axis"]["Mineral Physics Scatter"] = ax_scatter
            else:
                self.gui_elements["Temporary"]["Canvas"]["Mineral Physics Scatter"].get_tk_widget().grid()
            #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 1:   # Mineral Chemistry
            ## Cleaning
            for key, gui_items in self.gui_elements["Temporary"].items():
                if len(gui_items) > 0:
                    if key not in ["Canvas", "Figure", "Axis"]:
                        for gui_item in gui_items:
                            gui_item.grid_remove()
                        gui_items.clear()
                    else:
                        if key == "Canvas":
                            for key, gui_element in self.gui_elements["Temporary"]["Canvas"].items():
                                gui_element.get_tk_widget().grid_remove()
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
                #
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
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 2:   # Synthetic LA-ICP-MS
            ## Cleaning
            for key, gui_items in self.gui_elements["Temporary"].items():
                if len(gui_items) > 0:
                    if key not in ["Canvas", "Figure", "Axis"]:
                        for gui_item in gui_items:
                            gui_item.grid_remove()
                        gui_items.clear()
                    else:
                        if key == "Canvas":
                            for key, gui_element in self.gui_elements["Temporary"]["Canvas"].items():
                                gui_element.get_tk_widget().grid_remove()
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
                #
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
                fig_laicpms, ax_laicpms = plt.subplots(facecolor=self.colors_gebpy["Background"])
                #
                element_laicpms = list(intensity_data.keys())
                for element in element_laicpms:
                    # ax_laicpms.scatter(time_data, intensity_data[element], edgecolor="black", alpha=0.5)
                    ax_laicpms.plot(time_data, intensity_data[element], label=element)
                    #
                ax_laicpms.set_xlabel("Time (s)", fontsize=9)
                ax_laicpms.set_ylabel("Intensity (cps)", labelpad=0.5, fontsize=9)
                ax_laicpms.set_xlim(0, 60)
                ax_laicpms.set_xticks(np.arange(0, 60 + 5, 5))
                ax_laicpms.set_ylim(1, 10**6)
                ax_laicpms.set_yscale("log")
                #ax_laicpms.grid(True)
                plt.grid(which="major", axis="both", linestyle="-")
                plt.minorticks_on()
                plt.grid(which="minor", axis="both", linestyle="-", alpha=0.25)
                ax_laicpms.set_axisbelow(True)
                #
                fig_laicpms.tight_layout()
                plt.legend()
                #
                canvas_laicpms = FigureCanvasTkAgg(fig_laicpms, master=self.parent)
                canvas_laicpms.get_tk_widget().grid(
                    row=0, column=90, rowspan=65, columnspan=90, sticky="nesw")
                #
                self.gui_elements["Temporary"]["Canvas"]["LA ICP MS"] = canvas_laicpms
                self.gui_elements["Temporary"]["Figure"]["LA ICP MS"] = fig_laicpms
                self.gui_elements["Temporary"]["Axis"]["LA ICP MS"] = ax_laicpms
            else:
                self.gui_elements["Temporary"]["Canvas"]["LA ICP MS"].get_tk_widget().grid()
            #
    #
    def simulate_laicpms_experiment(self):
        total_ppm = 10**6
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
        intensity_data = {}
        for element in self.list_elements:
            index_start = 0
            index_sig = 0
            index_end = 0
            index_end2 = 0
            if element != "O":
                intensity_data[element] = []
                #
                mean_bg = np.random.randint(1, 100)
                error_bg = np.random.uniform(0.1, 0.5)*mean_bg
                mean_sig = np.mean(self.data_mineral["chemistry"][element])*10**6
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

        if var_name in self.oxide_minerals:         # Oxides
            self.data_mineral = Oxides(
                mineral=var_name, data_type=True, traces_list=self.trace_elements).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.carbonate_minerals:   # Carbonates
            self.data_mineral = Carbonates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.sulfide_minerals:   # Sulfides
            self.data_mineral = Sulfides(
                mineral=var_name, data_type=True, traces_list=self.trace_elements).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.sulfate_minerals:   # Sulfates
            self.data_mineral = Sulfates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.halide_minerals:   # Halides
            self.data_mineral = Halides(
                mineral=var_name, dict=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.phospide_minerals:   # Phospides
            self.data_mineral = Phospides(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.phosphate_minerals:   # Phosphates
            self.data_mineral = Phosphates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.tectosilicate_minerals:   # Tectosilicates
            self.data_mineral = Tectosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.nesosilicate_minerals:   # Nesosilicates
            self.data_mineral = Nesosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.sorosilicate_minerals:   # Sorosilicates
            self.data_mineral = Sorosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.inosilicate_minerals:   # Inosilicates
            self.data_mineral = Inosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.phyllosilicate_minerals:   # Phyllosilicates
            self.data_mineral = Phyllosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.cyclosilicate_minerals:   # Cyclosilicates
            self.data_mineral = Cyclosilicates(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        elif var_name in self.miscellaneous_minerals:   # Miscellaneous
            self.data_mineral = Organics(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        #
        for key, dataset in self.data_mineral.items():
            if key == "chemistry":
                self.list_elements = list(dataset.keys())
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
        self.gui_variables["Entry"]["Number Datapoints"] = tk.IntVar()
        self.gui_variables["Entry"]["Number Datapoints"].set(100)
        self.gui_variables["Entry"]["Porosity Min"] = tk.IntVar()
        self.gui_variables["Entry"]["Porosity Min"].set(0)
        self.gui_variables["Entry"]["Porosity Max"] = tk.IntVar()
        self.gui_variables["Entry"]["Porosity Max"].set(10)
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
        self.gui_elements_sub["Trace Elements"]["Static"]["Button"].extend([btn_simulation])
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
        self.selected_minerals = {}
        n_digits = 8
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
        ## Siliciclastic Rocks
        if var_name == "Sandstone":
            data = Sandstone(fluid="water", actualThickness=0).create_sandstone(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Conglomerate":
            data = Sandstone(fluid="water", actualThickness=0).create_conglomerate(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Mudstone":
            data = Sandstone(fluid="water", actualThickness=0).create_mudstone(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100], dominance="Mnt")
        elif var_name == "Shale":
            data = Sandstone(fluid="water", actualThickness=0).create_shale_alt(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        #
        ## Carbonate Rocks
        elif var_name == "Limestone":
            data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone_alternative(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        elif var_name == "Dolomite Rock":
            data = CarbonateRocks(fluid="water", actualThickness=0).create_dolomite_rock_alternative(
                number=self.gui_variables["Entry"]["Number Datapoints"].get(),
                porosity=[self.gui_variables["Entry"]["Porosity Min"].get()/100,
                          self.gui_variables["Entry"]["Porosity Max"].get()/100])
        #
        self.data_rock = {}
        categories = ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi", "fluid",
                      "mineralogy", "chemistry"]
        for category in categories:
            if category in ["rho", "rho_s", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi"]:
                self.data_rock[category] = data[category]
            elif category in ["mineralogy", "chemistry"]:
                self.data_rock[category] = data[category]
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
            text="Export Data")
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
            "Sphalerite", "Uraninite", "Wolframite", "Pollucite", "Pyrite", "Marmatite"]
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
                             "Pollucite", "Pyrite", "Marmatite"]
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
        self.selected_minerals = {}
        n_digits = 8
        #
        ## Preparation Simulation
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
                selected_minerals[mineral] = {
                    "Min": self.gui_variables["Entry"]["Mineralogy"]["Minimum"][mineral].get(),
                    "Max": self.gui_variables["Entry"]["Mineralogy"]["Maximum"][mineral].get()}
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
        fractions = {"Rock Forming": np.around(fraction_rock_forming/100, n_digits),
                     "Ore": np.around(fraction_ore/100, n_digits), "Clay": np.around(fraction_clay/100, n_digits)}
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
            amounts_raw = np.random.normal(loc=mean, scale=sigma, size=n_datapoints)
            #
            if mineral in self.rock_forming_minerals:
                fraction_factor = fractions["Rock Forming"]
            elif mineral in self.ore_minerals:
                fraction_factor = fractions["Ore"]
            elif mineral in self.clay_minerals:
                fraction_factor = fractions["Clay"]
            #
            amounts_ppm = np.array([int(fraction_factor[index]*value) for index, value in enumerate(amounts_raw)])
            amounts_ppm[amounts_ppm < 0] = 0
            #
            mineral_amounts[data_mineral["mineral"]] = list(amounts_ppm)
            mineral_data[data_mineral["mineral"]] = data_mineral
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
            #
            for mineral, dataset in mineral_data.items():
                if mineral not in self.data_rock["mineralogy"]:
                    self.data_rock["mineralogy"][mineral] = []
                #
                mineral_amount = mineral_amounts[mineral][n]*10**(-2)
                #
                rho_solid += mineral_amount*dataset["rho"][n]
                bulk_modulus += mineral_amount*dataset["K"][n]
                shear_modulus += mineral_amount*dataset["G"][n]
                gamma_ray += mineral_amount*dataset["GR"][n]
                photoelectricity += mineral_amount*dataset["PE"][n]
                #
                for element, value in dataset["chemistry"].items():
                    if element not in elements_list:
                        elements_list.append(element)
                        w_elements[element] = 0.0
                    #
                    if element not in self.data_rock["chemistry"]:
                        self.data_rock["chemistry"][element] = []
            #
            rho_solid = round(rho_solid, 3)
            bulk_modulus = round(bulk_modulus, 3)
            shear_modulus = round(shear_modulus, 3)
            gamma_ray = round(gamma_ray, 3)
            photoelectricity = round(photoelectricity, 3)
            #
            for mineral, dataset in mineral_amounts.items():
                mineral_amount = round(
                    (mineral_amounts[mineral][n]*mineral_data[mineral]["rho"][n])/(100*rho_solid), n_digits)
                w_minerals[mineral] = mineral_amount
                self.data_rock["mineralogy"][mineral].append(mineral_amount)
            #
            rho = round((1 - porosity)*rho_solid + porosity*data_water[2]/1000, 3)
            youngs_modulus = round((9*bulk_modulus*shear_modulus)/(3*bulk_modulus + shear_modulus), 3)
            poisson_ratio = round((3*bulk_modulus - 2*shear_modulus)/(6*bulk_modulus + 2*shear_modulus), 4)
            vP = round(((bulk_modulus*10**9 + 4/3*shear_modulus*10**9)/(rho))**0.5, 3)
            vS = round(((shear_modulus*10**9)/(rho))**0.5, 3)
            vPvS = round(vP/vS, 3)
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
            text="Export Data")
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
    def change_rb_analysis_rocks(self):
        start_column = 35
        #
        if self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 0:   # Rock Physics
            if self.last_rb_analysis_rock.get() != 0:
                ## Cleaning
                try:
                    self.canvas_scatter.get_tk_widget().grid_forget()
                except:
                    pass
                #
                for key, gui_items in self.gui_elements["Temporary"].items():
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
                ## Labels
                lbl_title = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=2, n_columns=45,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Rock Physics", font_option="sans 12 bold", relief=tk.GROOVE)
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
                self.gui_elements["Rockbuilder Temporary"]["Label"].extend(
                    [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
                #
                categories = [
                    "rho\n (kg/m\u00B3)", "vP\n (m/s)", "vS\n (m/s)", "vP/vS\n (1)", "K\n (GPa)", "G\n (GPa)", "E\n (GPa)",
                    "nu\n (1)", "GR\n (API)", "PE\n (barns/e\u207B)", "phi\n (%)"]
                categories_short = ["rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "phi"]
                for index, category in enumerate(categories):
                    lbl_category = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=category, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Label"].append(lbl_category)
                    #
                    ## Entries
                    #
                    var_entr_min = round(min(self.data_rock[categories_short[index]]), 3)
                    var_entr_max = round(max(self.data_rock[categories_short[index]]), 3)
                    var_entr_mean = round(np.mean(self.data_rock[categories_short[index]]), 3)
                    var_entr_error = round(np.std(self.data_rock[categories_short[index]], ddof=1), 3)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 9, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][categories_short[index]], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 18, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][categories_short[index]], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 27, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][categories_short[index]], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(3*index + 4), column_id=start_column + 36, n_rows=3, n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][categories_short[index]], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                #
                ## Diagram
                if "Rock Physics Scatter" not in self.gui_elements["Temporary"]["Canvas"]:
                    self.fig_scatter, self.ax_scatter = plt.subplots(
                        ncols=3, nrows=3, figsize=(9, 9), facecolor=self.colors_gebpy["Background"])
                    #
                    categories = [["phi", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "nu"]]
                    labels = [["$\\varphi$ (%)", "GR (API)", "PE (barns/e\u207B)"], ["vP (m/s)", "vS (m/s)", "vP/vS (1)"],
                              ["K (GPa)", "G (GPa)", "nu (1)"]]
                    for i, subcategories in enumerate(categories):
                        for j, key in enumerate(subcategories):
                            self.ax_scatter[i][j].scatter(
                                self.data_rock["rho"], self.data_rock[key], color=self.colors_gebpy["Accent"],
                                edgecolor="black", alpha=0.33)
                            #
                            self.ax_scatter[i][j].set_xlabel("Density - kg/m\u00B3", fontsize=9)
                            self.ax_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=9)
                            self.ax_scatter[i][j].grid(True)
                            self.ax_scatter[i][j].set_axisbelow(True)
                    #
                    self.fig_scatter.tight_layout()
                    #
                    self.canvas_scatter = FigureCanvasTkAgg(self.fig_scatter, master=self.parent)
                    self.canvas_scatter.get_tk_widget().grid(
                        row=0, column=90, rowspan=65, columnspan=90, sticky="nesw")
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Canvas"]["Rock Physics Scatter"] = self.canvas_scatter
                    self.gui_elements["Rockbuilder Temporary"]["Figure"]["Rock Physics Scatter"] = self.fig_scatter
                    self.gui_elements["Rockbuilder Temporary"]["Axis"]["Rock Physics Scatter"] = self.ax_scatter
                else:
                    self.gui_elements["Rockbuilder Temporary"]["Canvas"]["Rock Physics Scatter"].get_tk_widget().grid()
                #
                self.last_rb_analysis_rock.set(self.gui_variables["Radiobutton"]["Analysis Mode"].get())
                #
            else:
                self.canvas_scatter.get_tk_widget().grid()
        #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 1: # Mineral Composition
            if self.last_rb_analysis_rock.get() != 1:
                ## Cleaning
                try:
                    self.canvas_scatter.get_tk_widget().grid_forget()
                except:
                    pass
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
                ## Labels
                lbl_title = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=2, n_columns=45,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Mineral Composition", font_option="sans 12 bold", relief=tk.GROOVE)
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
                self.gui_elements["Rockbuilder Temporary"]["Label"].extend(
                    [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
                #
                self.list_minerals_rock.sort()
                for index, mineral in enumerate(self.list_minerals_rock):
                    lbl_element = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=mineral, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Label"].append(lbl_element)
                    #
                    ## Entries
                    #
                    if mineral == "Urn":
                        n_digits = 6
                    else:
                        n_digits = 2
                    #
                    var_entr_min = round(min(self.data_rock["mineralogy"][mineral])*10**2, n_digits)
                    var_entr_max = round(max(self.data_rock["mineralogy"][mineral])*10**2, n_digits)
                    var_entr_mean = round(np.mean(self.data_rock["mineralogy"][mineral])*10**2, n_digits)
                    var_entr_error = round(np.std(self.data_rock["mineralogy"][mineral], ddof=1)*10**2, n_digits)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 9, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][mineral], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 18, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][mineral], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 27, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][mineral], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 36, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][mineral], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                #
                self.last_rb_analysis_rock.set(self.gui_variables["Radiobutton"]["Analysis Mode"].get())
        #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 2: # Element Composition
            if self.last_rb_analysis_rock.get() != 2:
                ## Cleaning
                try:
                    self.canvas_scatter.get_tk_widget().grid_forget()
                except:
                    pass
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
                ## Labels
                lbl_title = SimpleElements(
                    parent=self.parent, row_id=0, column_id=start_column, n_rows=2, n_columns=45,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text="Element Composition", font_option="sans 12 bold", relief=tk.GROOVE)
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
                self.gui_elements["Rockbuilder Temporary"]["Label"].extend(
                    [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
                #
                self.list_elements_rock.sort()
                for index, element in enumerate(self.list_elements_rock):
                    lbl_element = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column, n_rows=2, n_columns=9,
                        bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                        text=element, font_option="sans 10 bold", relief=tk.GROOVE)
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Label"].append(lbl_element)
                    #
                    ## Entries
                    #
                    if element == "U":
                        n_digits = 6
                    else:
                        n_digits = 2
                    #
                    var_entr_min = round(min(self.data_rock["chemistry"][element])*10**2, n_digits)
                    var_entr_max = round(max(self.data_rock["chemistry"][element])*10**2, n_digits)
                    var_entr_mean = round(np.mean(self.data_rock["chemistry"][element])*10**2, n_digits)
                    var_entr_error = round(np.std(self.data_rock["chemistry"][element], ddof=1)*10**2, n_digits)
                    #
                    ## ppm
                    # var_entr_min = int(min(self.data_rock["chemistry"][element])*10**6)
                    # var_entr_max = int(max(self.data_rock["chemistry"][element])*10**6)
                    # var_entr_mean = int(np.mean(self.data_rock["chemistry"][element])*10**6)
                    # var_entr_error = round(np.std(self.data_rock["chemistry"][element], ddof=1)*10**6, 2)
                    #
                    entr_min = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 9, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Minimum"][element], var_entr_set=var_entr_min)
                    entr_max = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 18, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Maximum"][element], var_entr_set=var_entr_max)
                    entr_mean = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 27, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Mean"][element], var_entr_set=var_entr_mean)
                    entr_error = SimpleElements(
                        parent=self.parent, row_id=(2*index + 4), column_id=start_column + 36, n_rows=2,
                        n_columns=9,
                        bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                        var_entr=self.gui_variables["Entry"]["Error"][element], var_entr_set=var_entr_error)
                    #
                    self.gui_elements["Rockbuilder Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                #
                self.last_rb_analysis_rock.set(self.gui_variables["Radiobutton"]["Analysis Mode"].get())
    #
    ######################################
    ## G e n e r a l  F u n c t i o n s ##
    ######################################
    #
    def calculate_fractions(self, number=100, mode="simplified"):
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
            amounts_helper = {}
            for group, value in self.fractions_sum.items():
                index = 1
                amount_now = 0
                for mineral, dataset in self.selected_minerals.items():
                    if group == dataset["Group"]:
                        if group == "Rock Forming":
                            if index == n_rock_forming:
                                var_amount = amount_fraction_rf - amount_now
                            elif index == 1:
                                var_amount = amount_fraction_rf*(rd.uniform(
                                    self.selected_minerals[mineral]["Min"], self.selected_minerals[mineral]["Max"])/100)
                            else:
                                var_amount = amount_fraction_rf*(rd.uniform(
                                    self.selected_minerals[mineral]["Min"], (1 - amount_now)*100)/100)
                            #
                            amounts_helper[mineral] = round(var_amount, 4)
                            amount_now += round(var_amount, 4)
                            index += 1
                        elif group == "Ore":
                            if index == n_ore:
                                var_amount = amount_fraction_o - amount_now
                            else:
                                var_amount = amount_fraction_o*(rd.uniform(
                                    self.selected_minerals[mineral]["Min"], self.selected_minerals[mineral]["Max"])/100)
                            #
                            amounts_helper[mineral] = round(var_amount, 4)
                            amount_now += round(var_amount, 4)
                            index += 1
                        elif group == "Clay":
                            if index == n_clay:
                                var_amount = amount_fraction_c - amount_now
                            else:
                                var_amount = amount_fraction_c*(rd.uniform(
                                    self.selected_minerals[mineral]["Min"], self.selected_minerals[mineral]["Max"])/100)
                            #
                            amounts_helper[mineral] = round(var_amount, 4)
                            amount_now += round(var_amount, 4)
                            index += 1
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
        root = tk.Tk()
        GebPyGUI(root)
        root.mainloop()
    #
#
if __name__ == "__main__":
    root = tk.Tk()
    GebPyGUI(root)
    root.mainloop()