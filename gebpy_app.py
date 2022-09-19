#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gebpy_app.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		19.09.2022

#-----------------------------------------------

## MODULES
import tkinter as tk
import collections
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from modules.gui_elements import SimpleElements
from modules.oxides import Oxides

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
        gui_priority = ["Static", "Temporary"]
        gui_subwindows = ["Trace Elements"]
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
        self.colors_gebpy = {"Background": "#EFEFEF", "Navigation": "#252422", "Accent": "#EB5E28", "Option": "#CCC5B9",
                             "White": "#FFFFFF", "Black": "#000000"}
        #
        ### Variables
        self.gui_variables = {}
        for gui_element in gui_elements:
            self.gui_variables[gui_element] = {}
        #
        ### General Settings
        self.parent = parent
        self.parent.title("GebPy")
        self.parent.geometry("1800x1000+0+0")
        self.parent.resizable(False, False)
        self.parent["bg"] = self.colors_gebpy["Background"]
        #
        ## Geometry and Layout
        window_width = 1800
        window_heigth = 1000
        row_min = 10
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
        gebpy_logo = tk.PhotoImage(file="documents/readme_images/GebPy_Logo_new.png")
        gebpy_logo = gebpy_logo.subsample(5, 5)
        img = tk.Label(self.parent, image=gebpy_logo, bg=self.colors_gebpy["Navigation"])
        img.image = gebpy_logo
        img.grid(row=0, column=0, rowspan=8, columnspan=32, sticky="nesw")
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
            "Cyclosilicates", "Inosilicates", "Phyllosilicates", "Carbonates", "Phosphates", "Phospides"]
        mineral_groups.sort()
        for mineral_group in mineral_groups:
            if mineral_group == "Oxides":
                sub_oxides = tk.Menu(sub_mineral_groups, tearoff=0)
                self.oxide_minerals = [
                    "Quartz", "Magnetite", "Hematite", "Al-Spinel", "Ilmenite", "Cassiterite", "Chromite", "Corundum",
                    "Rutile", "Pyrolusite", "Magnesiochromite", "Zincochromite", "Cr-Spinel", "Cuprospinel",
                    "Jacobsite", "Magnesioferrite", "Trevorite", "Franklinite", "Ulvospinel", "Fe-Spinel", "Uraninite",
                    "Litharge", "Massicot", "Minium", "Plattnerite", "Scrutinyite", "Zincite", "Columbite", "Tantalite",
                    "Coltan", "Crocoite", "Wulfenite", "Goethite"]
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
                    "Ankerite", "Azurite", "Malachite", "Ikaite"]
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
                    "Vaesite", "Cattierite"]
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
                    "Muscovite", "Glauconite", "Notronite", "Saponite", "Talc", "Chrysotile", "Antigorite",
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
                    "Siliciclastic": ["Sandstone", "Shale", "Mudstone"],
                    "Carbonate": ["Limestone", "Dolomite Rock"]}
                sedimentary_rocks = collections.OrderedDict(sorted(sedimentary_rocks.items()))
                i = 1
                n = len(sedimentary_rocks)
                for rock_group, rock_list in sedimentary_rocks.items():
                    rock_list.sort()
                    for rock in rock_list:
                        sub_sedimentary.add_command(
                            label=rock)
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
                            label=rock)
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
            label="Rock Builder")
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
            parent=self.parent, row_id=94, column_id=16, n_rows=4, n_columns=15,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Quit PySILLS", command=self.parent.destroy)
        #
        self.gui_elements["Static"]["Button"].append(btn_quit)
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
        self.gui_variables["Option Menu"]["Trace Elements"] = tk.StringVar()
        self.gui_variables["Option Menu"]["Trace Elements"].set("Select Trace Element")
        self.gui_variables["Checkbox"]["Trace Elements"] = {}
        self.gui_variables["Entry"]["Minimum"] = {}
        self.gui_variables["Entry"]["Maximum"] = {}
        self.gui_variables["Entry"]["Mean"] = {}
        self.gui_variables["Entry"]["Error"] = {}
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
            parent=self.parent, row_id=8, column_id=0, n_rows=4, n_columns=32, bg=self.colors_gebpy["Accent"],
            fg=self.colors_gebpy["Navigation"]).create_label(
            text="Mineralogy", font_option="sans 14 bold", relief=tk.FLAT)
        lbl_name = SimpleElements(
            parent=self.parent, row_id=12, column_id=0, n_rows=2, n_columns=32, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text=name, font_option="sans 12 bold", relief=tk.FLAT)
        lbl_samples = SimpleElements(
            parent=self.parent, row_id=14, column_id=0, n_rows=2, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Number of Datapoints", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_traces = SimpleElements(
            parent=self.parent, row_id=16, column_id=0, n_rows=6, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Trace Elements", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].extend(
            [lbl_mineralogy, lbl_name, lbl_samples, lbl_traces])
        #
        ## Entries
        entr_samples = SimpleElements(
            parent=self.parent, row_id=14, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Number Samples"])
        #
        self.gui_elements["Static"]["Entry"].extend([entr_samples])
        #
        ## Radiobuttons
        rb_trace_without = SimpleElements(
            parent=self.parent, row_id=16, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Without Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=0,
            color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
        rb_trace_with = SimpleElements(
            parent=self.parent, row_id=18, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="With Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
        #
        self.gui_elements["Static"]["Radiobutton"].extend([rb_trace_without, rb_trace_with])
        #
        ## Buttons
        btn_generate_data = SimpleElements(
            parent=self.parent, row_id=22, column_id=16, n_rows=4, n_columns=15,
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
                parent=self.parent, row_id=0, column_id=start_column, n_rows=4, n_columns=45,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Mineral Physics", font_option="sans 12 bold", relief=tk.GROOVE)
            lbl_results = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Results", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_min = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 9, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Minimum", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_max = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 18, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Maximum", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_mean = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 27, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Mean", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_error = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 36, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Error", font_option="sans 10 bold", relief=tk.GROOVE)
            #
            self.gui_elements["Temporary"]["Label"].extend(
                [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
            #
            categories = [
                "M\n (kg/mol)", "V\n (A$^3$/mol)", "rho\n (kg/m3)", "vP\n (m/s)", "vS\n (m/s)", "vP/vS\n (1)", "K\n (GPa)",
                "G\n (GPa)", "E\n (GPa)", "nu\n (1)", "GR\n (API)", "PE\n (barns/e$^-$)"]
            categories_short = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]
            for index, category in enumerate(categories):
                lbl_category = SimpleElements(
                    parent=self.parent, row_id=2*(2*index + 4), column_id=start_column, n_rows=4, n_columns=9,
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
                    parent=self.parent, row_id=2*(2*index + 4), column_id=start_column + 9, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Minimum"][categories_short[index]], var_entr_set=var_entr_min)
                entr_max = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 18, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Maximum"][categories_short[index]], var_entr_set=var_entr_max)
                entr_mean = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 27, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mean"][categories_short[index]], var_entr_set=var_entr_mean)
                entr_error = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 36, n_rows=4, n_columns=9,
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
                    row=0, column=90, rowspan=100, columnspan=90, sticky="nesw")
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
                parent=self.parent, row_id=0, column_id=start_column, n_rows=4, n_columns=45,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Mineral Chemistry", font_option="sans 12 bold", relief=tk.GROOVE)
            lbl_results = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Results", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_min = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 9, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Minimum", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_max = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 18, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Maximum", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_mean = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 27, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Mean", font_option="sans 10 bold", relief=tk.GROOVE)
            lbl_error = SimpleElements(
                parent=self.parent, row_id=4, column_id=start_column + 36, n_rows=4, n_columns=9,
                bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_label(
                text="Error", font_option="sans 10 bold", relief=tk.GROOVE)
            #
            self.gui_elements["Temporary"]["Label"].extend(
                [lbl_title, lbl_results, lbl_min, lbl_max, lbl_mean, lbl_error])
            #
            self.list_elements.sort()
            for index, element in enumerate(self.list_elements):
                lbl_element = SimpleElements(
                    parent=self.parent, row_id=2*(2*index + 4), column_id=start_column, n_rows=4, n_columns=9,
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
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 9, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Minimum"][element], var_entr_set=var_entr_min)
                entr_max = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 18, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Maximum"][element], var_entr_set=var_entr_max)
                entr_mean = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 27, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mean"][element], var_entr_set=var_entr_mean)
                entr_error = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 36, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Error"][element], var_entr_set=var_entr_error)
                #
                self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                #
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
            #
            if self.btn_traces == None:
                self.btn_traces = SimpleElements(
                    parent=self.parent, row_id=20, column_id=16, n_rows=2, n_columns=15,
                    bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
                    text="Trace Elements",
                    command=lambda var_traces=data_mineral["trace elements"]: self.select_trace_elements(var_traces))
            #
            self.gui_elements["Temporary"]["Button"].append(self.btn_traces)
    #
    def select_trace_elements(self, var_traces):
        self.window_trace_elements = tk.Toplevel(self.parent)
        self.window_trace_elements.title("Trace Elements")
        self.window_trace_elements.geometry("300x300")
        self.window_trace_elements.resizable(False, False)
        self.window_trace_elements["bg"] = self.colors_gebpy["Background"]
        #
        ## Geometry and Layout
        window_width = 300
        window_heigth = 300
        row_min = 10
        n_rows = int(window_heigth / row_min)
        column_min = 10
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
        ## Labels
        lbl_title = SimpleElements(
            parent=self.window_trace_elements, row_id=0, column_id=0, n_rows=2, n_columns=16,
            bg=self.colors_gebpy["Accent"], fg=self.colors_gebpy["Navigation"]).create_label(
            text="Trace Elements", font_option="sans 12 bold", relief=tk.GROOVE)
        #
        if lbl_title not in self.gui_elements_sub["Trace Elements"]["Static"]["Label"]:
            self.gui_elements_sub["Trace Elements"]["Static"]["Label"].append(lbl_title)
        #
        var_traces.sort()
        for index, trace_element in enumerate(var_traces):
            if index < 10:
                ## Labels
                lbl_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2*(index + 1), column_id=0, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Accent"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=trace_element, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if trace_element not in self.gui_variables["Checkbox"]["Trace Elements"]:
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element].set(0)
                cb_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1), column_id=4, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Trace Elements"][trace_element])
                #
                if lbl_trace not in self.gui_elements_sub["Trace Elements"]["Static"]["Label"]:
                    self.gui_elements_sub["Trace Elements"]["Static"]["Label"].append(lbl_trace)
                    self.gui_elements_sub["Trace Elements"]["Static"]["Checkbox"].append(cb_trace)
                #
            else:
                ## Labels
                lbl_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1) - 20, column_id=8, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Accent"], fg=self.colors_gebpy["Navigation"]).create_label(
                    text=trace_element, font_option="sans 10 bold", relief=tk.FLAT)
                #
                ## Checkboxes
                if trace_element not in self.gui_variables["Checkbox"]["Trace Elements"]:
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element] = tk.IntVar()
                    self.gui_variables["Checkbox"]["Trace Elements"][trace_element].set(0)
                cb_trace = SimpleElements(
                    parent=self.window_trace_elements, row_id=2 * (index + 1) - 20, column_id=12, n_rows=2, n_columns=4,
                    bg=self.colors_gebpy["Background"], fg=self.colors_gebpy["Navigation"]).create_checkbox(
                    text="", var_cb=self.gui_variables["Checkbox"]["Trace Elements"][trace_element])
                #
                if lbl_trace not in self.gui_elements_sub["Trace Elements"]["Static"]["Label"]:
                    self.gui_elements_sub["Trace Elements"]["Static"]["Label"].append(lbl_trace)
                    self.gui_elements_sub["Trace Elements"]["Static"]["Checkbox"].append(cb_trace)
        #
        ## Buttons
        btn_select_all = SimpleElements(
            parent=self.window_trace_elements, row_id=24, column_id=0, n_rows=4, n_columns=8,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Select All", command=lambda var_cb=self.gui_variables["Checkbox"]["Trace Elements"]:
            self.select_all_checkboxes(var_cb))
        btn_unselect_all = SimpleElements(
            parent=self.window_trace_elements, row_id=24, column_id=8, n_rows=4, n_columns=8,
            bg=self.colors_gebpy["Option"], fg=self.colors_gebpy["Navigation"]).create_button(
            text="Unselect All", command=lambda var_cb=self.gui_variables["Checkbox"]["Trace Elements"]:
            self.unselect_all_checkboxes(var_cb))
        #
        if btn_select_all not in self.gui_elements_sub["Trace Elements"]["Static"]["Button"]:
            self.gui_elements_sub["Trace Elements"]["Static"]["Button"].extend([btn_select_all, btn_unselect_all])
    #
    def run_simulation_mineralogy(self, var_name):
        self.traces_list = []
        for key, var_cb in self.gui_variables["Checkbox"]["Trace Elements"].items():
            if var_cb.get() == 1:
                self.traces_list.append(key)

        if var_name in self.oxide_minerals:
            self.data_mineral = Oxides(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
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
            parent=self.parent, row_id=26, column_id=0, n_rows=4, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Analysis Mode", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].append(lbl_analysis)
        #
        rb_geophysics = SimpleElements(
            parent=self.parent, row_id=26, column_id=16, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Physics", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=0,
            color_bg=self.colors_gebpy["Background"], command=self.change_rb_analysis)
        rb_geochemistry = SimpleElements(
            parent=self.parent, row_id=28, column_id=16, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Chemistry", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis)
        #
        self.gui_elements["Static"]["Radiobutton"].extend([rb_geophysics, rb_geochemistry])
        #
        for key, gui_element in self.gui_elements["Temporary"].items():
            if key not in ["Canvas"]:
                for gui_item in gui_element:
                    gui_item.grid_remove()
            else:
                for key_2, gui_item in gui_element.items():
                    gui_item.get_tk_widget().grid_remove()
            gui_element.clear()
    #
    ######################################
    ## G e n e r a l  F u n c t i o n s ##
    ######################################
    #
    def select_all_checkboxes(self, var_cb):
        for key, cb_var in var_cb.items():
            cb_var.set(1)
    #
    def unselect_all_checkboxes(self, var_cb):
        for key, cb_var in var_cb.items():
            cb_var.set(0)
#
if __name__ == "__main__":
    root = tk.Tk()
    #
    GebPyGUI(parent=root)
    #
    root.mainloop()