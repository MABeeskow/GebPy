#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# File:         exploration.py
# Description:  Contains all necessary functions that are related to mineral exploration feature
# Author:       Maximilian Beeskow
# Last updated: 26.01.2025
# License:      GPL v3.0

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# external packages
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.patches as mpatches
from matplotlib.figure import Figure
# internal packages
try:
    from gebpy.modules.gui_elements import SimpleElements
    from gebpy.modules.siliciclastics import SiliciclasticRocks
    from gebpy.modules.carbonates import CarbonateRocks
    from gebpy.modules.igneous import Plutonic, Volcanic, UltraMafic, Pyroclastic
    from gebpy.modules.metamorphics import GranuliteFacies, GreenschistFacies, AmphiboliteFacies
    from gebpy.modules.ore import OreRocks
    from gebpy.modules.geophysics import Seismology
except:
    from modules.gui_elements import SimpleElements
    from modules.siliciclastics import SiliciclasticRocks
    from modules.carbonates import CarbonateRocks
    from modules.igneous import Plutonic, Volcanic, UltraMafic, Pyroclastic
    from modules.metamorphics import GranuliteFacies, GreenschistFacies, AmphiboliteFacies
    from modules.ore import OreRocks

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -


class ExplorationInterface:

    def __init__(self, parent):
        self.parent = parent
        self.colors = {
            "Background": "#ECEBEA", "Navigation": "#252422", "Accent": "#EB5E28", "Option": "#CCC5B9",
            "White": "#FFFFFF", "Black": "#000000", "Accent Blue": "#118AB2"}

        # Variables
        self.initialization = False

        self.var_rb_setup = tk.IntVar()
        self.var_rb_setup.set(0)

        self.var_entr_boreholes = tk.IntVar()
        self.var_entr_boreholes.set(2)
        self.var_entr_units = tk.IntVar()
        self.var_entr_units.set(4)
        self.var_entr_maximum_depth = tk.IntVar()
        self.var_entr_maximum_depth.set(100)
        self.var_entr_parts = tk.IntVar()
        self.var_entr_parts.set(10)
        self.dict_entr_top = {}
        self.dict_entr_bottom = {}
        self.var_entr_top = tk.StringVar()
        self.var_entr_top.set("0.0")
        self.var_entr_bottom = tk.StringVar()
        self.var_entr_bottom.set("25.0")

        self.var_opt_sedimentary = tk.StringVar()
        self.var_opt_sedimentary.set("Select sedimentary rock")
        self.var_opt_plutonic = tk.StringVar()
        self.var_opt_plutonic.set("Select plutonic rock")
        self.var_opt_volcanic = tk.StringVar()
        self.var_opt_volcanic.set("Select volcanic rock")
        self.var_opt_ore = tk.StringVar()
        self.var_opt_ore.set("Select ore rock")
        self.var_opt_igneous = tk.StringVar()
        self.var_opt_igneous.set("Select igneous rock")
        self.var_opt_metamorphic = tk.StringVar()
        self.var_opt_metamorphic.set("Select metamorphic rock")
        self.var_opt_ultramafic = tk.StringVar()
        self.var_opt_ultramafic.set("Select ultramafic rock")
        self.colors_gebpy = {
            "Background": "#EFEFEF", "Navigation": "#252422", "Accent": "#EB5E28", "Option": "#CCC5B9",
            "White": "#FFFFFF", "Black": "#000000", "Accent Blue": "#118AB2"}

    def create_subwindow_borehole_data(self):
        ## Window Settings
        window_width = 1200
        window_height = 800
        var_geometry = str(window_width) + "x" + str(window_height) + "+" + str(0) + "+" + str(0)

        row_min = 25
        self.n_rows = int(window_height/row_min)
        column_min = 20
        n_columns = int(window_width/column_min)

        str_title_window = "GebPy - Borehole data"
        self.subwindow_borehole_data = tk.Toplevel(self.parent)
        self.subwindow_borehole_data.title(str_title_window)
        self.subwindow_borehole_data.geometry(var_geometry)
        self.subwindow_borehole_data.resizable(False, False)
        self.subwindow_borehole_data["bg"] = self.colors["Navigation"]

        for x in range(n_columns):
            tk.Grid.columnconfigure(self.subwindow_borehole_data, x, weight=1)
        for y in range(self.n_rows):
            tk.Grid.rowconfigure(self.subwindow_borehole_data, y, weight=1)

        # Rows
        for i in range(0, self.n_rows):
            self.subwindow_borehole_data.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            self.subwindow_borehole_data.grid_columnconfigure(i, minsize=column_min)

        self.start_row = 0
        self.n_columns_setup = 11

        ## Frames
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=self.start_row, column_id=self.n_columns_setup + 1,
            n_rows=self.n_rows, n_columns=n_columns - self.n_columns_setup - 1, bg=self.colors["Background"],
            fg=self.colors["Background"]).create_frame()

        ## Labels
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=self.start_row, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Settings", relief=tk.FLAT, font_option="sans 12 bold", anchor_option=tk.W)

        ## Initialization
        self.container_borehole_lithology = {}
        self.container_lithology_data = {}
        if self.initialization == False:
            self.select_dataset_setup()
            self.update_settings(initialization=True)
            self.initialization = True

    def select_dataset_setup(self):
        if self.var_rb_setup.get() == 0:
            ## Labels
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 1, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Number of boreholes (#)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 3, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Number of units (#)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 5, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Maximum depth (m)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 7, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Samples per unit (#)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)

            ## Entries
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 2, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_boreholes)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 4, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_units)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 6, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_maximum_depth)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 8, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_parts)

            ## Buttons
            btn_01 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 10, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_button(
                text="Update settings", command=lambda initialization=False, update=True:
                self.update_settings(initialization, update))
            btn_01.configure(activebackground=self.colors["Accent"])

        elif self.var_rb_setup.get() == 1:
            print("Load dataset!")

    def update_settings(self, initialization=False, update=False, changed_borehole=False):
        """Updates the data for the borehole and unit data.
        Arguments
            initialization, boolean : specifies if the function runs the first time.
            update, boolean : specifies if it has only to update already existing variables.
            changed_borehole, boolean : specifies if only the ID of a borehole was changed.
        Outputs
            ---
        """
        ## Helper
        if initialization or update and changed_borehole == False:
            self.list_boreholes = np.arange(self.var_entr_boreholes.get()) + 1
            self.list_units = np.arange(self.var_entr_units.get()) + 1
            self.current_borehole_id = self.list_boreholes[0]
            self.current_unit_id = self.list_units[0]

            for index_borehole, borehole_id in enumerate(self.list_boreholes):
                self.dict_entr_top[borehole_id] = {}
                self.dict_entr_bottom[borehole_id] = {}
                average_thickness = round(self.var_entr_maximum_depth.get()/len(self.list_units), 2)
                current_bottom_depth = 0
                for index_unit, unit_id in enumerate(self.list_units):
                    self.dict_entr_top[borehole_id][unit_id] = tk.IntVar()
                    self.dict_entr_bottom[borehole_id][unit_id] = tk.IntVar()

                    if index_unit == 0:
                        self.dict_entr_top[borehole_id][unit_id].set(0)
                        current_bottom_depth += average_thickness
                        self.dict_entr_bottom[borehole_id][unit_id].set(round(current_bottom_depth))
                    else:
                        self.dict_entr_top[borehole_id][unit_id].set(round(current_bottom_depth))
                        current_bottom_depth += average_thickness
                        self.dict_entr_bottom[borehole_id][unit_id].set(round(current_bottom_depth))

        if initialization:
            ## Labels
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 12, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Current borehole", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            self.lbl_borehole_id = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 13, column_id=1, n_rows=1, n_columns=2,
                bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text=self.current_borehole_id, relief=tk.FLAT, font_option="sans 10 bold")
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 14, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Current unit", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            self.lbl_unit_id = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 15, column_id=1, n_rows=1, n_columns=2,
                bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text=self.current_unit_id, relief=tk.FLAT, font_option="sans 10 bold")
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 17, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Lithology", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 25, column_id=1, n_rows=1, n_columns=6,
                bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Top depth (m)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 26, column_id=1, n_rows=1, n_columns=6,
                bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Bottom depth (m)", relief=tk.FLAT, font_option="sans 10 bold", anchor_option=tk.W)

            ## Buttons
            btn_02 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 13, column_id=3, n_rows=1, n_columns=4,
                bg=self.colors["Background"], fg=self.colors["Navigation"]).create_button(
                text="Previous", command=lambda mode="previous": self.change_borehole(mode))
            btn_03 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 13, column_id=7, n_rows=1, n_columns=4,
                bg=self.colors["Background"], fg=self.colors["Navigation"]).create_button(
                text="Next", command=lambda mode="next": self.change_borehole(mode))
            btn_04 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 15, column_id=3, n_rows=1, n_columns=4,
                bg=self.colors["Background"], fg=self.colors["Navigation"]).create_button(
                text="Previous", command=lambda mode="previous": self.change_unit(mode))
            btn_05 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 15, column_id=7, n_rows=1, n_columns=4,
                bg=self.colors["Background"], fg=self.colors["Navigation"]).create_button(
                text="Next", command=lambda mode="next": self.change_unit(mode))
            btn_06 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 28, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_button(text="Load drilling data", command=self.load_drilling_data)
            btn_07 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 29, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_button(text="Show data", command=self.show_data)
            btn_08 = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 30, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_button(text="Export data", command=self.export_borehole_data)
            btn_02.configure(activebackground=self.colors["Accent"])
            btn_03.configure(activebackground=self.colors["Accent"])
            btn_04.configure(activebackground=self.colors["Accent"])
            btn_05.configure(activebackground=self.colors["Accent"])
            btn_06.configure(activebackground=self.colors["Accent"])
            btn_07.configure(activebackground=self.colors["Accent"])
            btn_08.configure(activebackground=self.colors["Accent"])

            ## Entries
            self.entr_top = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 25, column_id=7, n_rows=1, n_columns=4,
                bg=self.colors["Background"], fg=self.colors["Navigation"]).create_entry(
                var_entr=self.dict_entr_top[self.current_borehole_id][self.current_unit_id])
            self.entr_bottom = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 26, column_id=7, n_rows=1, n_columns=4,
                bg=self.colors["Background"], fg=self.colors["Navigation"]).create_entry(
                var_entr=self.dict_entr_bottom[self.current_borehole_id][self.current_unit_id])
            self.entr_top.bind(
                "<Return>", lambda event, mode="top": self.change_depth(mode, event))
            self.entr_bottom.bind(
                "<Return>", lambda event, mode="bottom": self.change_depth(mode, event))

            if self.current_unit_id == 1:
                self.entr_top.configure(state="disabled")
            else:
                self.entr_top.configure(state="normal")

            if self.current_unit_id == self.list_units[-1]:
                self.entr_bottom.configure(state="disabled")
            else:
                self.entr_bottom.configure(state="normal")

            ## Option Menus
            list_opt_sedimentary = [
                "Sandstone", "Shale", "Siltstone", "Mudstone", "Greywacke (Huckenholz)", "Conglomerate", "Limestone",
                "Dolostone", "Marl"]
            list_opt_plutonic = [
                "Granite", "Granodiorite", "Tonalite", "Gabbro", "Norite", "Diorite", "Monzodiorite", "Monzogabbro",
                "Monzonite", "Syenite", "Granitoid", "Quarzolite", "Foid-bearing Syenite", "Foid-bearing Monzonite",
                "Foid-bearing Monzodiorite", "Foid-bearing Monzogabbro", "Foid Monzosyenite", "Foid Monzodiorite",
                "Foid Monzogabbro", "Foidolite"]
            list_opt_volcanic = [
                "Rhyolite", "Dacite", "Trachyte", "Latite", "Andesite", "Basalt",
                "Foid-bearing Trachyte", "Foid-bearing Latite", "Foid-bearing Andesite", "Foid-bearing Basalt",
                "Phonolite", "Tephrite", "Foidite"]
            list_opt_ultramafic = [
                "Orthopyroxenite", "Clinopyroxenite", "Dunite", "Harzburgite", "Wehrlite", "Websterite", "Lherzolite",
                "Olivine-Websterite", "Olivine-Orthopyroxenite", "Olivine-Clinopyroxenite", "Peridotite", "Pyroxenite"]
            list_opt_metamorphic = [
                "Felsic Granulite", "Mafic Granulite", "Basaltic Greenschist", "Ultramafic Greenschist",
                "Pelitic Greenschist", "Greenstone", "Ortho-Amphibolite"]
            list_opt_ore = [
                "Itabirite", "Banded Iron Formation", "Compact Hematite", "Friable Hematite", "Goethite Hematite",
                "Al-rich Itabirite", "Compact Quartz Itabirite", "Friable Quartz Itabirite", "Goethite Itabirite"]
            opt_sedimentary = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 18, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_option_menu(
                var_opt=self.var_opt_sedimentary, var_opt_set=self.var_opt_sedimentary.get(),
                opt_list=list_opt_sedimentary, active_bg=self.colors["Accent"],
                command=lambda event, focus="sedimentary": self.select_lithology(focus, event))
            opt_plutonic = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 19, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_option_menu(
                var_opt=self.var_opt_plutonic, var_opt_set=self.var_opt_plutonic.get(), opt_list=list_opt_plutonic,
                active_bg=self.colors["Accent"], command=lambda event, focus="plutonic":
                self.select_lithology(focus, event))
            opt_volcanic = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 20, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_option_menu(
                var_opt=self.var_opt_volcanic, var_opt_set=self.var_opt_volcanic.get(), opt_list=list_opt_volcanic,
                active_bg=self.colors["Accent"], command=lambda event, focus="volcanic":
                self.select_lithology(focus, event))
            opt_ultramafic = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 21, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_option_menu(
                var_opt=self.var_opt_ultramafic, var_opt_set=self.var_opt_ultramafic.get(),
                opt_list=list_opt_ultramafic, active_bg=self.colors["Accent"],
                command=lambda event, focus="ultramafic": self.select_lithology(focus, event))
            opt_metamorphic = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 22, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_option_menu(
                var_opt=self.var_opt_metamorphic, var_opt_set=self.var_opt_metamorphic.get(),
                opt_list=list_opt_metamorphic, active_bg=self.colors["Accent"],
                command=lambda event, focus="metamorphic": self.select_lithology(focus, event))
            opt_ore = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row + 23, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_option_menu(
                var_opt=self.var_opt_ore, var_opt_set=self.var_opt_ore.get(),
                opt_list=list_opt_ore, active_bg=self.colors["Accent"],
                command=lambda event, focus="ore": self.select_lithology(focus, event))
            opt_sedimentary.configure(anchor=tk.W)
            opt_plutonic.configure(anchor=tk.W)
            opt_volcanic.configure(anchor=tk.W)
            opt_ultramafic.configure(anchor=tk.W)
            opt_metamorphic.configure(anchor=tk.W)
            opt_ore.configure(anchor=tk.W)

            ## Treeviews (Rock data)
            self.categories = ["rho (kg/m\u00B3)", "phi (%)", "vP (m/s)", "vS (m/s)", "vP/vS (1)", "K (GPa)", "G (GPa)",
                          "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
            list_categories = ["Category", "Minimum", "Maximum", "Mean", "Standard Deviation"]
            list_width = list(75*np.ones(len(list_categories)))
            list_width = [int(item) for item in list_width]
            list_width[0] = 90
            list_width[-1] = 135

            self.tv_lithology = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row, column_id=self.n_columns_setup + 1,
                n_rows=self.n_rows - 1, n_columns=2*self.n_columns_setup + 1, fg=self.colors["Black"],
                bg=self.colors["White"]).create_treeview(
                n_categories=len(list_categories), text_n=list_categories, width_n=list_width, individual=True)

            scb_v = ttk.Scrollbar(self.subwindow_borehole_data, orient="vertical")
            scb_h = ttk.Scrollbar(self.subwindow_borehole_data, orient="horizontal")
            self.tv_lithology.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
            scb_v.config(command=self.tv_lithology.yview)
            scb_h.config(command=self.tv_lithology.xview)
            scb_v.grid(
                row=0, column=3*self.n_columns_setup + 2, rowspan=self.n_rows - 1, columnspan=1, sticky="ns")
            scb_h.grid(
                row=self.n_rows - 1, column=self.n_columns_setup + 1, rowspan=1, columnspan=2*self.n_columns_setup + 1,
                sticky="ew")

            for index, category in enumerate(self.categories):
                entries = [category]
                n_digits = 3
                var_entr_min = round(0, n_digits)
                var_entr_max = round(0, n_digits)
                var_entr_mean = round(0, n_digits)
                var_entr_error = round(0, n_digits)

                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                self.tv_lithology.insert("", tk.END, values=entries)

            ## Treeviews (Borehole data)
            list_categories = ["Borehole", "Unit", "Top", "Bottom", "Lithology"]
            list_width = list(75*np.ones(len(list_categories)))
            list_width = [int(item) for item in list_width]
            list_width[0] = 75
            list_width[-1] = 150

            self.tv_borehole = SimpleElements(
                parent=self.subwindow_borehole_data, row_id=self.start_row, column_id=3*self.n_columns_setup + 3,
                n_rows=self.n_rows - 1, n_columns=2*self.n_columns_setup + 1, fg=self.colors["Black"],
                bg=self.colors["White"]).create_treeview(
                n_categories=len(list_categories), text_n=list_categories, width_n=list_width, individual=True)

            scb_v = ttk.Scrollbar(self.subwindow_borehole_data, orient="vertical")
            scb_h = ttk.Scrollbar(self.subwindow_borehole_data, orient="horizontal")
            self.tv_borehole.configure(xscrollcommand=scb_h.set, yscrollcommand=scb_v.set)
            scb_v.config(command=self.tv_borehole.yview)
            scb_h.config(command=self.tv_borehole.xview)
            scb_v.grid(
                row=0, column=5*self.n_columns_setup + 4, rowspan=self.n_rows - 1, columnspan=1, sticky="ns")
            scb_h.grid(
                row=self.n_rows - 1, column=3*self.n_columns_setup + 3, rowspan=1,
                columnspan=2*self.n_columns_setup + 1, sticky="ew")

            n = 0
            for i in self.list_boreholes:
                if i not in self.container_borehole_lithology:
                    self.container_borehole_lithology[i] = {}
                if n == 4:
                    entries = ["-", "-", "-", "-", "-"]
                    self.tv_borehole.insert("", tk.END, values=entries)
                    n = 0
                for j in self.list_units:
                    name = "undefined"

                    if j not in self.container_borehole_lithology[i]:
                        self.container_borehole_lithology[i][j] = tk.StringVar()
                        self.container_borehole_lithology[i][j].set(name)

                    var_entr_top = self.dict_entr_top[i][j].get()
                    var_entr_bottom = self.dict_entr_bottom[i][j].get()
                    var_entr_name = self.container_borehole_lithology[i][j].get()
                    entries = [i, j, var_entr_top, var_entr_bottom, var_entr_name]
                    self.tv_borehole.insert("", tk.END, values=entries)
                    n += 1

        if update:
            self.lbl_borehole_id.configure(text=self.current_borehole_id)
            self.lbl_unit_id.configure(text=self.current_unit_id)
            self.entr_top.configure(textvariable=self.dict_entr_top[self.current_borehole_id][self.current_unit_id])
            self.entr_bottom.configure(textvariable=self.dict_entr_bottom[self.current_borehole_id][
                self.current_unit_id])
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_igneous.set("Select igneous rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.update_borehole_table()

    def change_borehole(self, mode):
        """Changes the borehole.
        Arguments
            mode, str : specifies which borehole was selected.
        Outputs
            ---
        """
        if mode == "previous":
            if self.current_borehole_id == 1:
                self.current_borehole_id = self.list_boreholes[-1]
            else:
                self.current_borehole_id = self.current_borehole_id - 1
        elif mode == "next":
            if self.current_borehole_id == self.list_boreholes[-1]:
                self.current_borehole_id = 1
            else:
                self.current_borehole_id += 1

        self.lbl_borehole_id.configure(text=self.current_borehole_id)
        self.current_unit_id = self.list_units[0]

        try:
            self.fill_lithology_table(dataset=self.data, categories_long=self.categories)
        except:
            print("Please select first a lithology.")

        self.update_settings(initialization=False, update=True, changed_borehole=True)

    def change_unit(self, mode):
        """Changes the unit within a borehole.
        Arguments
            mode, str : specifies which unit was selected.
        Outputs
            ---
        """
        if mode == "previous":
            if self.current_unit_id == 1:
                self.current_unit_id = self.list_units[-1]
            else:
                self.current_unit_id -= 1

            self.entr_top.configure(textvariable=self.dict_entr_top[self.current_borehole_id][self.current_unit_id])
            self.entr_bottom.configure(textvariable=self.dict_entr_bottom[self.current_borehole_id][
                self.current_unit_id])

            if self.current_unit_id == 1:
                self.entr_top.configure(state="disabled")
            else:
                self.entr_top.configure(state="normal")

            if self.current_unit_id == self.list_units[-1]:
                self.entr_bottom.configure(state="disabled")
            else:
                self.entr_bottom.configure(state="normal")

        elif mode == "next":
            if self.current_unit_id == self.list_units[-1]:
                self.current_unit_id = 1
            else:
                self.current_unit_id += 1

            self.entr_top.configure(textvariable=self.dict_entr_top[self.current_borehole_id][self.current_unit_id])
            self.entr_bottom.configure(textvariable=self.dict_entr_bottom[self.current_borehole_id][
                self.current_unit_id])

            if self.current_unit_id == 1:
                self.entr_top.configure(state="disabled")
            else:
                self.entr_top.configure(state="normal")

            if self.current_unit_id == self.list_units[-1]:
                self.entr_bottom.configure(state="disabled")
            else:
                self.entr_bottom.configure(state="normal")

        self.lbl_unit_id.configure(text=self.current_unit_id)

        try:
            self.fill_lithology_table(dataset=self.data, categories_long=self.categories)
        except:
            print("Please select first a lithology.")
        
    def change_depth(self, mode, event):
        """Changes the depth of the top or bottom.
        Arguments
            mode, str : specifies if top or bottom was selected.
        Outputs
            ---
        """
        if mode == "top":
            if self.current_unit_id >= 2:
                new_depth = self.dict_entr_top[self.current_borehole_id][self.current_unit_id].get()
                current_top_depth_above = self.dict_entr_top[self.current_borehole_id][self.current_unit_id - 1].get()

                if new_depth > current_top_depth_above:
                    self.dict_entr_bottom[self.current_borehole_id][self.current_unit_id - 1].set(new_depth)
                else:
                    self.dict_entr_top[self.current_borehole_id][self.current_unit_id].set(current_top_depth_above + 1)
                    self.dict_entr_bottom[self.current_borehole_id][self.current_unit_id - 1].set(
                        current_top_depth_above + 1)

        elif mode == "bottom":
            if self.current_unit_id <= len(self.list_units):
                new_depth = self.dict_entr_bottom[self.current_borehole_id][self.current_unit_id].get()
                current_bottom_depth_below = self.dict_entr_bottom[self.current_borehole_id][
                    self.current_unit_id + 1].get()

                if new_depth < current_bottom_depth_below:
                    self.dict_entr_top[self.current_borehole_id][self.current_unit_id + 1].set(new_depth)
                else:
                    self.dict_entr_bottom[self.current_borehole_id][self.current_unit_id].set(
                        current_bottom_depth_below - 1)
                    self.dict_entr_top[self.current_borehole_id][self.current_unit_id + 1].set(
                        current_bottom_depth_below - 1)

        self.update_borehole_table()

    def select_lithology(self, focus, event):
        """Selects a lithology (e.g. sandstone, granite, ...).
        Arguments
            focus, str : specifies which lithology was selected.
        Outputs
            ---
        """
        if len(self.tv_lithology.get_children()) > 0:
            for item in self.tv_lithology.get_children():
                self.tv_lithology.delete(item)

        n_datapoints = self.var_entr_parts.get()
        self.categories = ["rho (kg/m\u00B3)", "phi (%)", "vP (m/s)", "vS (m/s)", "vP/vS (1)", "K (GPa)", "G (GPa)",
                      "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
        self.categories_short = ["rho", "phi", "vP", "vS", "vPvS", "K", "G", "E", "v", "GR", "PE"]
        borehole_id = self.current_borehole_id
        unit_id = self.current_unit_id

        if focus == "sedimentary":
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            # Siliciclastic rocks
            self.data = self.generate_rock_data(rockname=self.var_opt_sedimentary.get())

            self.container_borehole_lithology[borehole_id][unit_id].set(self.var_opt_sedimentary.get())

            if self.var_opt_sedimentary.get() not in self.container_lithology_data:
                self.container_lithology_data[self.var_opt_sedimentary.get()] = self.data.copy()
        elif focus == "plutonic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            self.data = self.generate_rock_data(rockname=self.var_opt_plutonic.get())

            self.container_borehole_lithology[borehole_id][unit_id].set(self.var_opt_plutonic.get())

            if self.var_opt_plutonic.get() not in self.container_lithology_data:
                self.container_lithology_data[self.var_opt_plutonic.get()] = self.data.copy()
        elif focus == "volcanic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            self.data = self.generate_rock_data(rockname=self.var_opt_volcanic.get())

            self.container_borehole_lithology[borehole_id][unit_id].set(self.var_opt_volcanic.get())

            if self.var_opt_volcanic.get() not in self.container_lithology_data:
                self.container_lithology_data[self.var_opt_volcanic.get()] = self.data.copy()
        elif focus == "ultramafic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            self.data = self.generate_rock_data(rockname=self.var_opt_ultramafic.get())

            self.container_borehole_lithology[borehole_id][unit_id].set(self.var_opt_ultramafic.get())

            if self.var_opt_ultramafic.get() not in self.container_lithology_data:
                self.container_lithology_data[self.var_opt_ultramafic.get()] = self.data.copy()
        elif focus == "metamorphic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_ore.set("Select ore rock")

            self.data = self.generate_rock_data(rockname=self.var_opt_metamorphic.get())

            self.container_borehole_lithology[borehole_id][unit_id].set(self.var_opt_metamorphic.get())

            if self.var_opt_metamorphic.get() not in self.container_lithology_data:
                self.container_lithology_data[self.var_opt_metamorphic.get()] = self.data.copy()
        elif focus == "ore":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")

            self.data = self.generate_rock_data(rockname=self.var_opt_ore.get())

            self.container_borehole_lithology[borehole_id][unit_id].set(self.var_opt_ore.get())

            if self.var_opt_ore.get() not in self.container_lithology_data:
                self.container_lithology_data[self.var_opt_ore.get()] = self.data.copy()

        try:
            self.fill_lithology_table(dataset=self.data, categories_long=self.categories)
        except:
            print("Please select first a lithology.")

        self.update_borehole_table()

    def update_borehole_table(self):
        if len(self.tv_borehole.get_children()) > 0:
            for item in self.tv_borehole.get_children():
                self.tv_borehole.delete(item)

        n = 0
        for i in self.list_boreholes:
            if i not in self.container_borehole_lithology:
                self.container_borehole_lithology[i] = {}
            if n == len(self.list_units):
                entries = ["-", "-", "-", "-", "-"]
                self.tv_borehole.insert("", tk.END, values=entries)
                n = 0
            for j in self.list_units:
                if j in self.container_borehole_lithology[i]:
                    name = self.container_borehole_lithology[i][j]
                else:
                    name = "undefined"

                if j not in self.container_borehole_lithology[i]:
                    self.container_borehole_lithology[i][j] = tk.StringVar()
                    self.container_borehole_lithology[i][j].set(name)

                var_entr_top = self.dict_entr_top[i][j].get()
                var_entr_bottom = self.dict_entr_bottom[i][j].get()
                var_entr_name = self.container_borehole_lithology[i][j].get()
                entries = [i, j, var_entr_top, var_entr_bottom, var_entr_name]
                self.tv_borehole.insert("", tk.END, values=entries)
                n += 1

    def fill_lithology_table(self, dataset, categories_long):
        borehole_id = self.current_borehole_id
        unit_id = self.current_unit_id
        rockname = self.container_borehole_lithology[borehole_id][unit_id].get()

        if rockname != "undefined":
            if rockname in self.container_lithology_data:
                dataset = self.container_lithology_data[rockname]
            else:
                dataset = self.generate_rock_data(rockname=rockname)
                self.container_lithology_data[rockname] = dataset

        if len(self.tv_lithology.get_children()) > 0:
            for item in self.tv_lithology.get_children():
                self.tv_lithology.delete(item)

        for index, category in enumerate(categories_long):
            entries = [category]

            if self.categories_short[index] == "vPvS":
                category_short = "vP/vS"
            elif self.categories_short[index] == "v":
                category_short = "nu"
            else:
                category_short = self.categories_short[index]

            n_digits = 3

            if rockname == "undefined":
                var_entr_min = round(0, n_digits)
                var_entr_max = round(0, n_digits)
                var_entr_mean = round(0, n_digits)
                var_entr_error = round(0, n_digits)
            else:
                var_entr_min = round(min(dataset[category_short]), n_digits)
                var_entr_max = round(max(dataset[category_short]), n_digits)
                var_entr_mean = round(np.mean(dataset[category_short]), n_digits)
                var_entr_error = round(np.std(dataset[category_short], ddof=1), n_digits)

            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            self.tv_lithology.insert("", tk.END, values=entries)

        entries = ["-", "-", "-", "-", "-"]
        self.tv_lithology.insert("", tk.END, values=entries)

        for mineral, data_values in dataset["mineralogy"].items():
            entries = [str(mineral) + str(" (%)")]

            n_digits = 2
            var_factor = 100

            if rockname == "undefined":
                var_entr_min = round(0, n_digits)
                var_entr_max = round(0, n_digits)
                var_entr_mean = round(0, n_digits)
                var_entr_error = round(0, n_digits)
            else:
                var_entr_min = round(var_factor*min(data_values), n_digits)
                var_entr_max = round(var_factor*max(data_values), n_digits)
                var_entr_mean = round(var_factor*np.mean(data_values), n_digits)
                var_entr_error = round(var_factor*np.std(data_values, ddof=1), n_digits)

            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            self.tv_lithology.insert("", tk.END, values=entries)

        entries = ["-", "-", "-", "-", "-"]
        self.tv_lithology.insert("", tk.END, values=entries)

        for element, data_values in dataset["chemistry"].items():
            entries = [str(element) + str(" (%)")]

            n_digits = 2
            var_factor = 100

            if rockname == "undefined":
                var_entr_min = round(0, n_digits)
                var_entr_max = round(0, n_digits)
                var_entr_mean = round(0, n_digits)
                var_entr_error = round(0, n_digits)
            else:
                var_entr_min = round(var_factor*min(data_values), n_digits)
                var_entr_max = round(var_factor*max(data_values), n_digits)
                var_entr_mean = round(var_factor*np.mean(data_values), n_digits)
                var_entr_error = round(var_factor*np.std(data_values, ddof=1), n_digits)

            entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
            self.tv_lithology.insert("", tk.END, values=entries)

        if "compounds" in dataset:
            entries = ["-", "-", "-", "-", "-"]
            self.tv_lithology.insert("", tk.END, values=entries)

            for compound, data_values in dataset["compounds"].items():
                entries = [str(compound) + str(" (%)")]

                n_digits = 2
                var_factor = 100

                if rockname == "undefined":
                    var_entr_min = round(0, n_digits)
                    var_entr_max = round(0, n_digits)
                    var_entr_mean = round(0, n_digits)
                    var_entr_error = round(0, n_digits)
                else:
                    var_entr_min = round(var_factor*min(data_values), n_digits)
                    var_entr_max = round(var_factor*max(data_values), n_digits)
                    var_entr_mean = round(var_factor*np.mean(data_values), n_digits)
                    var_entr_error = round(var_factor*np.std(data_values, ddof=1), n_digits)

                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
                self.tv_lithology.insert("", tk.END, values=entries)

    def export_borehole_data(self):
        if self.list_units.all() != None:
            report_file = filedialog.asksaveasfile(
                mode="w",
                initialfile="Report_Borehole-Data_" + str(len(self.list_boreholes)) + str(len(self.list_units)),
                defaultextension=".csv")
        else:
            report_file = filedialog.asksaveasfile(
                mode="w", initialfile="Report_Borehole-Data", defaultextension=".csv")

        ## General Data
        report_file.write("REPORT (BOREHOLE DATA)" + "\n")
        report_file.write("\n")

        ## Geophysical Data
        report_file.write("ROCK DATA" + "\n")
        raw_line = "BOREHOLE;UNIT;SAMPLE;TOP;BOTTOM;"

        for borehole_id in self.list_boreholes:
            if self.list_units.all() == None:
                list_units = self.dict_indices[borehole_id]
            else:
                list_units = self.list_units
            for unit_id in list_units:
                first_unit_id = list_units[0]
                if first_unit_id == 0:
                    correction = 1
                else:
                    correction = 0

                if unit_id in self.container_borehole_lithology[borehole_id]:
                    rockname = self.container_borehole_lithology[borehole_id][unit_id].get()

                    if rockname != "undefined":
                        dataset = self.generate_rock_data(rockname=rockname)
                        n_samples = len(dataset["rho"])
                    else:
                        n_samples = 0

                    list_keys = list(dataset.keys())
                    list_keys.remove("mineralogy")
                    list_keys.remove("chemistry")

                    if "compounds" in list_keys:
                        list_keys.remove("compounds")
                    if "fluid" in list_keys:
                        list_keys.remove("fluid")
                    if "phi_true" in list_keys:
                        list_keys.remove("phi_true")

                    list_keys = ["POISSON" if item == "nu" else item for item in list_keys]
                    list_minerals = list(dataset["mineralogy"].keys())
                    list_elements = list(dataset["chemistry"].keys())
                    list_compounds = list(dataset["compounds"].keys())

                    for key in list_keys:
                        raw_line += str(key)
                        raw_line += str(";")

                    for mineral in list_minerals:
                        raw_line += str(mineral)
                        raw_line += str(";")

                    for element in list_elements:
                        raw_line += str(element)
                        raw_line += str(";")

                    for compound in list_compounds:
                        raw_line += str(compound)
                        raw_line += str(";")

                    raw_line += str("\n")

                    report_file.write(raw_line)

                    raw_line = ""
                    for index in range(n_samples):
                        raw_line = (str(borehole_id) + ";" + str(unit_id + correction) + ";" + str(index + 1) + ";" +
                                    str(self.dict_entr_top[borehole_id][unit_id].get()) + ";" +
                                    str(self.dict_entr_bottom[borehole_id][unit_id].get()) + ";")
                        for key, values in dataset.items():
                            if key in list_keys:
                                if key in ["phi", "rho_s", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "GR", "PE"]:
                                    raw_line += str(round(values[index], 3))
                                else:
                                    raw_line += str(values)

                                raw_line += str(";")
                            elif key == "nu":
                                raw_line += str(round(values[index], 3))
                                raw_line += str(";")

                        for mineral in list_minerals:
                            if mineral not in ["Urn"]:
                                value = dataset["mineralogy"][mineral][index]
                                raw_line += str(round(value, 4))
                                raw_line += str(";")
                            else:
                                raw_line += str(round(value, 6))
                                raw_line += str(";")

                        for element, values in dataset["chemistry"].items():
                            if element in list_elements:
                                if element not in ["U"]:
                                    raw_line += str(round(values[index], 4))
                                    raw_line += str(";")
                                else:
                                    raw_line += str(round(values[index], 6))
                                    raw_line += str(";")

                        if "compounds" in dataset:
                            for compound, values in dataset["compounds"].items():
                                if compound in list_compounds:
                                    raw_line += str(round(values[index], 4))
                                    raw_line += str(";")

                        report_file.write(raw_line + "\n")

                    raw_line = "BOREHOLE;UNIT;SAMPLE;TOP;BOTTOM;"

    def generate_rock_data(self, rockname, test=False):
        if test == False:
            n_samples = self.var_entr_parts.get()
        else:
            n_samples = 1

        ## Sedimentary rocks
        # Siliciclastic rocks
        if rockname == "Sandstone":
            dataset = SiliciclasticRocks(fluid="water", actualThickness=0).create_sandstone(
                number=n_samples, porosity=[0.0, 0.3])
        elif rockname == "Conglomerate":
            dataset = SiliciclasticRocks(fluid="water", actualThickness=0).create_conglomerate(
                number=n_samples, porosity=[0.0, 0.3])
        elif rockname == "Siltstone":
            dataset = SiliciclasticRocks(fluid="water", actualThickness=0).create_siltstone(
                number=n_samples, porosity=[0.0, 0.1])
        elif rockname == "Mudstone":
            dataset = SiliciclasticRocks(fluid="water", actualThickness=0).create_mudstone_alt(
                number=n_samples, porosity=[0.0, 0.1])
        elif rockname == "Shale":
            dataset = SiliciclasticRocks(fluid="water", actualThickness=0).create_shale_alt(
                number=n_samples, porosity=[0.0, 0.1])
        elif rockname == "Greywacke (Huckenholz)":
            dataset = SiliciclasticRocks(fluid="water", actualThickness=0).create_greywacke_huckenholz(
                rock="Greywacke", number=n_samples, porosity=[0.0, 0.1])
        # Carbonate rocks
        elif rockname == "Limestone":
            dataset = CarbonateRocks(fluid="water", actualThickness=0).create_limestone(
                number=n_samples,  porosity=[0.0, 0.4])
        elif rockname == "Dolostone":
            dataset = CarbonateRocks(fluid="water", actualThickness=0).create_dolostone(
                number=n_samples, porosity=[0.0, 0.3])
        elif rockname == "Marl":
            dataset = CarbonateRocks(fluid="water", actualThickness=0).create_marl(
                number=n_samples, porosity=[0.0, 0.3])

        ## Igneous rocks (plutonic)
        elif rockname in [
            "Foid-bearing Syenite", "Foid-bearing Monzonite", "Foid-bearing Monzodiorite",
            "Foid-bearing Monzogabbro", "Foid Monzosyenite", "Foid Monzodiorite", "Foid Monzogabbro",
            "Foidolite"]:
            dataset = Plutonic(
                fluid="water", actualThickness=0, dict_output=True,
                porosity=[0.0, 0.1]).create_plutonic_rock_streckeisen(
                rock=rockname, number=n_samples, porosity=[0.0, 0.1], upper_streckeisen=False)
        elif rockname in [
                "Granite", "Granodiorite", "Tonalite", "Gabbro", "Norite", "Diorite", "Monzodiorite", "Monzogabbro",
                "Monzonite", "Syenite", "Granitoid", "Quarzolite"]:
            dataset = Plutonic(
                fluid="water", actualThickness=0, dict_output=True,
                porosity=[0.0, 0.1]).create_plutonic_rock_streckeisen(
                rock=rockname, number=n_samples, porosity=[0.0, 0.1])

        ## Igneous rocks (volcanic)
        elif rockname in [
            "Foid-bearing Trachyte", "Foid-bearing Latite", "Foid-bearing Andesite", "Foid-bearing Basalt",
            "Phonolite", "Tephrite", "Foidite"]:
            dataset = Volcanic(
                fluid="water", actualThickness=0, dict_output=True,
                porosity=[0.0, 0.1]).create_volcanic_rock_streckeisen(
                rock=rockname, number=n_samples, upper_streckeisen=False)
        elif rockname in ["Rhyolite", "Dacite", "Trachyte", "Latite", "Andesite", "Basalt"]:
            dataset = Volcanic(
                fluid="water", actualThickness=0, dict_output=True,
                porosity=[0.0, 0.1]).create_volcanic_rock_streckeisen(rock=rockname, number=n_samples)

        ## Ultramafic rocks
        elif rockname in [
            "Orthopyroxenite", "Clinopyroxenite", "Dunite", "Harzburgite", "Wehrlite", "Websterite", "Lherzolite",
            "Olivine-Websterite", "Olivine-Orthopyroxenite", "Olivine-Clinopyroxenite", "Peridotite", "Pyroxenite"]:
            dataset = UltraMafic(
                fluid="water", actualThickness=0, dict_output=True, porosity=[0.0, 0.1]).create_ultramafic_rock(
                rock=rockname, number=n_samples)

        ## Metamorphic rocks
        elif rockname == "Felsic Granulite":
            dataset = GranuliteFacies(fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_granulite(
                number=n_samples, classification="felsic")
        elif rockname == "Mafic Granulite":
            dataset = GranuliteFacies(fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_granulite(
                number=n_samples, classification="mafic")
        elif rockname == "Basaltic Greenschist":
            dataset = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenschist_basaltic_alt(
                number=n_samples)
        elif rockname == "Ultramafic Greenschist":
            dataset = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenschist_ultramafic_alt(
                number=n_samples)
        elif rockname == "Pelitic Greenschist":
            dataset = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenschist_pelitic_alt(
                number=n_samples)
        elif rockname == "Greenstone":
            dataset = GreenschistFacies(
                fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenstone(
                number=n_samples)
        # Amphibolite-Facies
        elif rockname == "Ortho-Amphibolite":
            dataset = AmphiboliteFacies(
                fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_amphibolite_ortho(
                number=n_samples)

        ## Ore rocks
        if rockname in [
            "Itabirite", "Banded Iron Formation", "Compact Hematite", "Friable Hematite", "Goethite Hematite",
            "Al-rich Itabirite", "Compact Quartz Itabirite", "Friable Quartz Itabirite", "Goethite Itabirite"]:
            dataset = OreRocks(fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_siliciclastic_itabirite(
                rock=rockname, number=n_samples)

        ## Evaporite rocks

        return dataset

    def load_drilling_data(self):
        filename = filedialog.askopenfilenames(
            parent=self.parent,
            filetypes=(("csv files", "*.csv"), ("txt files", "*.txt"), ("all files", "*.*")), initialdir=os.getcwd())

        df = pd.read_csv(filename[0])
        borehole_indices = self.get_borehole_indices(dataset=df["BOREHOLE"])
        helper_indices = {}
        helper_indices_raw = {}
        for index, borehole_id in enumerate(borehole_indices):
            list_indices = self.get_borehole_id_indices(dataset=df["BOREHOLE"], value=borehole_id)
            helper_indices[borehole_id] = list_indices
            helper_indices_raw[borehole_id] = list_indices
            if index > 0:
                diff_to_1 = list_indices[0] - 1
                helper_indices[borehole_id] = list(np.array(list_indices) - diff_to_1)
            else:
                helper_indices[borehole_id] = list(np.array(list_indices) + 1)

        self.list_units = None
        self.list_units = []
        self.dict_indices = helper_indices

        if len(self.tv_borehole.get_children()) > 0:
            for item in self.tv_borehole.get_children():
                self.tv_borehole.delete(item)

        n = 0
        columns_mandatory = ["BOREHOLE", "TOP", "BOTTOM", "LITHOLOGY"]
        self.helper_element_concentrations = {}
        for i in borehole_indices:
            if i not in self.helper_element_concentrations:
                self.helper_element_concentrations[i] = {}
            if i not in self.container_borehole_lithology:
                self.container_borehole_lithology[i] = {}
            if i not in self.dict_entr_top:
                self.dict_entr_top[i] = {}
            if i not in self.dict_entr_bottom:
                self.dict_entr_bottom[i] = {}
            for index, j in enumerate(helper_indices[i]):
                row_id = helper_indices_raw[i][index]
                j = index
                if j not in self.list_units:
                    self.list_units.append(j)
                columns_optionally = [col for col in df.columns if col not in columns_mandatory]
                if j not in self.helper_element_concentrations[i]:
                    self.helper_element_concentrations[i][j] = {}
                if j in self.container_borehole_lithology[i]:
                    name = self.container_borehole_lithology[i][j]
                else:
                    name = "undefined"

                if j not in self.container_borehole_lithology[i]:
                    self.container_borehole_lithology[i][j] = tk.StringVar()
                    self.container_borehole_lithology[i][j].set(name)
                if j not in self.dict_entr_top:
                    self.dict_entr_top[i][j] = tk.StringVar()
                if j not in self.dict_entr_bottom:
                    self.dict_entr_bottom[i][j] = tk.StringVar()

                n_columns = df.shape[1]     # Number of columns
                var_entr_top = df["TOP"].iloc[row_id]
                var_entr_bottom = df["BOTTOM"].iloc[row_id]
                var_entr_name = df["LITHOLOGY"].iloc[row_id]

                if n_columns > 4:
                    list_elements = [col[1:] for col in columns_optionally]
                    for index, element in enumerate(list_elements):
                        self.helper_element_concentrations[i][j][element] = df[columns_optionally[index]].iloc[row_id]

                rockname = self.convert_name(name=var_entr_name)
                self.container_borehole_lithology[i][j].set(rockname)
                self.dict_entr_top[i][j].set(round(var_entr_top, 1))
                self.dict_entr_bottom[i][j].set(round(var_entr_bottom, 1))
                entries = [i, j + 1, "{:.1f}".format(var_entr_top), "{:.1f}".format(var_entr_bottom), rockname]
                self.tv_borehole.insert("", tk.END, values=entries)
                n += 1

                if n == len(helper_indices[i]):
                    entries = ["-", "-", "-", "-", "-"]
                    self.tv_borehole.insert("", tk.END, values=entries)
                    n = 0

        self.list_units = np.array(self.list_units)

    def get_borehole_indices(self, dataset):
        """Creates a list containing the unique numbers of a given input list.
        Arguments
            dataset, list : contains the data.
        Outputs
            borehole_indices, list : contains a list of the unique numbers of the given input list
        """
        unique_numbers = dataset.unique()
        borehole_indices = unique_numbers.tolist()
        self.list_boreholes = borehole_indices

        return borehole_indices

    def get_borehole_id_indices(self, dataset, value):
        list_indices = dataset.index[dataset == value].tolist()

        return list_indices

    def convert_name(self, name):
        if name in ["BIF", "bif"]:
            rockname = "Banded Iron Formation"
        elif name in ["GSB", "GS", "gsb", "gs"]:
            rockname = "Greenstone"

        return rockname

    def create_final_dataset(self):
        data_borehole = {}
        self.color_lithology = {}
        for borehole_id in self.list_boreholes:
            data_borehole[borehole_id] = {
                "Lithology": [], "Top": [], "Bottom": [], "rho": [], "rho_s": [], "fluid": [], "vP": [], "vS": [],
                "K": [], "G": [], "E": [], "nu": [], "GR": [], "PE": [], "phi": [], "mineralogy": {}, "chemistry": {},
                "compounds": {}, "AI": [], "RC": []}
            list_all_minerals = []
            list_all_elements = []
            list_all_compounds = []

            if self.list_units.all() == None:
                list_units = self.dict_indices[borehole_id]
            else:
                list_units = self.list_units

            for unit_id in list_units:
                if unit_id in self.container_borehole_lithology[borehole_id]:
                    str_lithology = self.container_borehole_lithology[borehole_id][unit_id].get()

                    if str_lithology != "undefined":
                        if str_lithology not in self.color_lithology:
                            self.color_lithology[str_lithology] = {"Color": None}

                            ## Sedimentary rocks
                            # Siliciclastic rocks
                            if str_lithology in ["Sandstone", "Siltstone"]:
                                self.color_lithology[str_lithology] = {"Color": "tan"}
                            elif str_lithology in ["Shale", "Mudstone"]:
                                self.color_lithology[str_lithology] = {"Color": "olivedrab"}
                            elif str_lithology in ["Conglomerate"]:
                                self.color_lithology[str_lithology] = {"Color": "tan"}
                            elif "Greywacke" in str_lithology:
                                self.color_lithology[str_lithology] = {"Color": "tan"}
                            # Carbonate rocks
                            elif str_lithology in ["Limestone"]:
                                self.color_lithology[str_lithology] = {"Color": "skyblue"}
                            elif str_lithology in ["Limestone"]:
                                self.color_lithology[str_lithology] = {"Color": "lightcyan"}
                            elif str_lithology in ["Marl"]:
                                self.color_lithology[str_lithology] = {"Color": "moccasin"}
                            # Evaporitic rocks
                            elif str_lithology in ["Anhydrite"]:
                                self.color_lithology[str_lithology] = {"Color": "orchid"}
                            elif str_lithology in ["Rock Salt"]:
                                self.color_lithology[str_lithology] = {"Color": "lavender"}
                            elif str_lithology in ["Potash"]:
                                self.color_lithology[str_lithology] = {"Color": "yellowgreen"}
                            ## Igneous rocks
                            elif str_lithology in ["Granite", "alpha-Granite", "beta-Granite", "Gabbro", "Diorite"]:
                                self.color_lithology[str_lithology] = {"Color": "darkorange"}
                            ## Metamorphic rocks
                            elif str_lithology in ["Basaltic Greenschist", "Ultramafic Greenschist", "Pelitic Greenschist",
                                                   "Greenstone"]:
                                self.color_lithology[str_lithology] = {"Color": "darkgreen"}
                            elif str_lithology in ["Felsic Granulite", "Mafic Granulite"]:
                                self.color_lithology[str_lithology] = {"Color": "silver"}
                            elif str_lithology in ["Ortho-Amphibolite"]:
                                self.color_lithology[str_lithology] = {"Color": "teal"}
                            ## Ore rocks
                            elif str_lithology in [
                                "Itabirite", "Compact Hematite", "Friable Hematite", "Goethite Hematite",
                                "Al-rich Itabirite", "Compact Quartz Itabirite", "Friable Quartz Itabirite",
                                "Goethite Itabirite", "Banded Iron Formation"]:
                                self.color_lithology[str_lithology] = {"Color": "gray"}
                            elif str_lithology in ["Bauxite"]:
                                self.color_lithology[str_lithology] = {"Color": "peachpuff"}



                        dataset = self.generate_rock_data(rockname=str_lithology, test=True)
                        list_minerals = list(dataset["mineralogy"].keys())
                        list_elements = list(dataset["chemistry"].keys())
                        list_compounds = list(dataset["compounds"].keys())

                        for mineral in list_minerals:
                            if mineral not in list_all_minerals:
                                list_all_minerals.append(mineral)
                                data_borehole[borehole_id]["mineralogy"][mineral] = []

                        for element in list_elements:
                            if element not in list_all_elements:
                                list_all_elements.append(element)
                                data_borehole[borehole_id]["chemistry"][element] = []

                        for compound in list_compounds:
                            if compound not in list_all_compounds:
                                list_all_compounds.append(compound)
                                data_borehole[borehole_id]["compounds"][compound] = []

            for unit_id in list_units:
                if unit_id in self.container_borehole_lithology[borehole_id]:
                    str_lithology = self.container_borehole_lithology[borehole_id][unit_id].get()
                    val_top = float(self.dict_entr_top[borehole_id][unit_id].get())
                    val_bottom = float(self.dict_entr_bottom[borehole_id][unit_id].get())

                    if str_lithology != "undefined":
                        dataset = self.generate_rock_data(rockname=str_lithology)
                        n_samples = len(dataset["rho"])
                    else:
                        n_samples = 0

                    if n_samples > 0:
                        val_step = round((val_bottom - val_top)/n_samples, 2)
                        for index in range(n_samples):
                            val_sub_top = val_top + index*val_step
                            val_sub_bottom = val_top + (index + 1)*val_step
                            val_rho = round(dataset["rho"][index], 3)
                            val_rho_s = round(dataset["rho_s"][index], 3)
                            str_fluid = dataset["fluid"]
                            val_vP = round(dataset["vP"][index], 3)
                            val_vS = round(dataset["vS"][index], 3)
                            val_K = round(dataset["K"][index], 3)
                            val_G = round(dataset["G"][index], 3)
                            val_E = round(dataset["E"][index], 3)
                            val_poisson = dataset["nu"][index]
                            val_GR = round(dataset["GR"][index], 3)
                            val_PE = round(dataset["PE"][index], 3)
                            val_phi = dataset["phi"][index]
                            val_AI = round(val_rho*val_vP, 3)

                            if len(data_borehole[borehole_id]["AI"]) > 0:
                                val_last_AI = data_borehole[borehole_id]["AI"][-1]
                            else:
                                val_last_AI = 0

                            val_RC = (val_last_AI - val_AI)/(val_last_AI + val_AI)

                            list_minerals = list(dataset["mineralogy"].keys())
                            list_elements = list(dataset["chemistry"].keys())
                            list_compounds = list(dataset["compounds"].keys())

                            data_borehole[borehole_id]["Lithology"].append(str_lithology)
                            data_borehole[borehole_id]["Top"].append(val_sub_top)
                            data_borehole[borehole_id]["Bottom"].append(val_sub_bottom)
                            data_borehole[borehole_id]["rho"].append(val_rho)
                            data_borehole[borehole_id]["rho_s"].append(val_rho_s)
                            data_borehole[borehole_id]["fluid"].append(str_fluid)
                            data_borehole[borehole_id]["vP"].append(val_vP)
                            data_borehole[borehole_id]["vS"].append(val_vS)
                            data_borehole[borehole_id]["K"].append(val_K)
                            data_borehole[borehole_id]["G"].append(val_G)
                            data_borehole[borehole_id]["E"].append(val_E)
                            data_borehole[borehole_id]["nu"].append(val_poisson)
                            data_borehole[borehole_id]["GR"].append(val_GR)
                            data_borehole[borehole_id]["PE"].append(val_PE)
                            data_borehole[borehole_id]["phi"].append(val_phi)
                            data_borehole[borehole_id]["AI"].append(val_AI)

                            if len(data_borehole[borehole_id]["RC"]) > 0:
                                data_borehole[borehole_id]["RC"].append(val_RC)
                            else:
                                data_borehole[borehole_id]["RC"].append(0.0)

                            for mineral in list_all_minerals:
                                if mineral in list_minerals:
                                    val_mineral = dataset["mineralogy"][mineral][index]
                                    data_borehole[borehole_id]["mineralogy"][mineral].append(val_mineral)
                                else:
                                    data_borehole[borehole_id]["mineralogy"][mineral].append(0.0)

                            for element in list_all_elements:
                                if element in list_elements:
                                    val_element = dataset["chemistry"][element][index]
                                    data_borehole[borehole_id]["chemistry"][element].append(val_element)
                                else:
                                    data_borehole[borehole_id]["chemistry"][element].append(0.0)

                            for compound in list_all_compounds:
                                if compound in list_compounds:
                                    val_compound = dataset["compounds"][compound][index]
                                    data_borehole[borehole_id]["compounds"][compound].append(val_compound)
                                else:
                                    data_borehole[borehole_id]["compounds"][compound].append(0.0)
                    else:
                        pass

        return data_borehole

    def show_data(self):
        ## Window Settings
        window_width = 1200
        window_height = 800
        var_geometry = str(window_width) + "x" + str(window_height) + "+" + str(0) + "+" + str(0)

        row_min = 25
        self.n_rows_diagram = int(window_height/row_min)
        column_min = 20
        self.n_columns_diagram = int(window_width/column_min)

        str_title_window = "GebPy - Borehole data"
        self.subwindow_borehole_diagrams = tk.Toplevel(self.parent)
        self.subwindow_borehole_diagrams.title(str_title_window)
        self.subwindow_borehole_diagrams.geometry(var_geometry)
        self.subwindow_borehole_diagrams.resizable(False, False)
        self.subwindow_borehole_diagrams["bg"] = self.colors["Navigation"]

        for x in range(self.n_columns_diagram):
            tk.Grid.columnconfigure(self.subwindow_borehole_diagrams, x, weight=1)
        for y in range(self.n_rows_diagram):
            tk.Grid.rowconfigure(self.subwindow_borehole_diagrams, y, weight=1)

        # Rows
        for i in range(0, self.n_rows_diagram):
            self.subwindow_borehole_diagrams.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, self.n_columns_diagram):
            self.subwindow_borehole_diagrams.grid_columnconfigure(i, minsize=column_min)

        self.start_row = 0
        self.n_columns_setup = 11

        ## Frames
        SimpleElements(
            parent=self.subwindow_borehole_diagrams, row_id=self.start_row, column_id=self.n_columns_setup + 1,
            n_rows=self.n_rows, n_columns=self.n_columns_diagram - self.n_columns_setup - 1,
            bg=self.colors["Background"], fg=self.colors["Background"]).create_frame()

        ## Labels
        SimpleElements(
            parent=self.subwindow_borehole_diagrams, row_id=self.start_row, column_id=1, n_rows=1,
            n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Settings", relief=tk.FLAT, font_option="sans 12 bold", anchor_option=tk.W)

        self.data_borehole = self.create_final_dataset()
        self.show_well_log_diagram()

    def show_well_log_diagram(self, borehole_id=1):
        self.canvas = None
        max_thickness = max(self.data_borehole[borehole_id]["Bottom"])

        if max_thickness <= 100:
            step_depth = 10
        elif 100 < max_thickness <= 500:
            step_depth = 50
        elif 500 < max_thickness <= 1500:
            step_depth = 100
        elif max_thickness > 1500:
            step_depth = 200

        data_seismic = Seismology().create_seismic_trace_new(data_reflection=self.data_borehole[borehole_id]["RC"])

        self.fig, (self.ax1, self.ax2, self.ax3, self.ax4, self.ax5, self.ax6, self.ax7) = plt.subplots(
            1, 7, sharey="row", gridspec_kw={"wspace": 0.25}, figsize=(12, 24),
            facecolor=self.colors_gebpy["Background"])
        self.fig.subplots_adjust(wspace=0.25)
        # 1
        self.ax1.plot(self.data_borehole[borehole_id]["GR"], self.data_borehole[borehole_id]["Top"], color="#00549F",
                      linewidth=2)
        self.ax1.set_xlabel("GR [API]")
        self.ax1.set_ylabel("Depth [m]")

        if max(self.data_borehole[borehole_id]["GR"]) > 200:
            self.ax1.set_xlim(-0.5, max(self.data_borehole[borehole_id]["GR"]))
            self.ax1.set_xticks(np.arange(-0.5, max(self.data_borehole[borehole_id]["GR"]) + 200, 200))
            self.ax1.set_xscale("symlog")
        elif 50 < max(self.data_borehole[borehole_id]["GR"]) < 100:
            self.ax1.set_xlim(-0.5, max(self.data_borehole[borehole_id]["GR"]))
            self.ax1.set_xticks(np.arange(-0.5, max(self.data_borehole[borehole_id]["GR"]) + 20, 20))
            self.ax1.set_xscale("symlog")
        else:
            self.ax1.set_xlim(-1, max(self.data_borehole[borehole_id]["GR"]))
            self.ax1.set_xticks(np.arange(0, max(self.data_borehole[borehole_id]["GR"]) + 10, 10))

        self.ax1.set_ylim(0, max_thickness)
        self.ax1.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax1.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 2
        vP_edit = np.array(self.data_borehole[borehole_id]["vP"])/1000
        vS_edit = np.array(self.data_borehole[borehole_id]["vS"])/1000
        self.ax2.plot(vP_edit, self.data_borehole[borehole_id]["Top"], color="#00549F", linewidth=2)
        self.ax2.set_xlabel("$v_P$ [km/s]")
        self.ax2.set_xlim(0, max(vP_edit))
        self.ax2.set_xticks(np.arange(0, max(vP_edit) + 2.0, 2.0))
        self.ax2.xaxis.label.set_color("#00549F")
        self.ax2.set_ylim(0, max_thickness)
        self.ax2.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax2.grid(color="grey", linestyle="dashed")
        self.ax2_2 = self.ax2.twiny()
        self.ax2_2.plot(vS_edit, self.data_borehole[borehole_id]["Top"], color="#CC071E", linewidth=2)
        self.ax2_2.set_xlabel("$v_S$ [km/s]")
        self.ax2_2.set_xlim(0, max(vP_edit))
        self.ax2_2.set_xticks(np.arange(0, max(vP_edit) + 2.0, 2.0))
        self.ax2_2.minorticks_on()
        self.ax2_2.xaxis.label.set_color("#CC071E")
        self.ax2_2.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # # 3
        phi_edit = np.array(self.data_borehole[borehole_id]["phi"])*100
        self.ax3.plot(np.array(self.data_borehole[borehole_id]["rho"])/1000, self.data_borehole[borehole_id]["Top"],
                      color="#57AB27", linewidth=2)
        self.ax3.set_xlabel("$\\varrho$ [g/cm$^3$]")
        self.ax3.set_xlim(1.7, 3.2)
        self.ax3.set_xticks(np.around(np.linspace(1.7, 3.2, 4, endpoint=True), decimals=1))
        self.ax3.xaxis.label.set_color("#57AB27")
        self.ax3.set_ylim(0, max_thickness)
        self.ax3.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax3.grid(color="grey", linestyle="dashed")
        self.ax3_2 = self.ax3.twiny()
        self.ax3_2.plot(phi_edit, self.data_borehole[borehole_id]["Top"], color="#00549F", linewidth=2)
        self.ax3_2.set_xlabel("$\\varphi$ [1]")
        self.ax3_2.set_xlim(60, -30)
        self.ax3_2.set_xticks(np.around(np.linspace(60, -30, 6, endpoint=True), decimals=0))
        self.ax3_2.minorticks_on()
        self.ax3_2.xaxis.label.set_color("#00549F")
        self.ax3_2.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # # 4
        self.ax4.plot(self.data_borehole[borehole_id]["PE"], self.data_borehole[borehole_id]["Top"], color="#00549F",
                      linewidth=2)
        self.ax4.set_xlabel("PE [barns/electron]")

        if max(self.data_borehole[borehole_id]["PE"]) > 10:
            self.ax4.set_xlim(-0.5, max(self.data_borehole[borehole_id]["PE"]))
            self.ax4.set_xticks(np.arange(-0.5, len(str(int(max(self.data_borehole[borehole_id]["PE"])))), 3))
            self.ax4.set_xscale("symlog")
        else:
            self.ax4.set_xlim(0, max(self.data_borehole[borehole_id]["PE"]))
            self.ax4.set_xticks(np.arange(0, max(self.data_borehole[borehole_id]["PE"]) + 1, 1))

        self.ax4.set_ylim(0, max_thickness)
        self.ax4.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax4.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 5
        AI = np.asarray(self.data_borehole[borehole_id]["AI"])/10**6
        self.ax5.plot(AI, self.data_borehole[borehole_id]["Top"], color="#006165", linewidth=2)
        self.ax5.set_xlabel("AI (kNs/m$^3$)")
        self.ax5.set_xlim(0, 25)
        self.ax5.set_xticks(np.around(np.linspace(0, 25, 6, endpoint=True), decimals=0))
        self.ax5.set_ylim(0, max_thickness)
        self.ax5.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax5.grid(color="grey", linestyle="dashed")
        self.ax5.minorticks_on()
        plt.gca().invert_yaxis()
        plt.rc('axes', axisbelow=True)
        # 6
        # RC = np.asarray(self.rock_data["All Rocks"]["Physics"]["RC"])
        # self.ax6.plot(RC, self.stratigraphy_data["Top"], color="#006165", linewidth=2)
        # self.ax6.set_xlabel("RC (-)")
        # self.ax6.set_xlim(-0.5, 0.5)
        # self.ax6.set_xticks(np.around(np.linspace(-0.5, 0.5, 3, endpoint=True), decimals=1))
        # self.ax6.set_ylim(0, max_thickness)
        # self.ax6.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        # self.ax6.grid(color="grey", linestyle="dashed")
        # self.ax6.minorticks_on()
        # plt.gca().invert_yaxis()
        # plt.rc('axes', axisbelow=True)
        self.ax6.plot(data_seismic, self.data_borehole[borehole_id]["Top"], color="black", linewidth=1)
        self.ax6.fill_betweenx(self.data_borehole[borehole_id]["Top"], 0.0, data_seismic, where=(data_seismic > 0.0),
                               color="royalblue")
        self.ax6.fill_betweenx(self.data_borehole[borehole_id]["Top"], 0.0, data_seismic, where=(data_seismic < 0.0),
                               color="indianred")
        self.ax6.set_xlabel("Seismic trace")
        self.ax6.set_ylim(0, max_thickness)
        self.ax6.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax6.grid(color="grey", linestyle="dashed")
        self.ax6.minorticks_on()
        plt.gca().invert_yaxis()
        plt.rc('axes', axisbelow=True)
        # 7
        n_units = len(self.data_borehole[borehole_id]["Lithology"])
        legend_lithology = []

        for key, item in self.color_lithology.items():
            legend_lithology.append(
                mpatches.Patch(facecolor=item["Color"], edgecolor="black", hatch="", label=key))

        for index, val_top in enumerate(self.data_borehole[borehole_id]["Top"]):
            str_lithology = self.data_borehole[borehole_id]["Lithology"][index]
            val_bottom = self.data_borehole[borehole_id]["Bottom"][index]
            self.ax7.hist(
                x=np.linspace(val_top, val_bottom), bins=2, color=self.color_lithology[str_lithology]["Color"],
                orientation="horizontal")

        self.ax7.set_xlabel("Lithology")
        self.ax7.set_xlim(0, 5)
        self.ax7.set_xticks([])
        self.ax7.set_ylim(0, max_thickness)
        self.ax7.set_yticks(np.arange(0, max_thickness + step_depth, step_depth))
        self.ax7.margins(0.3, 0.0)
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        self.ax7.legend(handles=legend_lithology, loc="lower left", bbox_to_anchor=(-0.25, -0.125), shadow=False,
                        ncol=2, prop={'size': 7}, frameon=False)
        # #plt.tight_layout()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.subwindow_borehole_diagrams)
        self.canvas.get_tk_widget().grid(
            row=0, column=self.n_columns_setup + 1, rowspan=self.n_rows_diagram,
            columnspan=self.n_columns_diagram - self.n_columns_setup - 1, sticky="nesw")
        self.canvas.draw()