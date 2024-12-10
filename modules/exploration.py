#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# File:         exploration.py
# Description:  Contains all necessary functions that are related to mineral exploration feature
# Author:       Maximilian Beeskow
# Last updated: 10.12.2024
# License:      GPL v3.0

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# external packages
import tkinter as tk
from tkinter import ttk
import numpy as np
# internal packages
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

    def create_subwindow_borehole_data(self):
        ## Window Settings
        window_width = 900
        window_height = 750
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
                fg=self.colors["Navigation"]).create_button(text="Export data")
            btn_02.configure(activebackground=self.colors["Accent"])
            btn_03.configure(activebackground=self.colors["Accent"])
            btn_04.configure(activebackground=self.colors["Accent"])
            btn_05.configure(activebackground=self.colors["Accent"])
            btn_06.configure(activebackground=self.colors["Accent"])

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
                "Pelitic Greenschist", "Ortho-Amphibolite"]
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

            ## Treeviews
            categories = ["rho (kg/m\u00B3)", "phi (%)", "vP (m/s)", "vS (m/s)", "vP/vS (1)", "K (GPa)", "G (GPa)",
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

            for index, category in enumerate(categories):
                entries = [category]
                n_digits = 3
                var_entr_min = round(0, n_digits)
                var_entr_max = round(0, n_digits)
                var_entr_mean = round(0, n_digits)
                var_entr_error = round(0, n_digits)

                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])

                self.tv_lithology.insert("", tk.END, values=entries)

        if update:
            self.lbl_borehole_id.configure(text=self.current_borehole_id)
            self.lbl_unit_id.configure(text=self.current_unit_id)
            self.entr_top.configure(textvariable=self.dict_entr_top[self.current_borehole_id][self.current_unit_id])
            self.entr_bottom.configure(textvariable=self.dict_entr_bottom[self.current_borehole_id][
                self.current_unit_id])
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_igneous.set("Select igneous rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")

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
        categories = ["rho (kg/m\u00B3)", "phi (%)", "vP (m/s)", "vS (m/s)", "vP/vS (1)", "K (GPa)", "G (GPa)",
                      "E (GPa)", "nu (1)", "GR (API)", "PE (barns/e\u207B)"]
        categories_short = ["rho", "phi", "vP", "vS", "vPvS", "K", "G", "E", "v", "GR", "PE"]

        if focus == "sedimentary":
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            # Siliciclastic rocks
            if self.var_opt_sedimentary.get() == "Sandstone":
                data = SiliciclasticRocks(fluid="water", actualThickness=0).create_sandstone(
                    number=n_datapoints, porosity=[0.0, 0.3])
            elif self.var_opt_sedimentary.get() == "Conglomerate":
                data = SiliciclasticRocks(fluid="water", actualThickness=0).create_conglomerate(
                    number=n_datapoints, porosity=[0.0, 0.3])
            elif self.var_opt_sedimentary.get() == "Siltstone":
                data = SiliciclasticRocks(fluid="water", actualThickness=0).create_siltstone(
                    number=n_datapoints, porosity=[0.0, 0.1])
            elif self.var_opt_sedimentary.get() == "Mudstone":
                data = SiliciclasticRocks(fluid="water", actualThickness=0).create_mudstone_alt(
                    number=n_datapoints, porosity=[0.0, 0.1])
            elif self.var_opt_sedimentary.get() == "Shale":
                data = SiliciclasticRocks(fluid="water", actualThickness=0).create_shale_alt(
                    number=n_datapoints, porosity=[0.0, 0.1])
            elif self.var_opt_sedimentary.get() == "Greywacke (Huckenholz)":
                data = SiliciclasticRocks(fluid="water", actualThickness=0).create_greywacke_huckenholz(
                    rock="Greywacke", number=n_datapoints, porosity=[0.0, 0.1])
            # Carbonate rocks
            elif self.var_opt_sedimentary.get() == "Limestone":
                data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone(
                    number=n_datapoints,  porosity=[0.0, 0.4])
            elif self.var_opt_sedimentary.get() == "Dolostone":
                data = CarbonateRocks(fluid="water", actualThickness=0).create_dolostone(
                    number=n_datapoints, porosity=[0.0, 0.3])
            elif self.var_opt_sedimentary.get() == "Marl":
                data = CarbonateRocks(fluid="water", actualThickness=0).create_marl(
                    number=n_datapoints, porosity=[0.0, 0.3])

        elif focus == "plutonic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            if self.var_opt_plutonic.get() in [
                "Granite", "Granodiorite", "Tonalite", "Gabbro", "Norite", "Diorite", "Monzodiorite", "Monzogabbro",
                "Monzonite", "Syenite", "Granitoid", "Quarzolite", "Foid-bearing Syenite", "Foid-bearing Monzonite",
                "Foid-bearing Monzodiorite", "Foid-bearing Monzogabbro", "Foid Monzosyenite", "Foid Monzodiorite",
                "Foid Monzogabbro", "Foidolite"]:
                if self.var_opt_plutonic.get() in [
                    "Foid-bearing Syenite", "Foid-bearing Monzonite", "Foid-bearing Monzodiorite",
                    "Foid-bearing Monzogabbro", "Foid Monzosyenite", "Foid Monzodiorite", "Foid Monzogabbro",
                    "Foidolite"]:
                    data = Plutonic(
                        fluid="water", actualThickness=0, dict_output=True,
                        porosity=[0.0, 0.1]).create_plutonic_rock_streckeisen(
                        rock=self.var_opt_plutonic.get(), number=n_datapoints, porosity=[0.0, 0.1],
                        upper_streckeisen=False)
                else:
                    data = Plutonic(
                        fluid="water", actualThickness=0, dict_output=True,
                        porosity=[0.0, 0.1]).create_plutonic_rock_streckeisen(
                        rock=self.var_opt_plutonic.get(), number=n_datapoints, porosity=[0.0, 0.1])
            elif self.var_opt_igneous.get() in ["Rhyolite", "Dacite", "Trachyte", "Latite", "Andesite", "Basalt",
                "Foid-bearing Trachyte", "Foid-bearing Latite", "Foid-bearing Andesite", "Foid-bearing Basalt",
                "Phonolite", "Tephrite", "Foidite"]:
                if self.var_opt_igneous.get() in [
                    "Foid-bearing Trachyte", "Foid-bearing Latite", "Foid-bearing Andesite", "Foid-bearing Basalt",
                    "Phonolite", "Tephrite", "Foidite"]:
                    data = Volcanic(
                        fluid="water", actualThickness=0, dict_output=True,
                        porosity=[0.0, 0.1]).create_volcanic_rock_streckeisen(
                        rock=self.var_opt_igneous.get(), number=n_datapoints, upper_streckeisen=False)
                else:
                    data = Volcanic(
                        fluid="water", actualThickness=0, dict_output=True,
                        porosity=[0.0, 0.1]).create_volcanic_rock_streckeisen(
                        rock=self.var_opt_igneous.get(), number=n_datapoints)

        elif focus == "volcanic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            if self.var_opt_volcanic.get() in ["Rhyolite", "Dacite", "Trachyte", "Latite", "Andesite", "Basalt",
                "Foid-bearing Trachyte", "Foid-bearing Latite", "Foid-bearing Andesite", "Foid-bearing Basalt",
                "Phonolite", "Tephrite", "Foidite"]:
                if self.var_opt_volcanic.get() in [
                    "Foid-bearing Trachyte", "Foid-bearing Latite", "Foid-bearing Andesite", "Foid-bearing Basalt",
                    "Phonolite", "Tephrite", "Foidite"]:
                    data = Volcanic(
                        fluid="water", actualThickness=0, dict_output=True,
                        porosity=[0.0, 0.1]).create_volcanic_rock_streckeisen(
                        rock=self.var_opt_volcanic.get(), number=n_datapoints, upper_streckeisen=False)
                else:
                    data = Volcanic(
                        fluid="water", actualThickness=0, dict_output=True,
                        porosity=[0.0, 0.1]).create_volcanic_rock_streckeisen(
                        rock=self.var_opt_volcanic.get(), number=n_datapoints)

        elif focus == "ultramafic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")
            self.var_opt_ore.set("Select ore rock")

            if self.var_opt_ultramafic.get() in [
                "Orthopyroxenite", "Clinopyroxenite", "Dunite", "Harzburgite", "Wehrlite", "Websterite", "Lherzolite",
                "Olivine-Websterite", "Olivine-Orthopyroxenite", "Olivine-Clinopyroxenite", "Peridotite", "Pyroxenite"]:
                data = UltraMafic(
                    fluid="water", actualThickness=0, dict_output=True, porosity=[0.0, 0.1]).create_ultramafic_rock(
                    rock=self.var_opt_ultramafic.get(), number=n_datapoints)

        elif focus == "metamorphic":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_ore.set("Select ore rock")

            if self.var_opt_metamorphic.get() == "Felsic Granulite":
                data = GranuliteFacies(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_granulite(
                    number=n_datapoints, classification="felsic")
            elif self.var_opt_metamorphic.get() == "Mafic Granulite":
                data = GranuliteFacies(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_granulite(
                    number=n_datapoints, classification="mafic")
            elif self.var_opt_metamorphic.get() == "Basaltic Greenschist":
                data = GreenschistFacies(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenschist_basaltic_alt(
                    number=n_datapoints)
            elif self.var_opt_metamorphic.get() == "Ultramafic Greenschist":
                data = GreenschistFacies(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenschist_ultramafic_alt(
                    number=n_datapoints)
            elif self.var_opt_metamorphic.get() == "Pelitic Greenschist":
                data = GreenschistFacies(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_greenschist_pelitic_alt(
                    number=n_datapoints)
            # Amphibolite-Facies
            elif self.var_opt_metamorphic.get() == "Ortho-Amphibolite":
                data = AmphiboliteFacies(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_amphibolite_ortho(
                    number=n_datapoints)

        elif focus == "ore":
            self.var_opt_sedimentary.set("Select sedimentary rock")
            self.var_opt_plutonic.set("Select plutonic rock")
            self.var_opt_volcanic.set("Select volcanic rock")
            self.var_opt_ultramafic.set("Select ultramafic rock")
            self.var_opt_metamorphic.set("Select metamorphic rock")

            if self.var_opt_ore.get() in [
                "Itabirite", "Banded Iron Formation", "Compact Hematite", "Friable Hematite", "Goethite Hematite",
                "Al-rich Itabirite", "Compact Quartz Itabirite", "Friable Quartz Itabirite", "Goethite Itabirite"]:
                data = OreRocks(
                    fluid="water", actual_thickness=0, porosity=[0.0, 0.1]).create_siliciclastic_itabirite(
                    rock=self.var_opt_ore.get(), number=n_datapoints)

        self.fill_lithology_table(dataset=data, categories_long=categories, categories_short=categories_short)

    def fill_lithology_table(self, dataset, categories_long, categories_short):
        for index, category in enumerate(categories_long):
            entries = [category]

            if categories_short[index] == "vPvS":
                category_short = "vP/vS"
            elif categories_short[index] == "v":
                category_short = "nu"
            else:
                category_short = categories_short[index]

            n_digits = 3

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

                var_entr_min = round(var_factor*min(data_values), n_digits)
                var_entr_max = round(var_factor*max(data_values), n_digits)
                var_entr_mean = round(var_factor*np.mean(data_values), n_digits)
                var_entr_error = round(var_factor*np.std(data_values, ddof=1), n_digits)

                entries.extend([var_entr_min, var_entr_max, var_entr_mean, var_entr_error])
                self.tv_lithology.insert("", tk.END, values=entries)