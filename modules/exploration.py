#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# File:         exploration.py
# Description:  Contains all necessary functions that are related to mineral exploration
# Author:       Maximilian Beeskow
# Last updated: 01.02.2024
# License:      GPL v3.0

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# external packages
import tkinter as tk
from modules.gui_elements import SimpleElements

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -


class ExplorationInterface:

    def __init__(self, parent):
        self.parent = parent
        self.colors = {
            "Background": "#ECEBEA", "Navigation": "#252422", "Accent": "#EB5E28", "Option": "#CCC5B9",
            "White": "#FFFFFF", "Black": "#000000", "Accent Blue": "#118AB2"}

        # Variables
        self.var_rb_setup = tk.IntVar()
        self.var_rb_setup.set(0)

        self.var_entr_boreholes = tk.IntVar()
        self.var_entr_boreholes.set(2)
        self.var_entr_units = tk.IntVar()
        self.var_entr_units.set(4)
        self.var_entr_maximum_depth = tk.IntVar()
        self.var_entr_maximum_depth.set(100)
        self.var_entr_stepsize = tk.StringVar()
        self.var_entr_stepsize.set("0.5")

    def create_subwindow_borehole_data(self):
        ## Window Settings
        window_width = 900
        window_height = 600
        var_geometry = str(window_width) + "x" + str(window_height) + "+" + str(0) + "+" + str(0)

        row_min = 25
        n_rows = int(window_height/row_min)
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
        for y in range(n_rows):
            tk.Grid.rowconfigure(self.subwindow_borehole_data, y, weight=1)

        # Rows
        for i in range(0, n_rows):
            self.subwindow_borehole_data.grid_rowconfigure(i, minsize=row_min)
        # Columns
        for i in range(0, n_columns):
            self.subwindow_borehole_data.grid_columnconfigure(i, minsize=column_min)

        self.n_columns_setup = 10

        ## Frames
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=0, column_id=self.n_columns_setup + 1, n_rows=n_rows,
            n_columns=n_columns - self.n_columns_setup - 1, bg=self.colors["Background"],
            fg=self.colors["Background"]).create_frame()

        ## Labels
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=0, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Settings", relief=tk.FLAT, font_option="sans 16 bold", anchor_option=tk.W)

        ## Radiobuttons
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=1, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_radiobutton(
            text="Create new dataset", var_rb=self.var_rb_setup, value_rb=0, color_bg=self.colors["Navigation"],
            command=self.select_dataset_setup)
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=2, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_radiobutton(
            text="Load dataset", var_rb=self.var_rb_setup, value_rb=1, color_bg=self.colors["Navigation"],
            command=self.select_dataset_setup)
        
        ## Initialization
        self.select_dataset_setup()

    def select_dataset_setup(self):
        if self.var_rb_setup.get() == 0:
            ## Labels
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=4, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Number of boreholes (#)", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=6, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Number of units (#)", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=8, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Maximum depth (m)", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=10, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
                text="Stepsize (m)", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)

            ## Entries
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=5, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_boreholes)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=7, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_units)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=9, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_maximum_depth)
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=11, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Background"],
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_stepsize)

            ## Buttons
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=13, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"],
                fg=self.colors["Navigation"]).create_button(text="Update settings", command=self.update_settings)

        elif self.var_rb_setup.get() == 1:
            print("Load dataset!")

    def update_settings(self):
        ## Labels
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=15, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Current borehole", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=17, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Current unit", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)