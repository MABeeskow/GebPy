#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# File:         exploration.py
# Description:  Contains all necessary functions that are related to mineral exploration
# Author:       Maximilian Beeskow
# Last updated: 02.02.2024
# License:      GPL v3.0

# ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- -

# external packages
import tkinter as tk
import numpy as np
# internal packages
from modules.gui_elements import SimpleElements

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

        self.n_columns_setup = 11

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
        if self.initialization == False:
            self.select_dataset_setup()
            self.update_settings(initialization=True)
            self.initialization = True

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
                text="Parts per unit (#)", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)

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
                fg=self.colors["Navigation"]).create_entry(var_entr=self.var_entr_parts)

            ## Buttons
            SimpleElements(
                parent=self.subwindow_borehole_data, row_id=13, column_id=1, n_rows=1,
                n_columns=self.n_columns_setup - 1, bg=self.colors["Navigation"],
                fg=self.colors["Navigation"]).create_button(
                text="Update settings", command=lambda initialization=True: self.update_settings(initialization))

        elif self.var_rb_setup.get() == 1:
            print("Load dataset!")

    def update_settings(self, initialization=False):
        ## Helper
        if initialization:
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

        ## Labels
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=15, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Current borehole", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)
        self.lbl_borehole_id = SimpleElements(
            parent=self.subwindow_borehole_data, row_id=16, column_id=1, n_rows=1, n_columns=2,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text=self.current_borehole_id, relief=tk.FLAT, font_option="sans 14 bold")
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=17, column_id=1, n_rows=1, n_columns=self.n_columns_setup - 1,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Current unit", relief=tk.FLAT, font_option="sans 14 bold", anchor_option=tk.W)
        self.lbl_unit_id = SimpleElements(
            parent=self.subwindow_borehole_data, row_id=18, column_id=1, n_rows=1, n_columns=2,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text=self.current_unit_id, relief=tk.FLAT, font_option="sans 14 bold")
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=19, column_id=1, n_rows=1, n_columns=5,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Lithology", relief=tk.FLAT, font_option="sans 12 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=20, column_id=1, n_rows=1, n_columns=6,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Top depth (m)", relief=tk.FLAT, font_option="sans 12 bold", anchor_option=tk.W)
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=21, column_id=1, n_rows=1, n_columns=6,
            bg=self.colors["Navigation"], fg=self.colors["White"]).create_label(
            text="Bottom depth (m)", relief=tk.FLAT, font_option="sans 12 bold", anchor_option=tk.W)

        ## Buttons
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=16, column_id=3, n_rows=1, n_columns=4,
            bg=self.colors["Navigation"], fg=self.colors["Navigation"]).create_button(
            text="Previous", command=lambda mode="previous": self.change_borehole(mode))
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=16, column_id=7, n_rows=1, n_columns=4,
            bg=self.colors["Navigation"], fg=self.colors["Navigation"]).create_button(
            text="Next", command=lambda mode="next": self.change_borehole(mode))
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=18, column_id=3, n_rows=1, n_columns=4,
            bg=self.colors["Navigation"], fg=self.colors["Navigation"]).create_button(
            text="Previous", command=lambda mode="previous": self.change_unit(mode))
        SimpleElements(
            parent=self.subwindow_borehole_data, row_id=18, column_id=7, n_rows=1, n_columns=4,
            bg=self.colors["Navigation"], fg=self.colors["Navigation"]).create_button(
            text="Next", command=lambda mode="next": self.change_unit(mode))

        ## Entries
        self.entr_top = SimpleElements(
            parent=self.subwindow_borehole_data, row_id=20, column_id=7, n_rows=1, n_columns=4,
            bg=self.colors["Background"], fg=self.colors["Navigation"]).create_entry(
            var_entr=self.dict_entr_top[self.current_borehole_id][self.current_unit_id])
        self.entr_bottom = SimpleElements(
            parent=self.subwindow_borehole_data, row_id=21, column_id=7, n_rows=1, n_columns=4,
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

    def change_borehole(self, mode):
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
        self.update_settings()

    def change_unit(self, mode):
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