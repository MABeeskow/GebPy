#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gui.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.05.2023

#-----------------------------------------------

## MODULES
import os, sys, datetime
import tkinter as tk
import csv, re
from modules.gui_elements import SimpleElements as SE
from modules.sulfates import Sulfates
from modules.oxides import Oxides, RutileGroup, PericlaseGroup, WulfeniteGroup
from modules.sulfides import Sulfides
from modules.carbonates import Carbonates, CarbonateRocks
from modules.halides import Halides
from modules.silicates import Tectosilicates, Phyllosilicates, Nesosilicates, Sorosilicates, Inosilicates, Cyclosilicates
from modules.phospides import Phospides
from modules.phosphates import Phosphates
from modules.organics import Organics
from modules.minerals import feldspars
from modules.sedimentary_rocks import SedimentaryRocks
from modules.siliciclastics import sandstone, shale
from modules.metamorphics import MetamorphicRocks
from modules.carbonates import limestone, dolomite
from modules.sequences import SedimentaryBasin
from modules.ore import Ores
from modules.igneous import Plutonic, Volcanic, Pyroclastic
from modules.evaporites import Evaporites
from modules.sequences import DataProcessing as DP
from modules.geophysics import Elasticity as elast
from modules.series import Zechstein, Muschelkalk
from modules import fluids
import numpy as np
import random as rd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import NullFormatter
import matplotlib.patches as mpatches
import time
from modules.mineralogy import Mineralogy
import seaborn as sns

## GUI
class GebPyGUI(tk.Frame):
    #
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        #
        ## CONSTANTS
        self.color_menu = "#264653"
        self.color_border = "#7C9097"
        self.color_bg = "#E9ECED"
        self.color_fg_dark = "black"
        self.color_fg_light = "white"
        self.color_accent_01 = "#E76F51"    # Orange
        self.color_accent_02 = "#F0A794"    # Orange light
        self.color_accent_03 = "#E9C46A"    # Yellow
        self.color_accent_04 = "#F3DEAD"    # Yellow light
        #
        self.parent = parent
        self.parent.title("GebPy")
        self.parent.geometry("1700x960+10+10")
        self.parent.resizable(False, False)
        self.parent["bg"] = self.color_bg
        #
        self.lbl_w = {}
        self.lbl_w["physics"] = []
        self.lbl_w["chemistry"] = []
        self.entr_w = {}
        self.entr_w["physics"] = []
        self.entr_w["chemistry"] = []
        self.entr_w["custom"] = {}
        self.gui_elements = []
        self.container_gui = {}
        self.container_gui_specific = {}
        self.categories_sections = ["MINERALOGY", "PETROLOGY", "STRATIGRAPHY"]
        self.sub_categories_gui = ["LABEL", "BUTTON", "RADIOBUTTON", "CHECKBOX", "OPTION MENU", "ENTRY", "CANVAS_HISTO",
                                   "CANVAS_SCATTER", "CANVAS_GEOCHEM"]
        self.container_var = {}
        #
        for category in self.categories_sections:
            self.container_gui[category] = {}
            self.container_gui_specific[category] = {}
            self.container_var[category] = {}
            for sub_category in self.sub_categories_gui:
                self.container_gui[category][sub_category] = []
                self.container_gui_specific[category][sub_category] = []
                if sub_category == "ENTRY":
                    self.container_var[category][sub_category] = {}
                    self.container_var[category][sub_category]["MIN"] = {}
                    self.container_var[category][sub_category]["MAX"] = {}
                    self.container_var[category][sub_category]["MEAN"] = {}
                    self.container_var[category][sub_category]["STD"] = {}
                elif sub_category == "OPTION MENU":
                    self.container_var[category][sub_category] = {}
                else:
                    self.container_var[category][sub_category] = []
        #
        self.exp_data = []
        self.filename = []
        #
        for x in range(19):
            tk.Grid.columnconfigure(self.parent, x, weight=1)
        for y in range(48):
            tk.Grid.rowconfigure(self.parent, y, weight=1)
        #
        # Rows
        for i in range(0, 48):
            self.parent.grid_rowconfigure(i, minsize=20)
        # Columns
        self.parent.grid_columnconfigure(0, minsize=160)
        self.parent.grid_columnconfigure(1, minsize=110)
        self.parent.grid_columnconfigure(2, minsize=15)
        self.parent.grid_columnconfigure(3, minsize=125)
        for i in range(4, 8):
            self.parent.grid_columnconfigure(i, minsize=90)
        self.parent.grid_columnconfigure(8, minsize=15)
        for i in range(9, 18):
            self.parent.grid_columnconfigure(i, minsize=100)
        self.parent.grid_columnconfigure(18, minsize=15)
        #
        menu_frame = tk.Frame(self.parent, bg=self.color_menu, borderwidth=0, highlightthickness=0)
        menu_frame.grid(row=0, column=0, rowspan=48, columnspan=2, sticky="nesw")
        menu_frame2 = tk.Frame(self.parent, bg="#7C9097", borderwidth=0, highlightthickness=0)
        menu_frame2.grid(row=0, column=2, rowspan=48, sticky="nesw")
        #
        ## Logo
        gebpy_logo = tk.PhotoImage(file="../documents/readme_images/GebPy_Logo.png")
        gebpy_logo = gebpy_logo.subsample(6, 6)
        img = tk.Label(self.parent, image=gebpy_logo, bg=self.color_menu)
        img.image = gebpy_logo
        img.grid(row=0, column=0, rowspan=4, columnspan=2, sticky="nesw")
        #
        self.var_rb_main = tk.IntVar()
        self.var_rb_main.set(0)
        self.var_rb_mineralogy_mineral = tk.IntVar()
        self.var_rb_mineralogy_mineral.set(0)
        self.var_rb_mineralogy_plot = tk.IntVar()
        self.var_rb_mineralogy_plot.set(0)
        self.var_rb_mineralogy_traces = tk.IntVar()
        self.var_rb_mineralogy_traces.set(0)
        self.var_opt_mineral = tk.StringVar()
        self.var_opt_mineral.set("Select Mineral")
        self.var_opt_mineral_group = tk.StringVar()
        self.var_opt_mineral_group.set("Select Mineral Group")
        #
        ## Radiobuttons
        rb_01 = SE(
            parent=self.parent, row_id=4, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
            fg=self.color_bg).create_radiobutton(
            var_rb=self.var_rb_main, value_rb=0, text="Mineralogy", color_bg=self.color_menu,
            command=self.change_rb_main)
        rb_02 = SE(
            parent=self.parent, row_id=4, column_id=1, n_rows=1, n_columns=1, bg=self.color_menu,
            fg=self.color_bg).create_radiobutton(
            var_rb=self.var_rb_main, value_rb=1, text="Petrology", color_bg=self.color_menu,
            command=self.change_rb_main)
        rb_03 = SE(
            parent=self.parent, row_id=5, column_id=0, n_rows=1, n_columns=2, bg=self.color_menu,
            fg=self.color_bg).create_radiobutton(
            var_rb=self.var_rb_main, value_rb=2, text="Sequence Stratigraphy",
            color_bg=self.color_menu, command=self.change_rb_main)
        #
        self.change_rb_main(first_start=True)
        #
        btn_advstat = SE(parent=self.parent, row_id=42, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_03,
                         fg="black").create_button(text="Advanced Statistics")
        btn_exp = SE(parent=self.parent, row_id=44, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_03,
                     fg="black").create_button(text="Export Data",
                                               command=lambda var_btn="export data": self.pressed_button(var_btn))
        btn_exit = SE(parent=self.parent, row_id=46, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
                      fg=self.color_fg_dark).create_button(text="Quit", command=self.parent.quit)
    #
    def change_rb_main(self, first_start=False):
        if first_start == False:
            if self.var_rb_main.get() == 0:
                ## MINERALOGY
                ## CLEANING
                for category in self.categories_sections:
                    if category != "MINERALOGY":
                        for sub_category in self.sub_categories_gui:
                            for gui_item in self.container_gui[category][sub_category]:
                                if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                                    gui_item.grid_remove()
                                else:
                                    gui_item.get_tk_widget().grid_remove()
                #
                try:
                    self.leg.remove()
                except:
                    pass
                #
                ## RECONSTRUCTION
                for sub_category in self.sub_categories_gui:
                    for gui_item in self.container_gui["MINERALOGY"][sub_category]:
                        if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                            gui_item.grid()
                        else:
                            gui_item.get_tk_widget().grid()
                #
            elif self.var_rb_main.get() == 1:
                # PETROLOGY
                ## CLEANING
                for category in self.categories_sections:
                    if category != "PETROLOGY":
                        for sub_category in self.sub_categories_gui:
                            for gui_item in self.container_gui[category][sub_category]:
                                if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                                    gui_item.grid_remove()
                                else:
                                    gui_item.get_tk_widget().grid_remove()
                #
                try:
                    self.leg.remove()
                except:
                    pass
                #
                ## RECONSTRUCTION
                for sub_category in self.sub_categories_gui:
                    for gui_item in self.container_gui["PETROLOGY"][sub_category]:
                        if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                            gui_item.grid()
                        else:
                            gui_item.get_tk_widget().grid()
                #
                if len(self.container_gui["PETROLOGY"]["LABEL"]) == 0:
                    lb_01 = SE(
                        parent=self.parent, row_id=6, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
                        fg=self.color_fg_dark).create_label(text="Petrology", relief=tk.RAISED)
                    #
                    self.container_gui["PETROLOGY"]["LABEL"].append(lb_01)
                    #
                    var_opt_1_0 = tk.StringVar()
                    opt_list_1_0 = ["Siliciclastic Rocks", "Carbonate Rocks", "Plutonic Rocks", "Volcanic Rocks",
                                    "Pyroclastic Rocks", "Metamorphic Rocks", "Evaporite Rocks", "Ore Rocks"]
                    opt_list_1_0.sort()
                    self.opt_rocktype = SE(parent=self.parent, row_id=8, column_id=0, n_rows=2, n_columns=2,
                                           bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                        var_opt=var_opt_1_0, var_opt_set="Select Rock Type", opt_list=opt_list_1_0,
                        command=lambda var_opt=var_opt_1_0: self.select_opt(var_opt))
                    #
                    self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_rocktype)
                    #
                    btn_custseq = SE(
                        parent=self.parent, row_id=12, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                        fg=self.color_fg_dark).create_button(
                        text="Create Custom Rock", command=lambda var_btn="custom rock": self.pressed_button(var_btn))
                    #
                    self.container_gui["PETROLOGY"]["BUTTON"].append(btn_custseq)
                    #
            elif self.var_rb_main.get() == 2:
                ## STRATIGRAPHY
                ## CLEANING
                for category in self.categories_sections:
                    if category != "STRATIGRAPHY":
                        for sub_category in self.sub_categories_gui:
                            for gui_item in self.container_gui[category][sub_category]:
                                if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                                    gui_item.grid_remove()
                                else:
                                    gui_item.get_tk_widget().grid_remove()
                #
                try:
                    self.leg.remove()
                except:
                    pass
                #
                ## RECONSTRUCTION
                for sub_category in self.sub_categories_gui:
                    for gui_item in self.container_gui["STRATIGRAPHY"][sub_category]:
                        if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                            gui_item.grid()
                        else:
                            gui_item.get_tk_widget().grid()
                #
                if len(self.container_gui["STRATIGRAPHY"]["LABEL"]) == 0:
                    lb_01 = SE(
                        parent=self.parent, row_id=6, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
                        fg=self.color_fg_dark).create_label(text="Sequence Stratigraphy", relief=tk.RAISED)
                    #
                    self.container_gui["STRATIGRAPHY"]["LABEL"].append(lb_01)
                    #
                    var_opt_2_0 = tk.StringVar()
                    opt_list_2_0 = ["Zechstein", "Muschelkalk"]
                    self.opt_realseq = SE(parent=self.parent, row_id=8, column_id=0, n_rows=2, n_columns=2,
                                          bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                        var_opt=var_opt_2_0, var_opt_set="Select Real Sequences", opt_list=opt_list_2_0,
                        command=lambda var_opt=var_opt_2_0: self.select_opt(var_opt))
                    #
                    self.container_gui["STRATIGRAPHY"]["OPTION MENU"].append(self.opt_realseq)
                    #
                    btn_randseq = SE(
                        parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                        fg=self.color_fg_dark).create_button(
                        text="Create Random Sequence", command=lambda var_btn="random": self.pressed_button(var_btn))
                    #
                    self.container_gui["STRATIGRAPHY"]["BUTTON"].append(btn_randseq)
                    #
        else:
            ## MINERALOGY
            #
            self.container_var["MINERALOGY"]["OPTION MENU"]["GEOCHEM"] = tk.StringVar()
            #
            lb_01 = SE(parent=self.parent, row_id=6, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
                       fg=self.color_fg_dark).create_label(text="Mineralogy", relief=tk.RAISED)
            #
            self.container_gui["MINERALOGY"]["LABEL"].append(lb_01)
            #
            var_opt_0_0 = tk.StringVar()
            opt_list_0_0 = ["Oxides", "Sulfides", "Carbonates", "Halides", "Tectosilicates", "Phyllosilicates",
                            "Sulfates", "Nesosilicates", "Sorosilicates", "Inosilicates", "Phosphates", "Phosphides",
                            "Organics", "Cyclosilicates"]
            opt_list_0_0.sort()
            self.opt_mingroup = SE(parent=self.parent, row_id=8, column_id=0, n_rows=2, n_columns=2,
                                   bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=self.var_opt_mineral_group, var_opt_set="Select Mineral Group", opt_list=opt_list_0_0,
                command=lambda var_opt=self.var_opt_mineral_group: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_mingroup)
            #
            ## Labels
            lb_trace = SE(parent=self.parent, row_id=14, column_id=0, n_rows=1, n_columns=2, bg=self.color_accent_03,
                          fg=self.color_fg_dark).create_label(text="Trace Elements", relief=tk.RAISED)
            #
            self.container_gui["MINERALOGY"]["LABEL"].extend([lb_trace])
            #
            ## Radiobuttons
            rb_trace_01 = SE(
                parent=self.parent, row_id=15, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
                fg=self.color_bg).create_radiobutton(
                var_rb=self.var_rb_mineralogy_traces, value_rb=0, text="Without Traces",
                color_bg=self.color_menu, command=self.change_rb_trace)
            rb_trace_02 = SE(
                parent=self.parent, row_id=16, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
                fg=self.color_bg).create_radiobutton(
                var_rb=self.var_rb_mineralogy_traces, value_rb=1, text="With Traces",
                color_bg=self.color_menu, command=self.change_rb_trace)
            #
            self.container_gui["MINERALOGY"]["RADIOBUTTON"].extend([rb_trace_01, rb_trace_02])
            #
    #
    def pressed_button(self, var_btn):
        if var_btn == "random":
            Subsurface(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                       color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_btn, lbl_w=self.lbl_w,
                       entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_btn == "custom rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock="Custom", lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_btn == "export data":
            self.export_data(var_exp=self.exp_data)
        else:
            Subsurface(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                       color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_btn, lbl_w=self.lbl_w,
                       entr_w=self.entr_w, gui_elements=self.gui_elements)
    #
    def change_rb_comparison_plot(self, dataset):
        n_plot_histo = len(self.container_gui["MINERALOGY"]["CANVAS_HISTO"])
        n_plot_scatter = len(self.container_gui["MINERALOGY"]["CANVAS_SCATTER"])
        n_plot_geochem = len(self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"])
        #
        if self.var_rb_mineralogy_plot.get() == 0:      # Histogram
            if (n_plot_histo + n_plot_scatter + n_plot_geochem) == 0:
                self.plot_histogram_comparison(dataset=dataset)
            elif n_plot_histo > 0 and (n_plot_scatter + n_plot_geochem) == 0:
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid()
            elif n_plot_histo > 0 and n_plot_scatter > 0 and n_plot_geochem == 0:
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid()
            elif n_plot_histo > 0 and n_plot_scatter == 0 and n_plot_geochem > 0:
                self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid()
            elif n_plot_histo > 0 and n_plot_scatter > 0 and n_plot_geochem > 0:
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid()
            elif n_plot_histo == 0 and (n_plot_scatter + n_plot_geochem) > 0:
                try:
                    self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid_remove()
                except:
                    self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"][-1].get_tk_widget().grid_remove()
                self.plot_histogram_comparison(dataset=dataset)
        #
        elif self.var_rb_mineralogy_plot.get() == 1:    # Scatter
            if (n_plot_histo + n_plot_scatter + n_plot_geochem) == 0:
                self.plot_scatter_comparison(dataset=dataset)
            elif n_plot_scatter > 0 and (n_plot_histo + n_plot_geochem) == 0:
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid()
            elif n_plot_scatter > 0 and n_plot_histo > 0 and n_plot_geochem == 0:
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid()
            elif n_plot_scatter > 0 and n_plot_histo > 0 and n_plot_geochem > 0:
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid()
            elif n_plot_scatter == 0 and (n_plot_histo + n_plot_geochem) > 0:
                try:
                    self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid_remove()
                except:
                    self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"][-1].get_tk_widget().grid_remove()
                self.plot_scatter_comparison(dataset=dataset)
        #
        elif self.var_rb_mineralogy_plot.get() == 2:  # Geochemistry
            if (n_plot_histo + n_plot_scatter + n_plot_geochem) == 0:
                self.plot_geochem_comparison()
            elif n_plot_geochem == 0 and n_plot_histo > 0 and n_plot_scatter == 0:
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid_remove()
                self.plot_geochem_comparison()
            elif n_plot_geochem == 0 and n_plot_histo > 0 and n_plot_scatter > 0:
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid_remove()
                self.plot_geochem_comparison()
            elif n_plot_geochem > 0 and n_plot_histo > 0 and n_plot_scatter > 0:
                self.container_gui["MINERALOGY"]["CANVAS_HISTO"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_SCATTER"][-1].get_tk_widget().grid_remove()
                self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"][-1].get_tk_widget().grid()
    #
    def change_opt_comparison_geochem(self, event):
        self.plot_geochem_comparison()
    #
    def change_rb_comparison_key(self, dataset):
        dataset_keys = list(dataset.keys())
        #
        lbl_stat = SE(
            parent=self.parent, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
            fg=self.color_fg_dark).create_label(
            text="Statistics - "+str(dataset_keys[self.var_rb_mineralogy_mineral.get()]), relief=tk.RAISED)
        lb_param = SE(parent=self.parent, row_id=2, column_id=3, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Parameter", relief=tk.RAISED)
        lb_min = SE(parent=self.parent, row_id=2, column_id=4, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Minimum", relief=tk.RAISED)
        lb_max = SE(parent=self.parent, row_id=2, column_id=5, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Maximum", relief=tk.RAISED)
        lb_mean = SE(parent=self.parent, row_id=2, column_id=6, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Mean", relief=tk.RAISED)
        lb_std = SE(parent=self.parent, row_id=2, column_id=7, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Standard\n Deviation", relief=tk.RAISED)
        #
        self.container_gui["MINERALOGY"]["LABEL"].extend([lbl_stat, lb_param, lb_min, lb_max, lb_mean, lb_std])
        #
        self.fill_table_comparison(dataset=dataset[dataset_keys[self.var_rb_mineralogy_mineral.get()]])
    #
    def select_opt(self, var_opt):
        # Minerals
        ## OXIDES
        if var_opt in ["Magnetite", "Hematite", "Aluminium Spinels", "Ilmenite", "Cassiterite", "Chromite",
                       "Corundum", "Rutile", "Pyrolusite", "Magnesiochromite", "Zincochromite", "Chromium Spinels",
                       "Cuprospinel", "Jacobsite", "Magnesioferrite", "Trevorite", "Franklinite", "Ulvöspinel",
                       "Iron Spinels", "Uraninite", "Litharge", "Massicot", "Minium", "Plattnerite", "Scrutinyite",
                       "Zincite", "Columbite", "Tantalite", "Coltan", "Hematite-Group", "Rutile-Group",
                       "Periclase-Group", "Crocoite", "Wulfenite", "Wulfenite-Group", "Goethite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        #
        elif var_opt in ["Quartz"]:
            # #
            # ## Labels
            # lb_trace = SE(parent=self.parent, row_id=14, column_id=0, n_rows=1, n_columns=2, bg=self.color_accent_03,
            #               fg=self.color_fg_dark).create_label(text="Trace Elements", relief=tk.RAISED)
            # #
            # self.container_gui["MINERALOGY"]["LABEL"].extend([lb_trace])
            # #
            # ## Radiobuttons
            # rb_trace_01 = SE(
            #     parent=self.parent, row_id=15, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
            #     fg=self.color_bg).create_radiobutton(
            #     var_rb=self.var_rb_mineralogy_traces, value_rb=0, text="Without Traces",
            #     color_bg=self.color_menu, command=self.change_rb_trace)
            # rb_trace_02 = SE(
            #     parent=self.parent, row_id=16, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
            #     fg=self.color_bg).create_radiobutton(
            #     var_rb=self.var_rb_mineralogy_traces, value_rb=1, text="With Traces",
            #     color_bg=self.color_menu, command=self.change_rb_trace)
            # #
            # self.container_gui["MINERALOGY"]["RADIOBUTTON"].extend([rb_trace_01, rb_trace_02])
            # #
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
            #
        elif var_opt in ["Spinel Group", "Hematite Group", "Rutile Group", "Periclase Group", "Wulfenite Group"]:
            #
            ## CLEANING
            try:
                list_canvas = ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]
                for canvas in list_canvas:
                    if len(self.container_gui["MINERALOGY"][canvas]) > 0:
                        for canvas_item in self.container_gui["MINERALOGY"][canvas]:
                            canvas_item.get_tk_widget().grid_remove()
                        self.container_gui["MINERALOGY"][canvas].clear()
            except:
                pass
            #
            self.data_minerals = Mineralogy(keyword=var_opt).compare_minerals(number=100)
            #
            index = 0
            for key, value in self.data_minerals.items():
                rb_key = SE(
                    parent=self.parent, row_id=15 + index, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
                    fg=self.color_bg).create_radiobutton(
                    var_rb=self.var_rb_mineralogy_mineral, value_rb=index, text=key,
                    color_bg=self.color_menu, command=lambda dataset=self.data_minerals: self.change_rb_comparison_key(dataset))
                #
                if index == 0:
                    start_key = key
                    self.list_elements_geochem = list(value["chemistry"].keys())
                    self.container_var["MINERALOGY"]["OPTION MENU"]["GEOCHEM"].set(self.list_elements_geochem[0])
                    #
                    opt_geochem = SE(
                        parent=self.parent, row_id=19, column_id=1, n_rows=1, n_columns=1, bg=self.color_bg,
                        fg="black").create_option_menu(
                        var_opt=self.container_var["MINERALOGY"]["OPTION MENU"]["GEOCHEM"],
                        var_opt_set=self.container_var["MINERALOGY"]["OPTION MENU"]["GEOCHEM"].get(),
                        opt_list=self.list_elements_geochem,
                        active_bg=self.color_bg, command=lambda event: self.change_opt_comparison_geochem(event))
                    #
                    self.container_gui["MINERALOGY"]["OPTION MENU"].append(opt_geochem)
                #
                self.container_gui["MINERALOGY"]["RADIOBUTTON"].append(rb_key)
                #
                index += 1
            #
            plot_categories = ["Histogram", "Scatter", "Geochemistry"]
            for index, plot_category in enumerate(plot_categories):
                rb_plot = SE(
                    parent=self.parent, row_id=15 + index, column_id=1, n_rows=1, n_columns=1, bg=self.color_menu,
                    fg=self.color_bg).create_radiobutton(
                    var_rb=self.var_rb_mineralogy_plot, value_rb=index, text=plot_category,
                    color_bg=self.color_menu,
                    command=lambda dataset=self.data_minerals: self.change_rb_comparison_plot(dataset))
                #
                self.container_gui["MINERALOGY"]["RADIOBUTTON"].append(rb_plot)
                #
            #
            table_categories = ["Molar Mass", "Molar Volume", "Density", "vP", "vS", "vP/vS", "K", "G", "E", "Poisson",
                                "GR", "PE"]
            for index, table_category in enumerate(table_categories):
                lb_01 = SE(
                    parent=self.parent, row_id=4 + 2*index, column_id=3, n_rows=2, n_columns=1, bg=self.color_bg,
                    fg=self.color_fg_dark).create_label(text=table_category, relief=tk.RAISED)
                #
                self.container_gui["MINERALOGY"]["LABEL"].append(lb_01)
                #
                self.container_var["MINERALOGY"]["ENTRY"]["MIN"][table_category] = tk.StringVar()
                self.container_var["MINERALOGY"]["ENTRY"]["MAX"][table_category] = tk.StringVar()
                self.container_var["MINERALOGY"]["ENTRY"]["MEAN"][table_category] = tk.StringVar()
                self.container_var["MINERALOGY"]["ENTRY"]["STD"][table_category] = tk.StringVar()
                #
                entr_min = SE(
                    parent=self.parent, row_id=4 + 2*index, column_id=4, bg=self.color_bg, n_rows=2,
                    fg=self.color_fg_dark).create_entry(
                    var_entr=self.container_var["MINERALOGY"]["ENTRY"]["MIN"][table_category], var_entr_set=0.0)
                entr_max = SE(
                    parent=self.parent, row_id=4 + 2 * index, column_id=5, bg=self.color_bg, n_rows=2,
                    fg=self.color_fg_dark).create_entry(
                    var_entr=self.container_var["MINERALOGY"]["ENTRY"]["MAX"][table_category], var_entr_set=0.0)
                entr_mean = SE(
                    parent=self.parent, row_id=4 + 2 * index, column_id=6, bg=self.color_bg, n_rows=2,
                    fg=self.color_fg_dark).create_entry(
                    var_entr=self.container_var["MINERALOGY"]["ENTRY"]["MEAN"][table_category], var_entr_set=0.0)
                entr_std = SE(
                    parent=self.parent, row_id=4 + 2 * index, column_id=7, bg=self.color_bg, n_rows=2,
                    fg=self.color_fg_dark).create_entry(
                    var_entr=self.container_var["MINERALOGY"]["ENTRY"]["STD"][table_category], var_entr_set=0.0)
                #
                self.container_gui["MINERALOGY"]["ENTRY"].extend([entr_min, entr_max, entr_mean, entr_std])
                #
            #
            self.plot_histogram_comparison(dataset=self.data_minerals)
            self.fill_table_comparison(dataset=self.data_minerals[start_key])
            #
            lbl_stat = SE(
                parent=self.parent, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
                fg=self.color_fg_dark).create_label(
                text="Statistics - " + str(start_key), relief=tk.RAISED)
            lb_param = SE(parent=self.parent, row_id=2, column_id=3, n_rows=2, bg=self.color_bg,
                          fg=self.color_fg_dark).create_label(text="Parameter", relief=tk.RAISED)
            lb_min = SE(parent=self.parent, row_id=2, column_id=4, n_rows=2, bg=self.color_bg,
                        fg=self.color_fg_dark).create_label(text="Minimum", relief=tk.RAISED)
            lb_max = SE(parent=self.parent, row_id=2, column_id=5, n_rows=2, bg=self.color_bg,
                        fg=self.color_fg_dark).create_label(text="Maximum", relief=tk.RAISED)
            lb_mean = SE(parent=self.parent, row_id=2, column_id=6, n_rows=2, bg=self.color_bg,
                         fg=self.color_fg_dark).create_label(text="Mean", relief=tk.RAISED)
            lb_std = SE(parent=self.parent, row_id=2, column_id=7, n_rows=2, bg=self.color_bg,
                        fg=self.color_fg_dark).create_label(text="Standard\n Deviation", relief=tk.RAISED)
            #
            self.container_gui["MINERALOGY"]["LABEL"].extend([lbl_stat, lb_param, lb_min, lb_max, lb_mean, lb_std])
            #
        #
        ## SULFIDES
        elif var_opt in ["Pyrite", "Chalcopyrite", "Galena", "Acanthite", "Chalcocite", "Bornite", "Sphalerite",
                         "Pyrrhotite", "Millerite", "Pentlandite", "Covellite", "Cinnabar", "Realgar", "Orpiment",
                         "Stibnite", "Marcasite", "Molybdenite", "Fahlore", "Chalcopyrite-Group", "Gallite",
                         "Roquesite", "Lenaite", "Laforetite", "Vaesite", "Cattierite", "Pyrite-Group"]:
            #
            ## Labels
            lb_trace = SE(parent=self.parent, row_id=14, column_id=0, n_rows=1, n_columns=2, bg=self.color_accent_03,
                          fg=self.color_fg_dark).create_label(text="Trace Elements", relief=tk.RAISED)
            #
            self.container_gui["MINERALOGY"]["LABEL"].extend([lb_trace])
            #
            ## Radiobuttons
            rb_trace_01 = SE(
                parent=self.parent, row_id=15, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
                fg=self.color_bg).create_radiobutton(
                var_rb=self.var_rb_mineralogy_traces, value_rb=0, text="Without Traces",
                color_bg=self.color_menu, command=self.change_rb_trace)
            rb_trace_02 = SE(
                parent=self.parent, row_id=16, column_id=0, n_rows=1, n_columns=1, bg=self.color_menu,
                fg=self.color_bg).create_radiobutton(
                var_rb=self.var_rb_mineralogy_traces, value_rb=1, text="With Traces",
                color_bg=self.color_menu, command=self.change_rb_trace)
            #
            self.container_gui["MINERALOGY"]["RADIOBUTTON"].extend([rb_trace_01, rb_trace_02])
            #
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
            #
        ## CARBONATES
        elif var_opt in ["Calcite", "Dolomite", "Magnesite", "Siderite", "Rhodochrosite", "Aragonite", "Cerussite",
                         "Ankerite", "Azurite", "Malachite", "Ikaite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        ## HALIDES
        elif var_opt in ["Halite", "Fluorite", "Sylvite", "Carnallite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # SULFATES
        elif var_opt in ["Barite", "Celestite", "Anglesite", "Anhydrite", "Hanksite", "Gypsum", "Alunite", "Jarosite",
                         "Chalcanthite", "Kieserite", "Scheelite", "Hexahydrite", "Kainite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # TECTOSILICATES
        elif var_opt in ["Alkalifeldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # PHOSPHATES
        elif var_opt in ["Apatite-F", "Apatite-Cl", "Apatite-OH", "Apatite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # PHOSPHIDES
        elif var_opt in ["Allabogdanite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # PHYLLOSILICATES
        elif var_opt in ["Illite", "Kaolinite", "Montmorillonite", "Chamosite", "Clinochlore", "Pennantite", "Nimite",
                         "Chlorite", "Vermiculite", "Annite", "Phlogopite", "Eastonite", "Siderophyllite", "Biotite",
                         "Muscovite", "Glauconite", "Nontronite", "Saponite", "Talc", "Chrysotile", "Antigorite",
                         "Pyrophyllite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # NESOSILICATES
        elif var_opt in ["Zircon", "Titanite", "Thorite", "Andalusite", "Kyanite", "Sillimanite", "Topaz", "Staurolite",
                         "Fayalite", "Forsterite", "Tephroite", "Calcio-Olivine", "Liebenbergite", "Olivine", "Pyrope",
                         "Almandine", "Grossular", "Andradite", "Uvarovite", "Aluminium Garnet", "Calcium Garnet"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # SOROSILICATES
        elif var_opt in ["Epidote", "Zoisite", "Gehlenite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # INOSILICATES
        elif var_opt in ["Enstatite", "Ferrosilite", "Diopside", "Jadeite", "Aegirine", "Spodumene", "Wollastonite",
                         "Tremolite", "Actinolite", "Glaucophane", "Augite", "Riebeckite", "Arfvedsonite",
                         "Calcium Amphibole", "Sodium Amphibole", "Mg-Fe Pyroxene", "Calcium Pyroxene", "Donpeacorite",
                         "Orthopyroxene"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # CYCLOSILICATES
        elif var_opt in ["Beryl", "Benitoite", "Cordierite", "Sekaninaite", "Schorl", "Elbaite", "Liddicoatite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # ORGANICS
        elif var_opt in ["Organic Matter"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        ## ROCKS
        elif var_opt in ["Sandstone", "Shale", "Mudstone", "Conglomerate"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        #
        elif var_opt in ["Limestone", "Limestone (old)", "Dolomite Rock", "Marl"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Custom Carbonate Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        #
        # IGNEOUS ROCKS (Volcanic Rocks)
        elif var_opt in ["Rhyolite (Streckeisen)", "Trachyte (Streckeisen)", "Latite (Streckeisen)", 
                         "Andesite (Streckeisen)", "Basalt (Streckeisen)", "Dacite (Streckeisen)"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # IGNEOUS ROCKS (Plutonic Rocks)
        elif var_opt in ["Granite (Streckeisen)", "Granodiorite (Streckeisen)", "Tonalite (Streckeisen)",
                         "Gabbro (Streckeisen)", "Diorite (Streckeisen)", "Monzonite (Streckeisen)",
                         "Syenite (Streckeisen)", "Granitoid (Streckeisen)", "Quarzolite (Streckeisen)",
                         "Foid-bearing Syenite (Streckeisen)", "Foid-bearing Monzonite (Streckeisen)",
                         "Foid-bearing Monzodiorite (Streckeisen)", "Foid-bearing Monzogabbro (Streckeisen)",
                         "Foid Monzosyenite (Streckeisen)", "Foid Monzodiorite (Streckeisen)",
                         "Foid Monzogabbro (Streckeisen)", "Foidolite (Streckeisen)"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        # IGNEOUS ROCKS (Pyroclastic Rocks)
        elif var_opt in ["Pyroclastic Rock"]:
            Rocks(
                parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Felsic Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Intermediate Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Granite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Gabbro":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Diorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Granodiorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Syenite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Tonalite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Monzonite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Quartzolite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        elif var_opt == "Qz-rich Granitoid":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        ## Evaporite Rocks
        elif var_opt in ["Rock Salt", "Anhydrite (Rock)", "Potash"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        ## Ore Rocks
        elif var_opt in ["Kupferschiefer", "Compact Hematite Ore", "Banded Iron Formation"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        ## Metamorphic Rocks
        elif var_opt in ["Granulite", "Greenschist", "Greenschist (basaltic)", "Greenschist (ultramafic)",
                         "Greenschist (pelitic)", "Amphibolite (ortho)"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements, exp_data=self.exp_data, filename=self.filename)
        #
        # OXIDES
        elif var_opt == "Oxides":
            try:
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            #
            var_opt_0_1 = tk.StringVar()
            var_opt_0_2 = tk.StringVar()
            opt_list_0_1 = ["Quartz", "Magnetite", "Hematite", "Aluminium Spinels", "Iron Spinels", "Chromium Spinels",
                            "Corundum", "Ilmenite", "Rutile", "Pyrolusite", "Cassiterite", "Chromite",
                            "Magnesiochromite", "Zincochromite", "Cuprospinel", "Jacobsite", "Magnesioferrite",
                            "Trevorite", "Franklinite", "Ulvöspinel", "Uraninite", "Litharge", "Massicot", "Minium",
                            "Plattnerite", "Scrutinyite", "Zincite", "Columbite", "Tantalite", "Coltan",
                            "Hematite-Group", "Rutile-Group", "Periclase-Group", "Crocoite", "Wulfenite",
                            "Wulfenite-Group", "Goethite"]
            opt_list_0_1.sort()
            opt_oxide_comparison = ["Spinel Group", "Hematite Group", "Rutile Group", "Periclase Group",
                                    "Wulfenite Group"]
            opt_oxide_comparison.sort()
            #
            self.opt_oxide = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=self.var_opt_mineral, var_opt_set="Select Oxide Mineral", opt_list=opt_list_0_1, active_bg=self.color_accent_02,
                command=lambda var_opt=self.var_opt_mineral: self.select_opt(var_opt))
            self.opt_oxide_comp = SE(
                parent=self.parent, row_id=12, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_2, var_opt_set="Compare Oxide Minerals", opt_list=opt_oxide_comparison,
                active_bg=self.color_accent_02, command=lambda var_opt=var_opt_0_2: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].extend([self.opt_oxide, self.opt_oxide_comp])
            #
        # SULFIDES
        elif var_opt == "Sulfides":
            try:
                self.opt_oxide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            self.var_opt_sulfide = tk.StringVar()
            self.var_opt_sulfide_group = tk.StringVar()
            opt_list_0_2 = ["Pyrite", "Chalcopyrite", "Galena", "Acanthite", "Chalcocite", "Bornite", "Sphalerite",
                            "Pyrrhotite", "Millerite", "Pentlandite", "Covellite", "Cinnabar", "Realgar", "Orpiment",
                            "Stibnite", "Marcasite", "Molybdenite", "Fahlore", "Chalcopyrite-Group", "Gallite",
                            "Roquesite", "Lenaite", "Laforetite", "Vaesite", "Cattierite", "Pyrite-Group"]
            opt_list_0_2.sort()
            opt_list_sulfide_groups = ["Pyrite Group"]
            opt_list_sulfide_groups.sort()
            #
            self.opt_sulfide = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                  bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=self.var_opt_mineral, var_opt_set="Select Sulfide Mineral", opt_list=opt_list_0_2,
                command=lambda var_opt=self.var_opt_mineral: self.select_opt(var_opt))
            self.opt_sulfide_groups = SE(
                parent=self.parent, row_id=12, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                fg=self.color_fg_dark).create_option_menu(
                var_opt=self.var_opt_sulfide_group, var_opt_set="Select Sulfide Group", opt_list=opt_list_sulfide_groups,
                active_bg=self.color_accent_02)
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].extend([self.opt_sulfide, self.opt_sulfide_groups])
            #
        # CARBONATES
        elif var_opt == "Carbonates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_3 = tk.StringVar()
            opt_list_0_3 = ["Calcite", "Dolomite", "Magnesite", "Siderite", "Rhodochrosite", "Aragonite", "Cerussite",
                            "Ankerite", "Azurite", "Malachite", "Ikaite"]
            opt_list_0_3.sort()
            self.opt_carb = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                               bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_3, var_opt_set="Select Carbonate Mineral", opt_list=opt_list_0_3,
                command=lambda var_opt=var_opt_0_3: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_carb)
            #
        # HALOGENES
        elif var_opt == "Halides":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_4 = tk.StringVar()
            opt_list_0_4 = ["Halite", "Fluorite", "Sylvite", "Carnallite"]
            opt_list_0_4.sort()
            self.opt_halogene = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                   bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_4, var_opt_set="Select Halogene Mineral", opt_list=opt_list_0_4,
                command=lambda var_opt=var_opt_0_4: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_halogene)
            #
        # TECTOSILICATES
        elif var_opt == "Tectosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_5 = tk.StringVar()
            opt_list_0_5 = ["Alkalifeldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline"]
            opt_list_0_5.sort()
            self.opt_afs = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_5, var_opt_set="Select Tectosilicate Mineral", opt_list=opt_list_0_5,
                command=lambda var_opt=var_opt_0_5: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_afs)
            #
        elif var_opt == "Phyllosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_6 = tk.StringVar()
            opt_list_0_6 = ["Illite", "Kaolinite", "Montmorillonite", "Chamosite", "Clinochlore", "Pennantite",
                            "Nimite", "Chlorite", "Vermiculite", "Annite", "Phlogopite", "Eastonite", "Siderophyllite",
                            "Biotite", "Muscovite", "Glauconite", "Nontronite", "Saponite", "Talc", "Chrysotile",
                            "Antigorite", "Pyrophyllite"]
            opt_list_0_6.sort()
            self.opt_clays = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_6, var_opt_set="Select Phyllosilicate Mineral", opt_list=opt_list_0_6,
                command=lambda var_opt=var_opt_0_6: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_clays)
            #
        elif var_opt == "Sulfates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_7 = tk.StringVar()
            opt_list_0_7 = ["Barite", "Celestite", "Anglesite", "Anhydrite", "Hanksite", "Gypsum", "Alunite", "Jarosite",
                            "Chalcanthite", "Kieserite", "Scheelite", "Hexahydrite", "Kainite"]
            opt_list_0_7.sort()
            self.opt_sulfate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_7, var_opt_set="Select Sulfate Mineral", opt_list=opt_list_0_7,
                command=lambda var_opt=var_opt_0_7: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_sulfate)
            #
        # NESOSILICATES
        elif var_opt == "Nesosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_8 = tk.StringVar()
            opt_list_0_8 = ["Zircon", "Titanite", "Thorite", "Andalusite", "Kyanite", "Sillimanite", "Topaz", "Staurolite",
                            "Fayalite", "Forsterite", "Tephroite", "Calcio-Olivine", "Liebenbergite", "Olivine",
                            "Pyrope", "Almandine", "Grossular", "Andradite", "Uvarovite", "Aluminium Garnet",
                            "Calcium Garnet"]
            opt_list_0_8.sort()
            self.opt_nesosilicate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_8, var_opt_set="Select Nesosilicate Mineral", opt_list=opt_list_0_8,
                command=lambda var_opt=var_opt_0_8: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_nesosilicate)
            #
        # SOROSILICATES
        elif var_opt == "Sorosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_9 = tk.StringVar()
            opt_list_0_9 = ["Epidote", "Zoisite", "Gehlenite"]
            opt_list_0_9.sort()
            self.opt_sorosilicate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_9, var_opt_set="Select Sorosilicate Mineral", opt_list=opt_list_0_9,
                command=lambda var_opt=var_opt_0_9: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_sorosilicate)
            #
        # INOSILICATES
        elif var_opt == "Inosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_10 = tk.StringVar()
            opt_list_0_10 = ["Enstatite", "Ferrosilite", "Diopside", "Jadeite", "Aegirine", "Spodumene", "Wollastonite",
                             "Tremolite", "Actinolite", "Glaucophane", "Augite", "Riebeckite", "Arfvedsonite",
                             "Calcium Amphibole", "Sodium Amphibole", "Mg-Fe Pyroxene", "Calcium Pyroxene",
                             "Donpeacorite", "Orthopyroxene"]
            opt_list_0_10.sort()
            self.opt_inosilicate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_10, var_opt_set="Select Inosilicate Mineral", opt_list=opt_list_0_10,
                command=lambda var_opt=var_opt_0_10: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_inosilicate)
            #
        # PHOSPHATES
        elif var_opt == "Phosphates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_11 = tk.StringVar()
            opt_list_0_11 = ["Apatite", "Apatite-F", "Apatite-Cl", "Apatite-OH"]
            opt_list_0_11.sort()
            self.opt_phosphates = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_11, var_opt_set="Select Phosphate Mineral", opt_list=opt_list_0_11,
                command=lambda var_opt=var_opt_0_11: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_phosphates)
            #
        # PHOSPHIDES
        elif var_opt == "Phosphides":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_orgmat.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_12 = tk.StringVar()
            opt_list_0_12 = ["Allabogdanite"]
            opt_list_0_12.sort()
            self.opt_phosphides = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_12, var_opt_set="Select Phosphide Mineral", opt_list=opt_list_0_12,
                command=lambda var_opt=var_opt_0_12: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_phosphides)
            #
        # ORGANICS
        elif var_opt == "Organics":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_cyclo.grid_remove()
            except:
                pass
            var_opt_0_13 = tk.StringVar()
            opt_list_0_13 = ["Organic Matter"]
            opt_list_0_13.sort()
            self.opt_orgmat = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_13, var_opt_set="Select Organic Material", opt_list=opt_list_0_13,
                command=lambda var_opt=var_opt_0_13: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_orgmat)
            #
        # CYCLOSILICATES
        elif var_opt == "Cyclosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_sulfate.grid_remove()
                self.opt_nesosilicate.grid_remove()
                self.opt_sorosilicate.grid_remove()
                self.opt_inosilicate.grid_remove()
                self.opt_phosphates.grid_remove()
                self.opt_phosphides.grid_remove()
                self.opt_orgmat.grid_remove()
            except:
                pass
            var_opt_0_14 = tk.StringVar()
            opt_list_0_14 = ["Beryl", "Benitoite", "Cordierite", "Sekaninaite", "Schorl", "Elbaite", "Liddicoatite"]
            opt_list_0_14.sort()
            self.opt_cyclo = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                bg=self.color_accent_02,
                                fg=self.color_fg_dark).create_option_menu(var_opt=var_opt_0_14, var_opt_set="Select Cyclosilicate Mineral",
                                                                          opt_list=opt_list_0_14, command=lambda var_opt=var_opt_0_14: self.select_opt(var_opt))
            #
            self.container_gui["MINERALOGY"]["OPTION MENU"].append(self.opt_cyclo)
            #
        elif var_opt == "Siliciclastic Rocks":
            var_opt_1_1 = tk.StringVar()
            opt_list_1_1 = ["Sandstone", "Shale", "Conglomerate", "Mudstone"]
            opt_list_1_1.sort()
            self.opt_silic = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_1, var_opt_set="Select Rock", opt_list=opt_list_1_1,
                command=lambda var_opt=var_opt_1_1: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_silic)
            #
        elif var_opt == "Carbonate Rocks":
            var_opt_1_2 = tk.StringVar()
            opt_list_1_2 = ["Limestone", "Limestone (old)", "Dolomite Rock", "Marl"]
            opt_list_1_2.sort()
            self.opt_carb = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                               bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_2, var_opt_set="Select Rock", opt_list=opt_list_1_2,
                command=lambda var_opt=var_opt_1_2: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_carb)
            #
        elif var_opt == "Plutonic Rocks":
            var_opt_1_3 = tk.StringVar()
            opt_list_1_3 = [
                "Felsic Rock", "Intermediate Rock", "Granite", "Gabbro", "Diorite", "Granodiorite", "Monzonite",
                "Syenite", "Tonalite", "Quartzolite", "Qz-rich Granitoid", "Granite (Streckeisen)",
                "Granodiorite (Streckeisen)", "Tonalite (Streckeisen)", "Gabbro (Streckeisen)", "Diorite (Streckeisen)",
                "Monzonite (Streckeisen)", "Syenite (Streckeisen)", "Granitoid (Streckeisen)",
                "Quarzolite (Streckeisen)", "Foid-bearing Syenite (Streckeisen)", "Foid-bearing Monzonite (Streckeisen)",
                "Foid-bearing Monzodiorite (Streckeisen)", "Foid-bearing Monzogabbro (Streckeisen)",
                "Foid Monzosyenite (Streckeisen)", "Foid Monzodiorite (Streckeisen)",
                "Foid Monzogabbro (Streckeisen)", "Foidolite (Streckeisen)"]
            opt_list_1_3.sort()
            self.opt_ign = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_3, var_opt_set="Select Rock", opt_list=opt_list_1_3,
                command=lambda var_opt=var_opt_1_3: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_ign)
            #
        elif var_opt == "Volcanic Rocks":
            var_opt_1_7 = tk.StringVar()
            opt_list_1_7 = [
                "Rhyolite (Streckeisen)", "Trachyte (Streckeisen)", "Latite (Streckeisen)", "Andesite (Streckeisen)",
                "Basalt (Streckeisen)", "Dacite (Streckeisen)"]
            opt_list_1_7.sort()
            self.opt_ign = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_7, var_opt_set="Select Rock", opt_list=opt_list_1_7,
                command=lambda var_opt=var_opt_1_7: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_ign)
            #
        elif var_opt == "Pyroclastic Rocks":
            var_opt_1_8 = tk.StringVar()
            opt_list_1_8 = [
                "Pyroclastic Rock"]
            opt_list_1_8.sort()
            self.opt_ign = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_8, var_opt_set="Select Rock", opt_list=opt_list_1_8,
                command=lambda var_opt=var_opt_1_8: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_ign)
            #
        elif var_opt == "Evaporite Rocks":
            var_opt_1_4 = tk.StringVar()
            opt_list_1_4 = ["Rock Salt", "Anhydrite (Rock)", "Potash"]
            opt_list_1_4.sort()
            self.opt_evapr = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_4, var_opt_set="Select Rock", opt_list=opt_list_1_4,
                command=lambda var_opt=var_opt_1_4: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_evapr)
            #
        elif var_opt == "Ore Rocks":
            var_opt_1_5 = tk.StringVar()
            opt_list_1_5 = ["Kupferschiefer", "Compact Hematite Ore", "Banded Iron Formation"]
            opt_list_1_5.sort()
            self.opt_ore = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_5, var_opt_set="Select Rock", opt_list=opt_list_1_5,
                command=lambda var_opt=var_opt_1_5: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_ore)
            #
        elif var_opt == "Metamorphic Rocks":
            var_opt_1_6 = tk.StringVar()
            opt_list_1_6 = ["Granulite", "Greenschist", "Greenschist (basaltic)", "Greenschist (ultramafic)",
                            "Greenschist (pelitic)", "Amphibolite (ortho)"]
            opt_list_1_6.sort()
            self.opt_metamorph = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_6, var_opt_set="Select Rock", opt_list=opt_list_1_6,
                command=lambda var_opt=var_opt_1_6: self.select_opt(var_opt))
            #
            self.container_gui["PETROLOGY"]["OPTION MENU"].append(self.opt_metamorph)
            #
        elif var_opt == "Zechstein":
            thickness_complete = rd.randrange(900, 1500, 100)
            thickness_z5_random = int(rd.uniform(0.02, 0.06)*thickness_complete)
            thickness_z4_random = int(rd.uniform(0.18, 0.22)*thickness_complete)
            thickness_z3_random = int(rd.uniform(0.31, 0.35)*thickness_complete)
            thickness_z2_random = int(rd.uniform(0.30, 0.34)*thickness_complete)
            thickness_z1_random = int(thickness_complete - thickness_z5_random- thickness_z4_random
                                      - thickness_z3_random- thickness_z2_random)
            #
            # data_z5 = Zechstein(actual_thickness=0).create_zechstein_z5(top_z=0, thickness_z5=100)      # Ohre
            # data_z4 = Zechstein(actual_thickness=0).create_zechstein_z4(top_z=100, thickness_z4=200)    # Aller
            # data_z3 = Zechstein(actual_thickness=0).create_zechstein_z3(top_z=300, thickness_z3=300)    # Leine
            # data_z2 = Zechstein(actual_thickness=0).create_zechstein_z2(top_z=600, thickness_z2=300)    # Straßfurt
            # data_z1 = Zechstein(actual_thickness=300).create_zechstein_z1(top_z=900, thickness_z1=200)  # Werra
            data_z5 = Zechstein(actual_thickness=0).create_zechstein_z5(
                top_z=0,
                thickness_z5=thickness_z5_random)   # Ohre
            data_z4 = Zechstein(actual_thickness=0).create_zechstein_z4(
                top_z=thickness_z5_random,
                thickness_z4=thickness_z4_random)   # Aller
            data_z3 = Zechstein(actual_thickness=0).create_zechstein_z3(
                top_z=thickness_z5_random+thickness_z4_random,
                thickness_z3=thickness_z3_random)   # Leine
            data_z2 = Zechstein(actual_thickness=0).create_zechstein_z2(
                top_z=thickness_z5_random+thickness_z4_random+thickness_z3_random,
                thickness_z2=thickness_z2_random)   # Straßfurt
            data_z1 = Zechstein(actual_thickness=300).create_zechstein_z1(
                top_z=thickness_z5_random+thickness_z4_random+thickness_z3_random+thickness_z2_random,
                thickness_z1=thickness_z1_random)   # Werra
            #
            data_zechstein = data_z5 + data_z4 + data_z3 + data_z2 + data_z1
            Subsurface(
                parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_opt, lbl_w=self.lbl_w,
                entr_w=self.entr_w, gui_elements=self.gui_elements).create_real_world_sequences(
                name=var_opt, data_units=data_zechstein)
        #
        elif var_opt == "Muschelkalk":
            thickness_complete = rd.randrange(900, 1500, 100)
            thickness_muschelkalk_oberer_random = int(rd.uniform(0.3, 0.4) * thickness_complete)
            thickness_muschelkalk_mittlerer_random = int(rd.uniform(0.3, 0.4) * thickness_complete)
            thickness_muschelkalk_unterer_random = int(thickness_complete - thickness_muschelkalk_oberer_random
                                                       - thickness_muschelkalk_mittlerer_random)
            #
            data_muschelkalk_oberer = Muschelkalk(
                actual_thickness=0, resolution=10).create_muschelkalk_oberer(
                top_unit=0, thickness_unit=thickness_muschelkalk_oberer_random)  # Oberer Muschelkalk
            data_muschelkalk_mittlerer = Muschelkalk(
                actual_thickness=0, resolution=10).create_muschelkalk_mittlerer(
                top_unit=thickness_muschelkalk_oberer_random, thickness_unit=thickness_muschelkalk_mittlerer_random)  # Mittlerer Muschelkalk
            data_muschelkalk_unterer = Muschelkalk(
                actual_thickness=0).create_muschelkalk_unterer(
                top_unit=thickness_muschelkalk_oberer_random+thickness_muschelkalk_mittlerer_random,
                thickness_unit=thickness_muschelkalk_unterer_random)  # Unterer Muschelkalk
            #
            data_muschelkalk = data_muschelkalk_oberer + data_muschelkalk_mittlerer + data_muschelkalk_unterer
            Subsurface(
                parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_dark,
                color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_opt, lbl_w=self.lbl_w,
                entr_w=self.entr_w, gui_elements=self.gui_elements).create_real_world_sequences(
                name=var_opt, data_units=data_muschelkalk)
    #
    def change_rb_trace(self):
        # ## Cleaning
        # for sub_category in self.sub_categories_gui:
        #     if len(self.container_gui_specific["MINERALOGY"][sub_category]) > 0:
        #         for gui_item in self.container_gui_specific["MINERALOGY"][sub_category]:
        #             if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
        #                 gui_item.grid_remove()
        #             else:
        #                 gui_item.get_tk_widget().grid_remove()
        #
        if self.var_rb_mineralogy_traces.get() == 0:
            ## Cleaning
            for sub_category in self.sub_categories_gui:
                for gui_item in self.container_gui_specific["MINERALOGY"][sub_category]:
                    if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                        gui_item.grid_remove()
                    else:
                        gui_item.get_tk_widget().grid_remove()
                self.container_gui_specific["MINERALOGY"][sub_category].clear()
            #
        elif self.var_rb_mineralogy_traces.get() == 1:
            ## Cleaning
            for sub_category in self.sub_categories_gui:
                for gui_item in self.container_gui_specific["MINERALOGY"][sub_category]:
                    if sub_category not in ["CANVAS_HISTO", "CANVAS_SCATTER", "CANVAS_GEOCHEM"]:
                        gui_item.grid_remove()
                    else:
                        gui_item.get_tk_widget().grid_remove()
                self.container_gui_specific["MINERALOGY"][sub_category].clear()
            #
            if len(self.container_gui_specific["MINERALOGY"]["BUTTON"]) == 0:
                if self.var_opt_mineral_group.get() == "Oxides":
                    data_mineral = Oxides(mineral=self.var_opt_mineral.get(), data_type=True).get_data()
                elif self.var_opt_mineral_group.get() == "Sulfides":
                    data_mineral = Sulfides(mineral=self.var_opt_mineral.get(), data_type=True).get_data()
                btn_trace = SE(
                        parent=self.parent, row_id=15, column_id=1, n_rows=2, n_columns=1, bg=self.color_accent_02,
                        fg=self.color_fg_dark).create_button(
                        text="Select\n Trace Elements",
                        command=lambda traces_list=data_mineral["trace elements"]: self.select_trace_elements(traces_list))
                #
                self.container_gui_specific["MINERALOGY"]["BUTTON"].append(btn_trace)
                #
                lb_trace = SE(parent=self.parent, row_id=18, column_id=0, n_rows=1, n_columns=2, bg=self.color_accent_03,
                              fg=self.color_fg_dark).create_label(text="Test", relief=tk.RAISED)
                #
                self.container_gui_specific["MINERALOGY"]["LABEL"].extend([lb_trace])
    #
    def select_trace_elements(self, traces_list):
        var_opt_traces = tk.StringVar()
        traces_list.sort()
        opt_traces = SE(
            parent=self.parent, row_id=17, column_id=0, n_rows=1, n_columns=2, bg=self.color_accent_02,
            fg=self.color_fg_dark).create_option_menu(
            var_opt=var_opt_traces, var_opt_set="Select Trace Element", opt_list=traces_list)
        #
        self.container_gui_specific["MINERALOGY"]["OPTION MENU"].append(opt_traces)
    #
    def plot_histogram_comparison(self, dataset):
        labels = [["Density - kg/m$^3$", "Gamma Ray - API", "Photoelectricity - barns/e$^-$"],
                  ["vP - m/s", "vS - m/s", "vP/vS - 1"], ["K - GPa", "G - GPa", "Poisson - 1"]]
        labels_key = [["rho", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "nu"]]
        #
        sns.set_theme(style="darkgrid")
        self.fig_histo, ax_histo = plt.subplots(ncols=3, nrows=3, figsize=(11, 11), facecolor="#E9ECED")
        #
        for key, dataset_key in dataset.items():
            for i in range(3):
                for j in range(3):
                    sns.histplot(
                        data=dataset_key[labels_key[i][j]], ax=ax_histo[i][j], bins=10, color=dataset_key["color"],
                        element="step", label=key, alpha=0.25, linewidth=2)
                    #
                    ax_histo[i][j].set_yscale("log")
                    ax_histo[i][j].set_xlabel(labels[i][j], fontsize="small")
                    ax_histo[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize=9)
                    ax_histo[i][j].grid(True)
                    ax_histo[i][j].set_axisbelow(True)
        #
        self.fig_histo.suptitle("\n")
        handles, labels = ax_histo[2][2].get_legend_handles_labels()
        self.leg = self.fig_histo.legend(handles, labels, loc='upper center', ncol=len(self.data_minerals))
        self.leg.set_in_layout(False)
        self.fig_histo.tight_layout()
        #
        canvas_histo = FigureCanvasTkAgg(self.fig_histo, master=self.parent)
        canvas_histo.get_tk_widget().grid(row=0, column=9, rowspan=47, columnspan=9, sticky="nesw")
        #
        self.container_gui["MINERALOGY"]["CANVAS_HISTO"].append(canvas_histo)
    #
    def plot_scatter_comparison(self, dataset):
        labels = [["Molar Mass - kg/mol", "Gamma Ray - API", "Photoelectricity - barns/e$^-$"],
                  ["vP - m/s", "vS - m/s", "vP/vS - 1"], ["K - GPa", "G - GPa", "Poisson - 1"]]
        labels_key = [["M", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "nu"]]
        #
        sns.set_theme(style="darkgrid")
        self.fig_scatter, ax_scatter = plt.subplots(ncols=3, nrows=3, figsize=(11, 11), facecolor="#E9ECED")
        #
        for key, dataset_key in dataset.items():
            for i in range(3):
                for j in range(3):
                    sns.scatterplot(
                        x=dataset_key["rho"], y=dataset_key[labels_key[i][j]], ax=ax_scatter[i][j],
                        color=dataset_key["color"], label=key, legend=False)
                    #
                    ax_scatter[i][j].set_xlabel("Density - kg/m$^3$", fontsize=9)
                    ax_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=9)
                    ax_scatter[i][j].grid(True)
                    ax_scatter[i][j].set_axisbelow(True)
        #
        self.fig_scatter.suptitle("\n")
        handles, labels = ax_scatter[2][2].get_legend_handles_labels()
        self.leg = self.fig_scatter.legend(handles, labels, loc='upper center', ncol=len(self.data_minerals))
        self.leg.set_in_layout(False)
        self.fig_scatter.tight_layout()
        #
        canvas_histo = FigureCanvasTkAgg(self.fig_scatter, master=self.parent)
        canvas_histo.get_tk_widget().grid(row=0, column=9, rowspan=47, columnspan=9, sticky="nesw")
        #
        self.container_gui["MINERALOGY"]["CANVAS_SCATTER"].append(canvas_histo)
    #
    def plot_geochem_comparison(self):
        labels = [["Density - kg/m$^3$", "Gamma Ray - API", "Photoelectricity - barns/e$^-$"],
                  ["vP - m/s", "vS - m/s", "vP/vS - 1"], ["K - GPa", "G - GPa", "Poisson - 1"]]
        labels_key = [["rho", "GR", "PE"], ["vP", "vS", "vP/vS"], ["K", "G", "nu"]]
        #
        sns.set_theme(style="darkgrid")
        fig_geochem, ax_scatter = plt.subplots(ncols=3, nrows=3, figsize=(11, 11), facecolor="#E9ECED")
        #
        for key, dataset_key in self.data_minerals.items():
            for i in range(3):
                for j in range(3):
                    scatter_set = sns.scatterplot(
                        x=dataset_key["chemistry"][self.container_var["MINERALOGY"]["OPTION MENU"]["GEOCHEM"].get()],
                        y=dataset_key[labels_key[i][j]], ax=ax_scatter[i][j],
                        color=dataset_key["color"], label=key, legend=False)
                    #
                    ax_scatter[i][j].set_xlabel(
                        "w("+str(self.container_var["MINERALOGY"]["OPTION MENU"]["GEOCHEM"].get())+")", fontsize=9)
                    ax_scatter[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize=9)
                    ax_scatter[i][j].grid(True)
                    ax_scatter[i][j].set_axisbelow(True)
        #
        fig_geochem.suptitle("\n")
        handles, labels = ax_scatter[2][2].get_legend_handles_labels()
        self.leg = fig_geochem.legend(handles, labels, loc='upper center', ncol=len(self.data_minerals))
        self.leg.set_in_layout(False)
        fig_geochem.tight_layout()
        #
        canvas_geochem = FigureCanvasTkAgg(fig_geochem, master=self.parent)
        canvas_geochem.get_tk_widget().grid(row=0, column=9, rowspan=47, columnspan=9, sticky="nesw")
        #
        try:
            for canvas_item in self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"]:
                canvas_item.get_tk_widget().grid_remove()
        except:
            pass
        self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"].clear()
        self.container_gui["MINERALOGY"]["CANVAS_GEOCHEM"].append(canvas_geochem)
    #
    def fill_table_comparison(self, dataset):
        labels = ["M", "V", "rho", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE"]
        labels_key = ["Molar Mass", "Molar Volume", "Density", "vP", "vS", "vP/vS", "K", "G", "E", "Poisson", "GR",
                      "PE"]
        #
        for index, label in enumerate(labels):
            self.container_var["MINERALOGY"]["ENTRY"]["MIN"][labels_key[index]].set(
                round(np.min(dataset[label]), 3))
            self.container_var["MINERALOGY"]["ENTRY"]["MAX"][labels_key[index]].set(
                round(np.max(dataset[label]), 3))
            self.container_var["MINERALOGY"]["ENTRY"]["MEAN"][labels_key[index]].set(
                round(np.mean(dataset[label]), 3))
            self.container_var["MINERALOGY"]["ENTRY"]["STD"][labels_key[index]].set(
                round(np.std(dataset[label], ddof=1), 3))
    #
    def change_radiobutton_mode(self, var_rb_mode):
        if var_rb_mode.get() == 0:
            print(var_rb_mode.get())
        elif var_rb_mode.get() == 1:
            print(var_rb_mode.get())
        elif var_rb_mode.get() == 2:
            print(var_rb_mode.get())
    #
    def export_data(self, var_exp):
        header = list(var_exp[0].keys())
        header_bulk = list(var_exp[0].keys())
        elements = list(var_exp[0]["chemistry"].keys())
        n_elements = len(elements)
        header.remove("chemistry")
        header_bulk.remove("chemistry")
        if header[0] == "rock":
            minerals = list(var_exp[0]["mineralogy"].keys())
            n_minerals = len(minerals)
            header.remove("mineralogy")
            header_bulk.remove("mineralogy")
            header.extend(minerals)
        else:
            try:
                header.remove("state")
            except:
                pass
            header.remove("rho_e")
            header.remove("p")
            header.remove("U")
            header_bulk.remove("rho_e")
            header_bulk.remove("p")
            header_bulk.remove("U")
        header.extend(elements)
        #
        keys = re.split("\s+", self.filename[0])
        filename = "../outputs/Report_"
        for item in keys:
            filename = filename.__add__(str(item)+"_")
        filename = filename.__add__(str(self.filename[1]))
        filename = filename.__add__(".csv")
        #
        with open(filename, "w") as export_file:
            writer = csv.DictWriter(export_file, fieldnames=header)
            writer.writeheader()
            #
            if header[0] == "mineral":
                for item in var_exp:
                    row = ""
                    for bulk_key in header_bulk:
                        row = row.__add__(str(item[bulk_key])+";")
                    for element in elements:
                        row = row.__add__(str(item["chemistry"][element])+";")
                    #
                    export_file.write(row+"\n")
            else:
                for item in var_exp:
                    row = ""
                    for bulk_key in header_bulk:
                        row = row.__add__(str(item[bulk_key])+";")
                    for mineral in minerals:
                        row = row.__add__(str(item["mineralogy"][mineral])+";")
                    for element in elements:
                        row = row.__add__(str(item["chemistry"][element])+";")
                    #
                    export_file.write(row+"\n")
        export_file.close()
#
class Minerals:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, mineral, lbl_w, entr_w, gui_elements, exp_data, filename):
        #
        try:
            for lbl in lbl_w["physics"]:
                lbl.grid_forget()
            for entr in entr_w["physics"]:
                entr.grid_forget()
            lbl_w["physics"].clear()
            entr_w["physics"].clear()
        except:
            pass
        #
        try:
            for lbl in lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in entr_w["chemistry"]:
                entr.grid_forget()
            lbl_w["chemistry"].clear()
            entr_w["chemistry"].clear()
        except:
            pass
        #
        try:
            for gui_elmnt in gui_elements:
                gui_elmnt.grid_forget()
            gui_elmnt.clear()
        except:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        try:
            exp_data.clear()
            filename.clear()
        except:
            pass
        #
        self.parent_mineral = parent
        self.color_bg = color_bg
        self.color_fg = color_fg
        self.color_acc_01 = color_acc[0]
        self.color_acc_02 = color_acc[1]
        self.mineral = mineral
        self.var_rb = tk.IntVar()
        self.var_rb.set(0)
        self.var_rb_trace = tk.IntVar()
        self.var_rb_trace.set(0)
        self.var_entr = tk.IntVar()
        var_entr_start = 100
        self.lbl_w = lbl_w
        self.entr_w = entr_w
        self.gui_elements = gui_elements
        self.exp_data = exp_data
        self.filename = filename
        self.filename.extend([self.mineral, var_entr_start])
        self.container_gui = {}
        categories_gui = ["Label", "Button", "Entry", "Option Menu", "Checkbox", "Radiobutton"]
        categories_type = ["General", "Specific"]
        for category_type in categories_type:
            self.container_gui[category_type] = {}
            for category_gui in categories_gui:
                self.container_gui[category_type][category_gui] = []
        #
        ## Labels
        #
        lbl_stat = SE(parent=self.parent_mineral, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
           fg="black").create_label(text="Statistics - "+str(self.mineral), relief=tk.RAISED)
        lbl_plt = SE(parent=self.parent_mineral, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
           fg="black").create_label(text="Plots - "+str(self.mineral), relief=tk.RAISED)
        #
        self.gui_elements.extend([lbl_stat, lbl_plt])
        #
        lb_param = SE(
            parent=self.parent_mineral, row_id=2, column_id=3, n_rows=2, bg=self.color_bg,
            fg="black").create_label(text="Parameter", relief=tk.RAISED)
        lb_min = SE(
            parent=self.parent_mineral, row_id=2, column_id=4, n_rows=2, bg=self.color_bg,
            fg="black").create_label(text="Minimum", relief=tk.RAISED)
        lb_max = SE(
            parent=self.parent_mineral, row_id=2, column_id=5, n_rows=2, bg=self.color_bg,
            fg="black").create_label(text="Maximum", relief=tk.RAISED)
        lb_mean = SE(
            parent=self.parent_mineral, row_id=2, column_id=6, n_rows=2, bg=self.color_bg,
            fg="black").create_label(text="Mean", relief=tk.RAISED)
        lb_std = SE(
            parent=self.parent_mineral, row_id=2, column_id=7, n_rows=2, bg=self.color_bg,
            fg="black").create_label(text="Standard\n Deviation", relief=tk.RAISED)
        #
        self.gui_elements.extend([lb_param, lb_min, lb_max, lb_mean, lb_std])
        #
        lbl_01 = SE(parent=self.parent_mineral, row_id=4, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Molar mass\n (g/mol)", relief=tk.RAISED)
        lbl_02 = SE(parent=self.parent_mineral, row_id=6, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Density\n (kg/m3)", relief=tk.RAISED)
        lbl_03 = SE(parent=self.parent_mineral, row_id=8, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="P-wave velocity\n (m/s)", relief=tk.RAISED)
        lbl_04 = SE(parent=self.parent_mineral, row_id=10, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="S-wave velocity\n (m/s)", relief=tk.RAISED)
        lbl_05 = SE(parent=self.parent_mineral, row_id=12, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Velocity ratio\n (1)", relief=tk.RAISED)
        lbl_06 = SE(parent=self.parent_mineral, row_id=14, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Bulk modulus\n (GPa)", relief=tk.RAISED)
        lbl_07 = SE(parent=self.parent_mineral, row_id=16, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Shear modulus\n (GPa)", relief=tk.RAISED)
        lbl_08 = SE(parent=self.parent_mineral, row_id=18, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Poisson's ratio\n (1)", relief=tk.RAISED)
        lbl_09 = SE(parent=self.parent_mineral, row_id=20, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Gamma ray\n (API)", relief=tk.RAISED)
        lbl_10 = SE(parent=self.parent_mineral, row_id=22, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Photoelectricity\n (barns/electron)", relief=tk.RAISED)
        #
        lbl_11 = SE(parent=self.parent_mineral, row_id=29, column_id=0, bg=self.color_acc_01,
           fg="black").create_label(text="Number of samples", relief=tk.RAISED)
        #
        lbl_12 = SE(parent=self.parent_mineral, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                    fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        #
        self.lbl_w["physics"].extend([lbl_01, lbl_02, lbl_03, lbl_04, lbl_05, lbl_06, lbl_07, lbl_08, lbl_09, lbl_10,
                                      lbl_11, lbl_12])
        #
        ## Entry
        entr_01 = SE(parent=self.parent_mineral, row_id=29, column_id=1, bg=self.color_acc_02,
                     fg=color_fg).create_entry(var_entr=self.var_entr, var_entr_set=var_entr_start,
                                               command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        self.entr_w["physics"].append(entr_01)
        #
        ## Radiobuttons
        rb_01 = SE(parent=self.parent_mineral, row_id=30, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, value_rb=0, text="Histogram",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_mineral, row_id=31, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, value_rb=1, text="Scatter plot",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        #
        rb_03 = SE(parent=self.parent_mineral, row_id=32, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_trace, value_rb=2,
                                          text="Without Trace Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_trace: self.change_radiobutton(var_rb))
        rb_04 = SE(parent=self.parent_mineral, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_trace, value_rb=3,
                                          text="With Trace Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_trace: self.change_radiobutton(var_rb))
        #
        self.gui_elements.extend([rb_01, rb_02, rb_03, rb_04])
        #
        self.var_dict = True
        data_all = []
        for i in range(var_entr_start):
            # Oxides
            if self.mineral == "Quartz":
                data = Oxides(impurity="pure", data_type=True).create_quartz()
            elif self.mineral == "Magnetite":
                data = Oxides(impurity="pure", data_type=True).create_magnetite()
            elif self.mineral == "Goethite":
                data = Oxides(impurity="pure", data_type=True).create_goethite()
            elif self.mineral == "Crocoite":
                data = Oxides(impurity="pure", data_type=True).create_crocoite()
            elif self.mineral == "Wulfenite":
                data = Oxides(impurity="pure", data_type=True).create_wulfenite()
            elif self.mineral == "Hematite-Group":
                data = Oxides(impurity="pure", data_type=True).create_hematite_group()
            elif self.mineral == "Rutile-Group":
                data = RutileGroup().create_rutile_group()
            elif self.mineral == "Periclase-Group":
                data = PericlaseGroup().create_periclase_group()
            elif self.mineral == "Wulfenite-Group":
                data = WulfeniteGroup().create_wulfenite_group()
            elif self.mineral == "Cassiterite":
                data = Oxides(impurity="pure", data_type=True).create_cassiterite()
            elif self.mineral == "Pyrolusite":
                data = Oxides(impurity="pure", data_type=True).create_pyrolusite()
            elif self.mineral == "Chromite":
                data = Oxides(impurity="pure", data_type=True).create_chromite()
            elif self.mineral == "Magnesiochromite":
                data = Oxides(impurity="pure", data_type=True).create_magnesiochromite()
            elif self.mineral == "Zincochromite":
                data = Oxides(impurity="pure", data_type=True).create_zincochromite()
            elif self.mineral == "Rutile":
                data = Oxides(impurity="pure", data_type=True).create_rutile()
            elif self.mineral == "Hematite":
                data = Oxides(impurity="pure", data_type=True).create_hematite()
            elif self.mineral == "Ilmenite":
                data = Oxides(impurity="pure", data_type=True).create_ilmenite()
            elif self.mineral == "Litharge":
                data = Oxides(impurity="pure", data_type=True).create_litharge()
            elif self.mineral == "Massicot":
                data = Oxides(impurity="pure", data_type=True).create_massicot()
            elif self.mineral == "Minium":
                data = Oxides(impurity="pure", data_type=True).create_minium()
            elif self.mineral == "Plattnerite":
                data = Oxides(impurity="pure", data_type=True).create_plattnerite()
            elif self.mineral == "Scrutinyite":
                data = Oxides(impurity="pure", data_type=True).create_scrutinyite()
            elif self.mineral == "Zincite":
                data = Oxides(impurity="pure", data_type=True).create_zincite()
            elif self.mineral == "Cuprospinel":
                data = Oxides(impurity="pure", data_type=True).create_cuprospinel()
            elif self.mineral == "Jacobsite":
                data = Oxides(impurity="pure", data_type=True).create_jacobsite()
            elif self.mineral == "Uraninite":
                data = Oxides(impurity="pure", data_type=True).create_uraninite()
            elif self.mineral == "Magnesioferrite":
                data = Oxides(impurity="pure", data_type=True).create_magnesioferrite()
            elif self.mineral == "Aluminium Spinels":
                data = Oxides(impurity="pure", data_type=True).create_aluminium_spinel()
            elif self.mineral == "Chromium Spinels":
                data = Oxides(impurity="pure", data_type=True).create_chromium_spinel()
            elif self.mineral == "Iron Spinels":
                data = Oxides(impurity="pure", data_type=True).create_iron_spinel()
            elif self.mineral == "Columbite":
                data = Oxides(impurity="pure", data_type=True).create_columbite()
            elif self.mineral == "Tantalite":
                data = Oxides(impurity="pure", data_type=True).create_tantalite()
            elif self.mineral == "Coltan":
                data = Oxides(impurity="pure", data_type=True).create_coltan()
            elif self.mineral == "Trevorite":
                data = Oxides(impurity="pure", data_type=True).create_trevorite()
            elif self.mineral == "Franklinite":
                data = Oxides(impurity="pure", data_type=True).create_franklinite()
            elif self.mineral == "Ulvöspinel":
                data = Oxides(impurity="pure", data_type=True).create_ulvoespinel()
            elif self.mineral == "Pyrite":
                data = Sulfides(impurity="pure", data_type=True).create_pyrite()
            elif self.mineral == "Vaesite":
                data = Sulfides(impurity="pure", data_type=True).create_vaesite()
            elif self.mineral == "Cattierite":
                data = Sulfides(impurity="pure", data_type=True).create_cattierite()
            elif self.mineral == "Pyrite-Group":
                data = Sulfides(impurity="pure", data_type=True).create_pyrite_group()
            elif self.mineral == "Chalcopyrite":
                data = Sulfides(impurity="pure", data_type=True).create_chalcopyrite()
            elif self.mineral == "Galena":
                data = Sulfides(impurity="pure", data_type=True).create_galena()
            elif self.mineral == "Acanthite":
                data = Sulfides(impurity="pure", data_type=True).create_acanthite()
            elif self.mineral == "Chalcocite":
                data = Sulfides(impurity="pure", data_type=True).create_chalcocite()
            elif self.mineral == "Bornite":
                data = Sulfides(impurity="pure", data_type=True).create_bornite()
            elif self.mineral == "Sphalerite":
                data = Sulfides(impurity="pure", data_type=True).create_sphalerite()
            elif self.mineral == "Pyrrhotite":
                # data = Sulfides(impurity="pure", data_type=True).create_pyrrhotite()
                data = Sulfides(impurity="pure", data_type=True).create_pyrrhotite_group()
            elif self.mineral == "Millerite":
                data = Sulfides(impurity="pure", data_type=True).create_millerite()
            elif self.mineral == "Pentlandite":
                data = Sulfides(impurity="pure", data_type=True).create_pentlandite()
            elif self.mineral == "Covellite":
                data = Sulfides(impurity="pure", data_type=True).create_covellite()
            elif self.mineral == "Cinnabar":
                data = Sulfides(impurity="pure", data_type=True).create_cinnabar()
            elif self.mineral == "Stibnite":
                data = Sulfides(impurity="pure", data_type=True).create_stibnite()
            elif self.mineral == "Molybdenite":
                data = Sulfides(impurity="pure", data_type=True).create_molybdenite()
            elif self.mineral == "Realgar":
                data = Sulfides(impurity="pure", data_type=True).create_realgar()
            elif self.mineral == "Orpiment":
                data = Sulfides(impurity="pure", data_type=True).create_orpiment()
            elif self.mineral == "Marcasite":
                data = Sulfides(impurity="pure", data_type=True).create_marcasite()
            elif self.mineral == "Fahlore":
                data = Sulfides(impurity="pure", data_type=True).create_fahlore()
            elif self.mineral == "Chalcopyrite-Group":
                data = Sulfides(impurity="pure", data_type=True).create_chalcopyrite_group()
            elif self.mineral == "Gallite":
                data = Sulfides(impurity="pure", data_type=True).create_gallite()
            elif self.mineral == "Roquesite":
                data = Sulfides(impurity="pure", data_type=True).create_roquesite()
            elif self.mineral == "Lenaite":
                data = Sulfides(impurity="pure", data_type=True).create_lenaite()
            elif self.mineral == "Laforetite":
                data = Sulfides(impurity="pure", data_type=True).create_laforetite()
            # Carbonates
            elif self.mineral == "Calcite":
                data = Carbonates(impurity="pure", data_type=True).create_calcite()
            elif self.mineral == "Corundum":
                data = Oxides(impurity="pure", data_type=True).create_corundum()
            elif self.mineral == "Dolomite":
                data = Carbonates(impurity="pure", data_type=True).create_dolomite()
            elif self.mineral == "Magnesite":
                data = Carbonates(impurity="pure", data_type=True).create_magnesite()
            elif self.mineral == "Siderite":
                data = Carbonates(impurity="pure", data_type=True).create_siderite()
            elif self.mineral == "Rhodochrosite":
                data = Carbonates(impurity="pure", data_type=True).create_rhodochrosite()
            elif self.mineral == "Aragonite":
                data = Carbonates(impurity="pure", data_type=True).create_aragonite()
            elif self.mineral == "Cerussite":
                data = Carbonates(impurity="pure", data_type=True).create_cerussite()
            elif self.mineral == "Ankerite":
                data = Carbonates(impurity="pure", data_type=True).create_ankerite()
            elif self.mineral == "Azurite":
                data = Carbonates(impurity="pure", data_type=True).create_azurite()
            elif self.mineral == "Malachite":
                data = Carbonates(impurity="pure", data_type=True).create_malachite()
            elif self.mineral == "Ikaite":
                data = Carbonates(impurity="pure", data_type=True).create_ikaite()
            # Halides
            elif self.mineral == "Halite":
                data = Halides(impurity="pure", dict=True).create_halite()
            elif self.mineral == "Fluorite":
                data = Halides(impurity="pure", dict=True).create_fluorite()
            elif self.mineral == "Sylvite":
                data = Halides(impurity="pure", dict=True).create_sylvite()
            elif self.mineral == "Carnallite":
                data = Halides(impurity="pure", dict=True).create_carnallite()
            # Phyllosilicates
            elif self.mineral == "Illite":
                data = Phyllosilicates(data_type=True).create_illite()
            elif self.mineral == "Kaolinite":
                data = Phyllosilicates(data_type=True).create_kaolinite()
            elif self.mineral == "Montmorillonite":
                data = Phyllosilicates(data_type=True).create_montmorillonite()
            elif self.mineral == "Chamosite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_chamosite()
            elif self.mineral == "Clinochlore":
                data = Phyllosilicates(impurity="pure", data_type=True).create_clinochlore()
            elif self.mineral == "Pennantite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_pennantite()
            elif self.mineral == "Nimite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_nimite()
            elif self.mineral == "Chlorite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            elif self.mineral == "Vermiculite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_vermiculite()
            elif self.mineral == "Annite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_annite()
            elif self.mineral == "Phlogopite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_phlogopite()
            elif self.mineral == "Eastonite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_eastonite()
            elif self.mineral == "Siderophyllite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_siderophyllite()
            elif self.mineral == "Biotite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            elif self.mineral == "Muscovite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
            elif self.mineral == "Glauconite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_glauconite()
            elif self.mineral == "Nontronite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_nontronite()
            elif self.mineral == "Saponite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_saponite()
            elif self.mineral == "Talc":
                data = Phyllosilicates(impurity="pure", data_type=True).create_talc()
            elif self.mineral == "Chrysotile":
                data = Phyllosilicates(impurity="pure", data_type=True).create_chrysotile()
            elif self.mineral == "Antigorite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_antigorite()
            elif self.mineral == "Pyrophyllite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_pyrophyllite()
            # Tectosilicates
            elif self.mineral == "Alkalifeldspar":
                data = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            elif self.mineral == "Scapolite":
                data = Tectosilicates(impurity="pure", data_type=True).create_scapolite()
            elif self.mineral == "Danburite":
                data = Tectosilicates(impurity="pure", data_type=True).create_danburite()
            elif self.mineral == "Nepheline":
                data = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
            # Sulfates
            elif self.mineral == "Barite":
                data = Sulfates(impurity="pure", data_type=True).create_barite()
            elif self.mineral == "Anhydrite":
                data = Sulfates(impurity="pure", data_type=True).create_anhydrite()
            elif self.mineral == "Gypsum":
                data = Sulfates(impurity="pure", data_type=True).create_gypsum()
            elif self.mineral == "Celestite":
                data = Sulfates(impurity="pure", data_type=True).create_celestine()
            elif self.mineral == "Anglesite":
                data = Sulfates(impurity="pure", data_type=True).create_anglesite()
            elif self.mineral == "Hanksite":
                data = Sulfates(impurity="pure", data_type=True).create_hanksite()
            elif self.mineral == "Jarosite":
                data = Sulfates(impurity="pure", data_type=True).create_jarosite()
            elif self.mineral == "Alunite":
                data = Sulfates(impurity="pure", data_type=True).create_alunite()
            elif self.mineral == "Chalcanthite":
                data = Sulfates(impurity="pure", data_type=True).create_chalcanthite()
            elif self.mineral == "Kieserite":
                data = Sulfates(impurity="pure", data_type=True).create_kieserite()
            elif self.mineral == "Scheelite":
                data = Sulfates(impurity="pure", data_type=True).create_scheelite()
            elif self.mineral == "Hexahydrite":
                data = Sulfates(impurity="pure", data_type=True).create_hexahydrite()
            elif self.mineral == "Kainite":
                data = Sulfates(impurity="pure", data_type=True).create_kainite()
            # Nesosilicates
            elif self.mineral == "Zircon":
                data = Nesosilicates(data_type=True).create_zircon()
            elif self.mineral == "Thorite":
                data = Nesosilicates(data_type=True).create_thorite()
            elif self.mineral == "Titanite":
                data = Nesosilicates(data_type=True).create_titanite()
            elif self.mineral == "Andalusite":
                data = Nesosilicates(data_type=True).create_andalusite()
            elif self.mineral == "Kyanite":
                data = Nesosilicates(data_type=True).create_kyanite()
            elif self.mineral == "Sillimanite":
                data = Nesosilicates(data_type=True).create_sillimanite()
            elif self.mineral == "Topaz":
                data = Nesosilicates(data_type=True).create_topaz()
            elif self.mineral == "Staurolite":
                data = Nesosilicates(data_type=True).create_staurolite()
            elif self.mineral == "Forsterite":
                data = Nesosilicates(data_type=True).create_forsterite()
            elif self.mineral == "Fayalite":
                data = Nesosilicates(data_type=True).create_fayalite()
            elif self.mineral == "Tephroite":
                data = Nesosilicates(data_type=True).create_tephroite()
            elif self.mineral == "Calcio-Olivine":
                data = Nesosilicates(data_type=True).create_calcio_olivine()
            elif self.mineral == "Liebenbergite":
                data = Nesosilicates(data_type=True).create_liebenbergite()
            elif self.mineral == "Olivine":
                data = Nesosilicates(data_type=True).create_olivine()
            elif self.mineral == "Pyrope":
                data = Nesosilicates(data_type=True).create_pyrope()
            elif self.mineral == "Almandine":
                data = Nesosilicates(data_type=True).create_almandine()
            elif self.mineral == "Grossular":
                data = Nesosilicates(data_type=True).create_grossular()
            elif self.mineral == "Andradite":
                data = Nesosilicates(data_type=True).create_andradite()
            elif self.mineral == "Uvarovite":
                data = Nesosilicates(data_type=True).create_uvarovite()
            elif self.mineral == "Aluminium Garnet":
                data = Nesosilicates(data_type=True).create_aluminium_garnet()
            elif self.mineral == "Calcium Garnet":
                data = Nesosilicates(data_type=True).create_calcium_garnet()
            # Sorosilicates
            elif self.mineral == "Epidote":
                data = Sorosilicates(data_type=True).create_epidote()
            elif self.mineral == "Zoisite":
                data = Sorosilicates(data_type=True).create_zoisite()
            elif self.mineral == "Gehlenite":
                data = Sorosilicates(data_type=True).create_gehlenite()
            # Inosilicates
            elif self.mineral == "Enstatite":
                data = Inosilicates(data_type=True).create_enstatite()
            elif self.mineral == "Ferrosilite":
                data = Inosilicates(data_type=True).create_ferrosilite()
            elif self.mineral == "Donpeacorite":
                data = Inosilicates(data_type=True).create_donpeacorite()
            elif self.mineral == "Orthopyroxene":
                data = Inosilicates(data_type=True).create_orthopyroxene()
            elif self.mineral == "Diopside":
                data = Inosilicates(data_type=True).create_diopside()
            elif self.mineral == "Jadeite":
                data = Inosilicates(data_type=True).create_jadeite()
            elif self.mineral == "Aegirine":
                data = Inosilicates(data_type=True).create_aegirine()
            elif self.mineral == "Spodumene":
                data = Inosilicates(data_type=True).create_spodumene()
            elif self.mineral == "Wollastonite":
                data = Inosilicates(data_type=True).create_wollastonite()
            elif self.mineral == "Tremolite":
                data = Inosilicates(data_type=True).create_tremolite()
            elif self.mineral == "Actinolite":
                data = Inosilicates(data_type=True).create_actinolite()
            elif self.mineral == "Glaucophane":
                data = Inosilicates(data_type=True).create_glaucophane()
            elif self.mineral == "Augite":
                data = Inosilicates(data_type=True).create_augite()
            elif self.mineral == "Riebeckite":
                data = Inosilicates(data_type=True).create_riebeckite()
            elif self.mineral == "Arfvedsonite":
                data = Inosilicates(data_type=True).create_arfvedsonite()
            elif self.mineral == "Calcium Amphibole":
                data = Inosilicates(data_type=True).create_calcium_amphibole()
            elif self.mineral == "Sodium Amphibole":
                data = Inosilicates(data_type=True).create_sodium_amphibole()
            elif self.mineral == "Mg-Fe Pyroxene":
                data = Inosilicates(data_type=True).create_mg_fe_pyroxene()
            elif self.mineral == "Calcium Pyroxene":
                data = Inosilicates(data_type=True).create_calium_pyroxene()
            # Cyclosilicates
            elif self.mineral == "Beryl":
                data = Cyclosilicates(data_type=True).create_beryl()
            elif self.mineral == "Benitoite":
                data = Cyclosilicates(data_type=True).create_benitoite()
            elif self.mineral == "Cordierite":
                data = Cyclosilicates(data_type=True).create_cordierite()
            elif self.mineral == "Sekaninaite":
                data = Cyclosilicates(data_type=True).create_sekaninaite()
            elif self.mineral == "Schorl":
                data = Cyclosilicates(data_type=True).create_schorl()
            elif self.mineral == "Elbaite":
                data = Cyclosilicates(data_type=True).create_elbaite()
            elif self.mineral == "Liddicoatite":
                data = Cyclosilicates(data_type=True).create_liddicoatite()
            # Phosphates
            elif self.mineral == "Apatite-F":
                data = Phosphates(data_type=True).create_aptite_f()
            elif self.mineral == "Apatite-Cl":
                data = Phosphates(data_type=True).create_aptite_cl()
            elif self.mineral == "Apatite-OH":
                data = Phosphates(data_type=True).create_aptite_oh()
            elif self.mineral == "Apatite":
                data = Phosphates(data_type=True).create_aptite()
            # Phosphides
            elif self.mineral == "Allabogdanite":
                data = Phospides(data_type=True).create_allabogdanite()
            # Organics
            elif self.mineral == "Organic Matter":
                data = Organics(data_type=True).create_organics_matter()
            #
            self.color_mineral = "#7C9097"
            #
            data_all.append(data)
        self.exp_data.extend(data_all)
        #
        if self.var_dict == False:
            elements = np.array(data_all[0][-1])[:, 0]
            self.element_list = []
            for element in elements:
                self.element_list.append(DP(dataset=data_all).extract_element_amounts(type="mineral", element=element)*100)
            #
            self.rho_b = DP(dataset=data_all).extract_densities(type="mineral", keyword="bulk")
            self.molar_mass = DP(dataset=data_all).extract_molar_mass()
            self.volume = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="V")
            self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="bulk")
            self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="shear")
            self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="poisson")
            self.vP = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vP")
            self.vS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vS")
            self.vPvS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vPvS")
            self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="mineral")
            self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="mineral")
            #
            if self.mineral == "Quartz":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Si")*100
            elif self.mineral == "Alkalifeldspar":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="K")*100
            elif self.mineral == "Plagioclase":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Ca")*100
        else:
            self.molar_mass = DP(dataset=data_all).extract_data(keyword="M")
            self.volume = DP(dataset=data_all).extract_data(keyword="V")
            self.rho_b = DP(dataset=data_all).extract_data(keyword="rho")
            self.vP = DP(dataset=data_all).extract_data(keyword="vP")
            self.vS = DP(dataset=data_all).extract_data(keyword="vS")
            self.vPvS = DP(dataset=data_all).extract_data(keyword="vP/vS")
            self.bulk_mod = DP(dataset=data_all).extract_data(keyword="K")
            self.shear_mod = DP(dataset=data_all).extract_data(keyword="G")
            self.youngs_mod = DP(dataset=data_all).extract_data(keyword="E")
            self.poisson = DP(dataset=data_all).extract_data(keyword="nu")
            self.gamma_ray = DP(dataset=data_all).extract_data(keyword="GR")
            self.photoelectricity = DP(dataset=data_all).extract_data(keyword="PE")
            self.chemistry = DP(dataset=data_all).extract_data(keyword="chemistry")
            elements = np.array(list(data_all[0]["chemistry"].keys()))
            self.element_list = {}
            for element in elements:
                self.element_list[element] = []
                for index, chem_dict in enumerate(self.chemistry, start=0):
                    self.element_list[element].append(abs(chem_dict[element]*100))
            # if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
            #     self.w_element = self.element_list["Fe"]
            # elif self.mineral in ["Quartz"]:
            #     self.w_element = self.element_list["Si"]
            # elif self.mineral in ["Galena"]:
            #     self.w_element = self.element_list["Pb"]
            # elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
            #     self.w_element = self.element_list["Cu"]
            # elif self.mineral in ["Cassiterite"]:
            #     self.w_element = self.element_list["Sn"]
            # elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
            #     self.w_element = self.element_list["Ca"]
            # elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
            #     self.w_element = self.element_list["Mg"]
            # elif self.mineral in ["Halite"]:
            #     self.w_element = self.element_list["Na"]
            # elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
            #     self.w_element = self.element_list["Cr"]
            # elif self.mineral in ["Pyrolusite"]:
            #     self.w_element = self.element_list["Mn"]
            # elif self.mineral in ["Ilmenite", "Rutile"]:
            #     self.w_element = self.element_list["Ti"]
            # elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
            #     self.w_element = self.element_list["K"]
            # elif self.mineral in ["Illite", "Corundum"]:
            #     self.w_element = self.element_list["Al"]
        #
        self.var_opt_chem = tk.StringVar()
        self.list_elements = list(self.chemistry[0].keys())
        opt_list_chem = ["No Selection"]
        opt_list_chem.extend(self.list_elements)
        self.opt_chem = SE(parent=self.parent_mineral, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                           fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element",
                                                          opt_list=opt_list_chem, active_bg=self.color_acc_02,
                                                          command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
        self.gui_elements.append(self.opt_chem)
        #
        self.results = [self.molar_mass, self.rho_b, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod,
                        self.poisson, self.gamma_ray, self.photoelectricity]
        #
        self.entr_list_min = []
        self.entr_list_max = []
        self.entr_list_mean = []
        self.entr_list_std = []
        for i in range(10+len(elements)):
            self.entr_list_min.append(tk.IntVar())
            self.entr_list_max.append(tk.IntVar())
            self.entr_list_mean.append(tk.IntVar())
            self.entr_list_std.append(tk.IntVar())
        #
        i = 0
        for element in elements:
            lbl_w = SE(parent=self.parent_mineral, row_id=25+i, column_id=3, bg=self.color_bg,
               fg="black").create_label(text=str(element), relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl_w)
            if self.var_dict == False:
                self.results.append(self.element_list[i])
            else:
                self.results.append(self.element_list[element])
            i += 1
        ## Entry Table
        for i in range(10+len(elements)):
            if i < 10:
                entr_min = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]), 3))
                entr_max = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]), 3))
                entr_mean = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                             var_entr_set=round(np.mean(self.results[i]), 3))
                entr_std = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                             var_entr_set=round(np.std(self.results[i], ddof=1), 3))
                self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
            elif i >= 10:
                entr_min = SE(parent=self.parent_mineral, row_id=15+i, column_id=4, n_rows=1, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                             var_entr_set=round(np.min(self.results[i]), 2))
                entr_max = SE(parent=self.parent_mineral, row_id=15+i, column_id=5, n_rows=1, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                             var_entr_set=round(np.max(self.results[i]), 2))
                entr_mean = SE(parent=self.parent_mineral, row_id=15+i, column_id=6, n_rows=1, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                              var_entr_set=round(np.mean(self.results[i]), 2))
                entr_std = SE(parent=self.parent_mineral, row_id=15+i, column_id=7, n_rows=1, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                             var_entr_set=round(np.std(self.results[i], ddof=1), 2))
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        self.labels = [["Densitiy $\\varrho$ (g/ccm)", "Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)"],
                       ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                       ["Molar mass $M$ (g/mol)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
        self.labels_scatter = [["Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)", "Velocity Ratio $v_P/v_S$ (1)"],
                               ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                               ["Molar mass $M$ (g/mol)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
        self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
        self.data_plot_scatter = [[self.vP/1000, self.vS/1000, self.vP/self.vS], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
        self.create_3x3_histo(parent=self.parent_mineral, data=self.data_plot, row_id=2, column_id=9, n_rows=45,
                              n_columns=9, color=self.color_mineral, labels=self.labels)
    #
    def create_3x3_histo(self, parent, data, row_id, column_id, n_rows, n_columns, color, labels):
        #
        self.canvas = None
        self.fig, self.ax = plt.subplots(ncols=3, nrows=3, figsize=(9, 9), facecolor="#E9ECED")
        #self.fig, self.ax = plt.subplots(ncols=3, nrows=3, facecolor="#E9ECED")
        #
        for i in range(3):
            for j in range(3):
                self.ax[i][j].axvline(x=np.mean(data[i][j]), color="#E76F51", linewidth=3, linestyle="dashed")
                self.ax[i][j].hist(data[i][j], bins=15, color=color, edgecolor="black")
                #
                self.ax[i][j].set_xlabel(labels[i][j], fontsize="small")
                self.ax[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize="small")
                self.ax[i][j].grid(True)
                self.ax[i][j].set_axisbelow(True)
        self.fig.tight_layout()
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def create_3x3_scatter(self, parent, data_x, data, row_id, column_id, n_rows, n_columns, color, labels, xlabel):
        #
        self.canvas = None
        self.fig, self.ax = plt.subplots(ncols=3, nrows=3, figsize=(9, 9), facecolor="#E9ECED")
        #
        for i in range(3):
            for j in range(3):
                if self.var_opt_chem.get() in ["No Selection", "Select Element"]:
                    self.ax[i][j].scatter(data_x, data[i][j], color=color, edgecolor="black", alpha=0.5)
                else:
                    if self.var_opt_chem.get() in self.list_elements:
                        plot = self.ax[i][j].scatter(data_x, data[i][j], c=self.data_c, cmap="viridis",
                                                     edgecolor="black", alpha=1)
                #
                self.ax[i][j].set_xlabel(xlabel, fontsize="small")
                self.ax[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize="small")
                self.ax[i][j].grid(True)
                self.ax[i][j].set_axisbelow(True)
        self.fig.tight_layout()
        if self.var_opt_chem.get() not in ["No Selection", "Select Element"]:
            self.fig.subplots_adjust(bottom=0.125)
            cbar_ax = self.fig.add_axes([0.15, 0.05, 0.8, 0.01])
            cbar = self.fig.colorbar(plot, format="%.3f", cax=cbar_ax, orientation="horizontal",
                                     ticks=np.linspace(min(self.data_c), max(self.data_c), 10, endpoint=True))
            cbar.set_label(self.var_opt_chem.get()+" (%)", rotation=0)
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def enter_samples(self, var_entr, event):
        try:
            for lbl in self.lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in self.entr_w["chemistry"]:
                entr.grid_forget()
            self.lbl_w["chemistry"].clear()
            self.entr_w["chemistry"].clear()
        except:
            pass
        try:
            self.fig.clf()
            self.ax.cla()
            self.canvas.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas:
                self.canvas.destroy()
        except AttributeError:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        try:
            self.exp_data.clear()
        except:
            pass
        #
        del self.filename[-1]
        self.filename.append(var_entr.get())
        #
        data_all = []
        for i in range(var_entr.get()):
            # Oxides
            if self.mineral == "Quartz":
                data = Oxides(impurity="pure", data_type=True).create_quartz()
            elif self.mineral == "Magnetite":
                data = Oxides(impurity="pure").create_magnetite(dict=True)
            elif self.mineral == "Goethite":
                data = Oxides(impurity="pure", data_type=True).create_goethite()
            elif self.mineral == "Crocoite":
                data = Oxides(impurity="pure", data_type=True).create_crocoite()
            elif self.mineral == "Wulfenite":
                data = Oxides(impurity="pure", data_type=True).create_wulfenite()
            elif self.mineral == "Hematite-Group":
                data = Oxides(impurity="pure", data_type=True).create_hematite_group()
            elif self.mineral == "Rutile-Group":
                data = RutileGroup().create_rutile_group()
            elif self.mineral == "Periclase-Group":
                data = PericlaseGroup().create_periclase_group()
            elif self.mineral == "Wulfenite-Group":
                data = WulfeniteGroup().create_wulfenite_group()
            elif self.mineral == "Cassiterite":
                data = Oxides(impurity="pure", data_type=True).create_cassiterite()
            elif self.mineral == "Chromite":
                data = Oxides(impurity="pure", data_type=True).create_chromite()
            elif self.mineral == "Magnesiochromite":
                data = Oxides(impurity="pure", data_type=True).create_magnesiochromite()
            elif self.mineral == "Zincochromite":
                data = Oxides(impurity="pure", data_type=True).create_zincochromite()
            elif self.mineral == "Trevorite":
                data = Oxides(impurity="pure", data_type=True).create_trevorite()
            elif self.mineral == "Hematite":
                data = Oxides(impurity="pure").create_hematite(dict=True)
            elif self.mineral == "Pyrolusite":
                data = Oxides(impurity="pure", data_type=True).create_pyrolusite()
            elif self.mineral == "Rutile":
                data = Oxides(impurity="pure", data_type=True).create_rutile()
            elif self.mineral == "Ilmenite":
                data = Oxides(impurity="pure", data_type=True).create_ilmenite()
            elif self.mineral == "Litharge":
                data = Oxides(impurity="pure", data_type=True).create_litharge()
            elif self.mineral == "Massicot":
                data = Oxides(impurity="pure", data_type=True).create_massicot()
            elif self.mineral == "Minium":
                data = Oxides(impurity="pure", data_type=True).create_minium()
            elif self.mineral == "Plattnerite":
                data = Oxides(impurity="pure", data_type=True).create_plattnerite()
            elif self.mineral == "Scrutinyite":
                data = Oxides(impurity="pure", data_type=True).create_scrutinyite()
            elif self.mineral == "Zincite":
                data = Oxides(impurity="pure", data_type=True).create_zincite()
            elif self.mineral == "Uraninite":
                data = Oxides(impurity="pure", data_type=True).create_uraninite()
            elif self.mineral == "Corundum":
                data = Oxides(impurity="pure", data_type=True).create_corundum()
            elif self.mineral == "Aluminium Spinels":
                data = Oxides(impurity="pure", data_type=True).create_aluminium_spinel()
            elif self.mineral == "Chromium Spinels":
                data = Oxides(impurity="pure", data_type=True).create_chromium_spinel()
            elif self.mineral == "Iron Spinels":
                data = Oxides(impurity="pure", data_type=True).create_iron_spinel()
            elif self.mineral == "Columbite":
                data = Oxides(impurity="pure", data_type=True).create_columbite()
            elif self.mineral == "Tantalite":
                data = Oxides(impurity="pure", data_type=True).create_tantalite()
            elif self.mineral == "Coltan":
                data = Oxides(impurity="pure", data_type=True).create_coltan()
            elif self.mineral == "Cuprospinel":
                data = Oxides(impurity="pure", data_type=True).create_cuprospinel()
            elif self.mineral == "Jacobsite":
                data = Oxides(impurity="pure", data_type=True).create_jacobsite()
            elif self.mineral == "Magnesioferrite":
                data = Oxides(impurity="pure", data_type=True).create_magnesioferrite()
            elif self.mineral == "Franklinite":
                data = Oxides(impurity="pure", data_type=True).create_franklinite()
            elif self.mineral == "Ulvöspinel":
                data = Oxides(impurity="pure", data_type=True).create_ulvoespinel()
            elif self.mineral == "Pyrite":
                data = Sulfides(impurity="pure", data_type=True).create_pyrite()
            elif self.mineral == "Vaesite":
                data = Sulfides(impurity="pure", data_type=True).create_vaesite()
            elif self.mineral == "Cattierite":
                data = Sulfides(impurity="pure", data_type=True).create_cattierite()
            elif self.mineral == "Pyrite-Group":
                data = Sulfides(impurity="pure", data_type=True).create_pyrite_group()
            elif self.mineral == "Chalcopyrite":
                data = Sulfides(impurity="pure", data_type=True).create_chalcopyrite()
            elif self.mineral == "Galena":
                data = Sulfides(impurity="pure", data_type=True).create_galena()
            elif self.mineral == "Acanthite":
                data = Sulfides(impurity="pure", data_type=True).create_acanthite()
            elif self.mineral == "Chalcocite":
                data = Sulfides(impurity="pure", data_type=True).create_chalcocite()
            elif self.mineral == "Bornite":
                data = Sulfides(impurity="pure", data_type=True).create_bornite()
            elif self.mineral == "Sphalerite":
                data = Sulfides(impurity="pure", data_type=True).create_sphalerite()
            elif self.mineral == "Pyrrhotite":
                # data = Sulfides(impurity="pure", data_type=True).create_pyrrhotite()
                data = Sulfides(impurity="pure", data_type=True).create_pyrrhotite_group()
            elif self.mineral == "Millerite":
                data = Sulfides(impurity="pure", data_type=True).create_millerite()
            elif self.mineral == "Pentlandite":
                data = Sulfides(impurity="pure", data_type=True).create_pentlandite()
            elif self.mineral == "Covellite":
                data = Sulfides(impurity="pure", data_type=True).create_covellite()
            elif self.mineral == "Cinnabar":
                data = Sulfides(impurity="pure", data_type=True).create_cinnabar()
            elif self.mineral == "Stibnite":
                data = Sulfides(impurity="pure", data_type=True).create_stibnite()
            elif self.mineral == "Molybdenite":
                data = Sulfides(impurity="pure", data_type=True).create_molybdenite()
            elif self.mineral == "Fahlore":
                data = Sulfides(impurity="pure", data_type=True).create_fahlore()
            elif self.mineral == "Chalcopyrite-Group":
                data = Sulfides(impurity="pure", data_type=True).create_chalcopyrite_group()
            elif self.mineral == "Gallite":
                data = Sulfides(impurity="pure", data_type=True).create_gallite()
            elif self.mineral == "Roquesite":
                data = Sulfides(impurity="pure", data_type=True).create_roquesite()
            elif self.mineral == "Lenaite":
                data = Sulfides(impurity="pure", data_type=True).create_lenaite()
            elif self.mineral == "Laforetite":
                data = Sulfides(impurity="pure", data_type=True).create_laforetite()
            # Carbonates
            elif self.mineral == "Realgar":
                data = Sulfides(impurity="pure", data_type=True).create_realgar()
            elif self.mineral == "Orpiment":
                data = Sulfides(impurity="pure", data_type=True).create_orpiment()
            elif self.mineral == "Marcasite":
                data = Sulfides(impurity="pure", data_type=True).create_marcasite()
            elif self.mineral == "Calcite":
                data = Carbonates(impurity="pure", data_type=True).create_calcite()
            elif self.mineral == "Dolomite":
                data = Carbonates(impurity="pure", data_type=True).create_dolomite()
            elif self.mineral == "Magnesite":
                data = Carbonates(impurity="pure", data_type=True).create_magnesite()
            elif self.mineral == "Siderite":
                data = Carbonates(impurity="pure", data_type=True).create_siderite()
            elif self.mineral == "Rhodochrosite":
                data = Carbonates(impurity="pure", data_type=True).create_rhodochrosite()
            elif self.mineral == "Aragonite":
                data = Carbonates(impurity="pure", data_type=True).create_aragonite()
            elif self.mineral == "Cerussite":
                data = Carbonates(impurity="pure", data_type=True).create_cerussite()
            elif self.mineral == "Ankerite":
                data = Carbonates(impurity="pure", data_type=True).create_ankerite()
            elif self.mineral == "Azurite":
                data = Carbonates(impurity="pure", data_type=True).create_azurite()
            elif self.mineral == "Malachite":
                data = Carbonates(impurity="pure", data_type=True).create_malachite()
            elif self.mineral == "Ikaite":
                data = Carbonates(impurity="pure", data_type=True).create_ikaite()
            # Halides
            elif self.mineral == "Halite":
                data = Halides(impurity="pure", dict=True).create_halite()
            elif self.mineral == "Fluorite":
                data = Halides(impurity="pure", dict=True).create_fluorite()
            elif self.mineral == "Sylvite":
                data = Halides(impurity="pure", dict=True).create_sylvite()
            elif self.mineral == "Carnallite":
                data = Halides(impurity="pure", dict=True).create_carnallite()
            # Phyllosilicates
            elif self.mineral == "Illite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_illite()
            elif self.mineral == "Kaolinite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
            elif self.mineral == "Montmorillonite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
            elif self.mineral == "Chamosite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_chamosite()
            elif self.mineral == "Clinochlore":
                data = Phyllosilicates(impurity="pure", data_type=True).create_clinochlore()
            elif self.mineral == "Pennantite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_pennantite()
            elif self.mineral == "Nimite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_nimite()
            elif self.mineral == "Chlorite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
            elif self.mineral == "Vermiculite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_vermiculite()
            elif self.mineral == "Annite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_annite()
            elif self.mineral == "Phlogopite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_phlogopite()
            elif self.mineral == "Eastonite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_eastonite()
            elif self.mineral == "Siderophyllite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_siderophyllite()
            elif self.mineral == "Biotite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
            elif self.mineral == "Muscovite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
            elif self.mineral == "Glauconite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_glauconite()
            elif self.mineral == "Nontronite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_nontronite()
            elif self.mineral == "Saponite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_saponite()
            elif self.mineral == "Talc":
                data = Phyllosilicates(impurity="pure", data_type=True).create_talc()
            elif self.mineral == "Chrysotile":
                data = Phyllosilicates(impurity="pure", data_type=True).create_chrysotile()
            elif self.mineral == "Antigorite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_antigorite()
            elif self.mineral == "Pyrophyllite":
                data = Phyllosilicates(impurity="pure", data_type=True).create_pyrophyllite()
            # Tectosilicates
            elif self.mineral == "Alkalifeldspar":
                data = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
            elif self.mineral == "Scapolite":
                data = Tectosilicates(impurity="pure", data_type=True).create_scapolite()
            elif self.mineral == "Danburite":
                data = Tectosilicates(impurity="pure", data_type=True).create_danburite()
            elif self.mineral == "Nepheline":
                data = Tectosilicates(impurity="pure", data_type=True).create_nepheline()
            # Sulfates
            elif self.mineral == "Barite":
                data = Sulfates(impurity="pure", data_type=True).create_barite()
            elif self.mineral == "Anhydrite":
                data = Sulfates(impurity="pure", data_type=True).create_anhydrite()
            elif self.mineral == "Gypsum":
                data = Sulfates(impurity="pure", data_type=True).create_gypsum()
            elif self.mineral == "Celestite":
                data = Sulfates(impurity="pure", data_type=True).create_celestine()
            elif self.mineral == "Anglesite":
                data = Sulfates(impurity="pure", data_type=True).create_anglesite()
            elif self.mineral == "Hanksite":
                data = Sulfates(impurity="pure", data_type=True).create_hanksite()
            elif self.mineral == "Jarosite":
                data = Sulfates(impurity="pure", data_type=True).create_jarosite()
            elif self.mineral == "Alunite":
                data = Sulfates(impurity="pure", data_type=True).create_alunite()
            elif self.mineral == "Chalcanthite":
                data = Sulfates(impurity="pure", data_type=True).create_chalcanthite()
            elif self.mineral == "Kieserite":
                data = Sulfates(impurity="pure", data_type=True).create_kieserite()
            elif self.mineral == "Scheelite":
                data = Sulfates(impurity="pure", data_type=True).create_scheelite()
            elif self.mineral == "Hexahydrite":
                data = Sulfates(impurity="pure", data_type=True).create_hexahydrite()
            elif self.mineral == "Kainite":
                data = Sulfates(impurity="pure", data_type=True).create_kainite()
            # Nesosilicates
            elif self.mineral == "Zircon":
                data = Nesosilicates(data_type=True).create_zircon()
            elif self.mineral == "Thorite":
                data = Nesosilicates(data_type=True).create_thorite()
            elif self.mineral == "Titanite":
                data = Nesosilicates(data_type=True).create_titanite()
            elif self.mineral == "Andalusite":
                data = Nesosilicates(data_type=True).create_andalusite()
            elif self.mineral == "Kyanite":
                data = Nesosilicates(data_type=True).create_kyanite()
            elif self.mineral == "Sillimanite":
                data = Nesosilicates(data_type=True).create_sillimanite()
            elif self.mineral == "Topaz":
                data = Nesosilicates(data_type=True).create_topaz()
            elif self.mineral == "Staurolite":
                data = Nesosilicates(data_type=True).create_staurolite()
            elif self.mineral == "Forsterite":
                data = Nesosilicates(data_type=True).create_forsterite()
            elif self.mineral == "Fayalite":
                data = Nesosilicates(data_type=True).create_fayalite()
            elif self.mineral == "Tephroite":
                data = Nesosilicates(data_type=True).create_tephroite()
            elif self.mineral == "Calcio-Olivine":
                data = Nesosilicates(data_type=True).create_calcio_olivine()
            elif self.mineral == "Liebenbergite":
                data = Nesosilicates(data_type=True).create_liebenbergite()
            elif self.mineral == "Olivine":
                data = Nesosilicates(data_type=True).create_olivine()
            elif self.mineral == "Pyrope":
                data = Nesosilicates(data_type=True).create_pyrope()
            elif self.mineral == "Almandine":
                data = Nesosilicates(data_type=True).create_almandine()
            elif self.mineral == "Grossular":
                data = Nesosilicates(data_type=True).create_grossular()
            elif self.mineral == "Andradite":
                data = Nesosilicates(data_type=True).create_andradite()
            elif self.mineral == "Uvarovite":
                data = Nesosilicates(data_type=True).create_uvarovite()
            elif self.mineral == "Aluminium Garnet":
                data = Nesosilicates(data_type=True).create_aluminium_garnet()
            elif self.mineral == "Calcium Garnet":
                data = Nesosilicates(data_type=True).create_calcium_garnet()
            # Sorosilicates
            elif self.mineral == "Epidote":
                data = Sorosilicates(data_type=True).create_epidote()
            elif self.mineral == "Zoisite":
                data = Sorosilicates(data_type=True).create_zoisite()
            elif self.mineral == "Gehlenite":
                data = Sorosilicates(data_type=True).create_gehlenite()
            # Inosilicates
            elif self.mineral == "Enstatite":
                data = Inosilicates(data_type=True).create_enstatite()
            elif self.mineral == "Ferrosilite":
                data = Inosilicates(data_type=True).create_ferrosilite()
            elif self.mineral == "Donpeacorite":
                data = Inosilicates(data_type=True).create_donpeacorite()
            elif self.mineral == "Orthopyroxene":
                data = Inosilicates(data_type=True).create_orthopyroxene()
            elif self.mineral == "Diopside":
                data = Inosilicates(data_type=True).create_diopside()
            elif self.mineral == "Jadeite":
                data = Inosilicates(data_type=True).create_jadeite()
            elif self.mineral == "Aegirine":
                data = Inosilicates(data_type=True).create_aegirine()
            elif self.mineral == "Spodumene":
                data = Inosilicates(data_type=True).create_spodumene()
            elif self.mineral == "Wollastonite":
                data = Inosilicates(data_type=True).create_wollastonite()
            elif self.mineral == "Tremolite":
                data = Inosilicates(data_type=True).create_tremolite()
            elif self.mineral == "Actinolite":
                data = Inosilicates(data_type=True).create_actinolite()
            elif self.mineral == "Glaucophane":
                data = Inosilicates(data_type=True).create_glaucophane()
            elif self.mineral == "Augite":
                data = Inosilicates(data_type=True).create_augite()
            elif self.mineral == "Riebeckite":
                data = Inosilicates(data_type=True).create_riebeckite()
            elif self.mineral == "Arfvedsonite":
                data = Inosilicates(data_type=True).create_arfvedsonite()
            elif self.mineral == "Calcium Amphibole":
                data = Inosilicates(data_type=True).create_calcium_amphibole()
            elif self.mineral == "Sodium Amphibole":
                data = Inosilicates(data_type=True).create_sodium_amphibole()
            elif self.mineral == "Mg-Fe Pyroxene":
                data = Inosilicates(data_type=True).create_mg_fe_pyroxene()
            elif self.mineral == "Calcium Pyroxene":
                data = Inosilicates(data_type=True).create_calium_pyroxene()
            # Cyclosilicates
            elif self.mineral == "Beryl":
                data = Cyclosilicates(data_type=True).create_beryl()
            elif self.mineral == "Benitoite":
                data = Cyclosilicates(data_type=True).create_benitoite()
            elif self.mineral == "Cordierite":
                data = Cyclosilicates(data_type=True).create_cordierite()
            elif self.mineral == "Sekaninaite":
                data = Cyclosilicates(data_type=True).create_sekaninaite()
            elif self.mineral == "Schorl":
                data = Cyclosilicates(data_type=True).create_schorl()
            elif self.mineral == "Elbaite":
                data = Cyclosilicates(data_type=True).create_elbaite()
            elif self.mineral == "Liddicoatite":
                data = Cyclosilicates(data_type=True).create_liddicoatite()
            # Phosphates
            elif self.mineral == "Apatite-F":
                data = Phosphates(data_type=True).create_aptite_f()
            elif self.mineral == "Apatite-Cl":
                data = Phosphates(data_type=True).create_aptite_cl()
            elif self.mineral == "Apatite-OH":
                data = Phosphates(data_type=True).create_aptite_oh()
            elif self.mineral == "Apatite":
                data = Phosphates(data_type=True).create_aptite()
            # Phosphides
            elif self.mineral == "Allabogdanite":
                data = Phospides(data_type=True).create_allabogdanite()
            # Organics
            elif self.mineral == "Organic Matter":
                data = Organics(data_type=True).create_organics_matter()
            #
            self.color_mineral = "#7C9097"
            #
            data_all.append(data)
        self.exp_data.extend(data_all)
        #
        if self.var_dict == False:
            elements = np.array(data_all[0][-1])[:, 0]
            self.element_list = []
            for element in elements:
                self.element_list.append(DP(dataset=data_all).extract_element_amounts(type="mineral", element=element)*100)
            #
            self.rho_b = DP(dataset=data_all).extract_densities(type="mineral", keyword="bulk")
            self.molar_mass = DP(dataset=data_all).extract_molar_mass()
            self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="bulk")
            self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="shear")
            self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="poisson")
            self.vP = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vP")
            self.vS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vS")
            self.vPvS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vPvS")
            self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="mineral")
            self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="mineral")
            #
            if self.mineral == "Quartz":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Si")*100
            elif self.mineral == "Alkalifeldspar":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="K")*100
            elif self.mineral == "Plagioclase":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Ca")*100
        else:
            self.molar_mass = DP(dataset=data_all).extract_data(keyword="M")
            self.volume = DP(dataset=data_all).extract_data(keyword="V")
            self.rho_b = DP(dataset=data_all).extract_data(keyword="rho")
            self.vP = DP(dataset=data_all).extract_data(keyword="vP")
            self.vS = DP(dataset=data_all).extract_data(keyword="vS")
            self.vPvS = DP(dataset=data_all).extract_data(keyword="vP/vS")
            self.bulk_mod = DP(dataset=data_all).extract_data(keyword="K")
            self.shear_mod = DP(dataset=data_all).extract_data(keyword="G")
            self.youngs_mod = DP(dataset=data_all).extract_data(keyword="E")
            self.poisson = DP(dataset=data_all).extract_data(keyword="nu")
            self.gamma_ray = DP(dataset=data_all).extract_data(keyword="GR")
            self.photoelectricity = DP(dataset=data_all).extract_data(keyword="PE")
            self.chemistry = DP(dataset=data_all).extract_data(keyword="chemistry")
            elements = np.array(list(data_all[0]["chemistry"].keys()))
            self.element_list = {}
            for element in elements:
                self.element_list[element] = []
                for index, chem_dict in enumerate(self.chemistry, start=0):
                    self.element_list[element].append(abs(chem_dict[element]*100))
            # if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
            #     self.w_element = self.element_list["Fe"]
            # elif self.mineral in ["Quartz"]:
            #     self.w_element = self.element_list["Si"]
            # elif self.mineral in ["Galena"]:
            #     self.w_element = self.element_list["Pb"]
            # elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
            #     self.w_element = self.element_list["Cu"]
            # elif self.mineral in ["Cassiterite"]:
            #     self.w_element = self.element_list["Sn"]
            # elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
            #     self.w_element = self.element_list["Ca"]
            # elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
            #     self.w_element = self.element_list["Mg"]
            # elif self.mineral in ["Halite"]:
            #     self.w_element = self.element_list["Na"]
            # elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
            #     self.w_element = self.element_list["Cr"]
            # elif self.mineral in ["Ilmenite", "Rutile"]:
            #     self.w_element = self.element_list["Ti"]
            # elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
            #     self.w_element = self.element_list["K"]
            # elif self.mineral in ["Illite", "Corundum"]:
            #     self.w_element = self.element_list["Al"]
            # elif self.mineral in ["Pyrolusite"]:
            #     self.w_element = self.element_list["Mn"]
        #
        self.list_elements = list(self.chemistry[0].keys())
        opt_list_chem = ["No Selection"]
        opt_list_chem.extend(self.list_elements)
        self.opt_chem = SE(parent=self.parent_mineral, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                           fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element",
                                                          opt_list=opt_list_chem, active_bg=self.color_acc_02,
                                                          command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
        self.gui_elements.append(self.opt_chem)
        #
        self.var_opt_chem.set("Select Element")
        #
        self.results = [self.molar_mass, self.rho_b, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod,
                        self.poisson, self.gamma_ray, self.photoelectricity]
        #
        self.entr_list_min = []
        self.entr_list_max = []
        self.entr_list_mean = []
        self.entr_list_std = []
        for i in range(10+len(elements)):
            self.entr_list_min.append(tk.IntVar())
            self.entr_list_max.append(tk.IntVar())
            self.entr_list_mean.append(tk.IntVar())
            self.entr_list_std.append(tk.IntVar())
        #
        i = 0
        for element in elements:
            lbl_w = SE(parent=self.parent_mineral, row_id=25+i, column_id=3, bg=self.color_bg,
                       fg="black").create_label(text=str(element), relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl_w)
            if self.var_dict == False:
                self.results.append(self.element_list[i])
            else:
                self.results.append(self.element_list[element])
            i += 1
        ## Entry Table
        for i in range(10+len(elements)):
            if i < 10:
                entr_min = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]), 3))
                entr_max = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]), 3))
                entr_mean = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                             var_entr_set=round(np.mean(self.results[i]), 3))
                entr_std = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                   fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                             var_entr_set=round(np.std(self.results[i], ddof=1), 3))
                self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
            elif i >= 10:
                entr_min = SE(parent=self.parent_mineral, row_id=15+i, column_id=4, n_rows=1, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                             var_entr_set=round(np.min(self.results[i]), 2))
                entr_max = SE(parent=self.parent_mineral, row_id=15+i, column_id=5, n_rows=1, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                             var_entr_set=round(np.max(self.results[i]), 2))
                entr_mean = SE(parent=self.parent_mineral, row_id=15+i, column_id=6, n_rows=1, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                              var_entr_set=round(np.mean(self.results[i]), 2))
                entr_std = SE(parent=self.parent_mineral, row_id=15+i, column_id=7, n_rows=1, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                             var_entr_set=round(np.std(self.results[i], ddof=1), 2))
                #
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
        self.data_plot_scatter = [[self.vP/1000, self.vS/1000, self.vP/self.vS], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
        self.create_3x3_histo(parent=self.parent_mineral, data=self.data_plot, row_id=2, column_id=9, n_rows=45, n_columns=9,
                              color=self.color_mineral, labels=self.labels)
        #
        self.var_rb.set(0)
    #
    def change_radiobutton(self, var_rb):
        if var_rb.get() == 0:
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except AttributeError:
                pass
            #
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            try:
                plt.close("all")
            except:
                pass
            #
            #self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson], [self.molar_mass, self.gamma_ray, self.photoelectricity]]
            self.create_3x3_histo(parent=self.parent_mineral, data=self.data_plot, row_id=2, column_id=9, n_rows=45, n_columns=9,
                                  color=self.color_mineral, labels=self.labels)
            #
        elif var_rb.get() == 1:
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except AttributeError:
                pass
            #
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            try:
                plt.close("all")
            except:
                pass
            #
            data_x = np.array(self.rho_b)/1000
            xlabel = "Density $\\varrho$ g/ccm"
            self.create_3x3_scatter(parent=self.parent_mineral, data_x=data_x, data=self.data_plot_scatter, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_mineral, labels=self.labels_scatter,
                                    xlabel=xlabel)
            #
        elif var_rb.get() == 2:
            if var_rb.get() != self.var_actual_trace:
                try:
                    self.fig.clf()
                    self.ax.cla()
                    self.canvas.get_tk_widget().pack_forget()
                    for lbl in self.lbl_w["chemistry"]:
                        lbl.grid_forget()
                    for entr in self.entr_w["chemistry"]:
                        entr.grid_forget()
                    self.lbl_w["chemistry"].clear()
                    self.entr_w["chemistry"].clear()
                except AttributeError:
                    pass
                #
                try:
                    if self.canvas:
                        self.canvas.destroy()
                except AttributeError:
                    pass
                #
                try:
                    plt.close("all")
                except:
                    pass
                #
                self.var_actual_trace = var_rb.get()
                data_all = []
                for i in range(self.var_entr.get()):
                    # Oxides
                    if self.mineral == "Quartz":
                        data = Oxides(impurity="pure", data_type=True).create_quartz()
                    elif self.mineral == "Magnetite":
                        data = Oxides(impurity="pure", data_type=True).create_magnetite()
                    elif self.mineral == "Goethite":
                        data = Oxides(impurity="pure", data_type=True).create_goethite()
                    elif self.mineral == "Cassiterite":
                        data = Oxides(impurity="pure", data_type=True).create_cassiterite()
                    elif self.mineral == "Chromite":
                        data = Oxides(impurity="pure", data_type=True).create_chromite()
                    elif self.mineral == "Magnesiochromite":
                        data = Oxides(impurity="pure", data_type=True).create_magnesiochromite()
                    elif self.mineral == "Zincochromite":
                        data = Oxides(impurity="pure", data_type=True).create_zincochromite()
                    elif self.mineral == "Trevorite":
                        data = Oxides(impurity="pure", data_type=True).create_trevorite()
                    elif self.mineral == "Hematite":
                        data = Oxides(impurity="pure").create_hematite(dict=True)
                    elif self.mineral == "Pyrolusite":
                        data = Oxides(impurity="pure", data_type=True).create_pyrolusite()
                    elif self.mineral == "Rutile":
                        data = Oxides(impurity="pure", data_type=True).create_rutile()
                    elif self.mineral == "Ilmenite":
                        data = Oxides(impurity="pure", data_type=True).create_ilmenite()
                    elif self.mineral == "Litharge":
                        data = Oxides(impurity="pure", data_type=True).create_litharge()
                    elif self.mineral == "Massicot":
                        data = Oxides(impurity="pure", data_type=True).create_massicot()
                    elif self.mineral == "Minium":
                        data = Oxides(impurity="pure", data_type=True).create_minium()
                    elif self.mineral == "Plattnerite":
                        data = Oxides(impurity="pure", data_type=True).create_plattnerite()
                    elif self.mineral == "Scrutinyite":
                        data = Oxides(impurity="pure", data_type=True).create_scrutinyite()
                    elif self.mineral == "Zincite":
                        data = Oxides(impurity="pure", data_type=True).create_zincite()
                    elif self.mineral == "Uraninite":
                        data = Oxides(impurity="pure", data_type=True).create_uraninite()
                    elif self.mineral == "Corundum":
                        data = Oxides(impurity="pure", data_type=True).create_corundum()
                    elif self.mineral == "Aluminium Spinels":
                        data = Oxides(impurity="pure", data_type=True).create_aluminium_spinel()
                    elif self.mineral == "Chromium Spinels":
                        data = Oxides(impurity="pure", data_type=True).create_chromium_spinel()
                    elif self.mineral == "Iron Spinels":
                        data = Oxides(impurity="pure", data_type=True).create_iron_spinel()
                    elif self.mineral == "Cuprospinel":
                        data = Oxides(impurity="pure", data_type=True).create_cuprospinel()
                    elif self.mineral == "Jacobsite":
                        data = Oxides(impurity="pure", data_type=True).create_jacobsite()
                    elif self.mineral == "Magnesioferrite":
                        data = Oxides(impurity="pure", data_type=True).create_magnesioferrite()
                    elif self.mineral == "Franklinite":
                        data = Oxides(impurity="pure", data_type=True).create_franklinite()
                    elif self.mineral == "Ulvöspinel":
                        data = Oxides(impurity="pure", data_type=True).create_ulvoespinel()
                    elif self.mineral == "Pyrite":
                        data = Sulfides(impurity="pure", data_type=True).create_pyrite()
                    elif self.mineral == "Chalcopyrite":
                        data = Sulfides(impurity="pure", data_type=True).create_chalcopyrite()
                    elif self.mineral == "Galena":
                        data = Sulfides(impurity="pure", data_type=True).create_galena()
                    elif self.mineral == "Acanthite":
                        data = Sulfides(impurity="pure", data_type=True).create_acanthite()
                    elif self.mineral == "Chalcocite":
                        data = Sulfides(impurity="pure", data_type=True).create_chalcocite()
                    elif self.mineral == "Bornite":
                        data = Sulfides(impurity="pure", data_type=True).create_bornite()
                    elif self.mineral == "Sphalerite":
                        data = Sulfides(impurity="pure", data_type=True).create_sphalerite()
                    elif self.mineral == "Pyrrhotite":
                        data = Sulfides(impurity="pure", data_type=True).create_pyrrhotite()
                    elif self.mineral == "Millerite":
                        data = Sulfides(impurity="pure", data_type=True).create_millerite()
                    elif self.mineral == "Pentlandite":
                        data = Sulfides(impurity="pure", data_type=True).create_pentlandite()
                    elif self.mineral == "Covellite":
                        data = Sulfides(impurity="pure", data_type=True).create_covellite()
                    elif self.mineral == "Cinnabar":
                        data = Sulfides(impurity="pure", data_type=True).create_cinnabar()
                    elif self.mineral == "Stibnite":
                        data = Sulfides(impurity="pure", data_type=True).create_stibnite()
                    elif self.mineral == "Molybdenite":
                        data = Sulfides(impurity="pure", data_type=True).create_molybdenite()
                    elif self.mineral == "Fahlore":
                        data = Sulfides(impurity="pure", data_type=True).create_fahlore()
                    elif self.mineral == "Realgar":
                        data = Sulfides(impurity="pure", data_type=True).create_realgar()
                    elif self.mineral == "Orpiment":
                        data = Sulfides(impurity="pure", data_type=True).create_orpiment()
                    elif self.mineral == "Marcasite":
                        data = Sulfides(impurity="pure", data_type=True).create_marcasite()
                    elif self.mineral == "Calcite":
                        data = Carbonates(impurity="pure", data_type=True).create_calcite()
                    elif self.mineral == "Dolomite":
                        data = Carbonates(impurity="pure", data_type=True).create_dolomite()
                    elif self.mineral == "Magnesite":
                        data = Carbonates(impurity="pure", data_type=True).create_magnesite()
                    elif self.mineral == "Siderite":
                        data = Carbonates(impurity="pure", data_type=True).create_siderite()
                    elif self.mineral == "Rhodochrosite":
                        data = Carbonates(impurity="pure", data_type=True).create_rhodochrosite()
                    elif self.mineral == "Aragonite":
                        data = Carbonates(impurity="pure", data_type=True).create_aragonite()
                    elif self.mineral == "Cerussite":
                        data = Carbonates(impurity="pure", data_type=True).create_cerussite()
                    elif self.mineral == "Ankerite":
                        data = Carbonates(impurity="pure", data_type=True).create_ankerite()
                    elif self.mineral == "Azurite":
                        data = Carbonates(impurity="pure", data_type=True).create_azurite()
                    elif self.mineral == "Malachite":
                        data = Carbonates(impurity="pure", data_type=True).create_malachite()
                    # Halides
                    elif self.mineral == "Halite":
                        data = Halides(impurity="pure", dict=True).create_halite()
                    elif self.mineral == "Fluorite":
                        data = Halides(impurity="pure", dict=True).create_fluorite()
                    elif self.mineral == "Sylvite":
                        data = Halides(impurity="pure", dict=True).create_sylvite()
                    elif self.mineral == "Carnallite":
                        data = Halides(impurity="pure", dict=True).create_carnallite()
                    # Phyllosilicates
                    elif self.mineral == "Illite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_illite()
                    elif self.mineral == "Kaolinite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
                    elif self.mineral == "Montmorillonite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
                    elif self.mineral == "Chamosite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_chamosite()
                    elif self.mineral == "Clinochlore":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_clinochlore()
                    elif self.mineral == "Pennantite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_pennantite()
                    elif self.mineral == "Nimite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_nimite()
                    elif self.mineral == "Chlorite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
                    elif self.mineral == "Vermiculite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_vermiculite()
                    elif self.mineral == "Annite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_annite()
                    elif self.mineral == "Phlogopite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_phlogopite()
                    elif self.mineral == "Eastonite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_eastonite()
                    elif self.mineral == "Siderophyllite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_siderophyllite()
                    elif self.mineral == "Biotite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
                    elif self.mineral == "Muscovite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
                    elif self.mineral == "Glauconite":
                        data = Phyllosilicates(impurity="pure", data_type=True).create_glauconite()
                    # Tectosilicates
                    elif self.mineral == "Alkalifeldspar":
                        data = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
                    elif self.mineral == "Plagioclase":
                        data = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
                    elif self.mineral == "Scapolite":
                        data = Tectosilicates(impurity="pure", data_type=True).create_scapolite()
                    elif self.mineral == "Barite":
                        data = Sulfates(impurity="pure", data_type=True).create_barite()
                    elif self.mineral == "Anhydrite":
                        data = Sulfates(impurity="pure", data_type=True).create_anhydrite()
                    elif self.mineral == "Gypsum":
                        data = Sulfates(impurity="pure", data_type=True).create_gypsum()
                    elif self.mineral == "Celestite":
                        data = Sulfates(impurity="pure", data_type=True).create_celestine()
                    elif self.mineral == "Anglesite":
                        data = Sulfates(impurity="pure", data_type=True).create_anglesite()
                    elif self.mineral == "Hanksite":
                        data = Sulfates(impurity="pure", data_type=True).create_hanksite()
                    elif self.mineral == "Jarosite":
                        data = Sulfates(impurity="pure", data_type=True).create_jarosite()
                    elif self.mineral == "Alunite":
                        data = Sulfates(impurity="pure", data_type=True).create_alunite()
                    elif self.mineral == "Chalcanthite":
                        data = Sulfates(impurity="pure", data_type=True).create_chalcanthite()
                    elif self.mineral == "Kieserite":
                        data = Sulfates(impurity="pure", data_type=True).create_kieserite()
                    elif self.mineral == "Scheelite":
                        data = Sulfates(impurity="pure", data_type=True).create_scheelite()
                    elif self.mineral == "Hexahydrite":
                        data = Sulfates(impurity="pure", data_type=True).create_hexahydrite()
                    elif self.mineral == "Kainite":
                        data = Sulfates(impurity="pure", data_type=True).create_kainite()
                    # Nesosilicates
                    elif self.mineral == "Zircon":
                        data = Nesosilicates(data_type=True).create_zircon()
                    elif self.mineral == "Thorite":
                        data = Nesosilicates(data_type=True).create_thorite()
                    elif self.mineral == "Titanite":
                        data = Nesosilicates(data_type=True).create_titanite()
                    elif self.mineral == "Andalusite":
                        data = Nesosilicates(data_type=True).create_andalusite()
                    elif self.mineral == "Kyanite":
                        data = Nesosilicates(data_type=True).create_kyanite()
                    elif self.mineral == "Sillimanite":
                        data = Nesosilicates(data_type=True).create_sillimanite()
                    elif self.mineral == "Topaz":
                        data = Nesosilicates(data_type=True).create_topaz()
                    elif self.mineral == "Staurolite":
                        data = Nesosilicates(data_type=True).create_staurolite()
                    elif self.mineral == "Forsterite":
                        data = Nesosilicates(data_type=True).create_forsterite()
                    elif self.mineral == "Fayalite":
                        data = Nesosilicates(data_type=True).create_fayalite()
                    elif self.mineral == "Tephroite":
                        data = Nesosilicates(data_type=True).create_tephroite()
                    elif self.mineral == "Calcio-Olivine":
                        data = Nesosilicates(data_type=True).create_calcio_olivine()
                    elif self.mineral == "Liebenbergite":
                        data = Nesosilicates(data_type=True).create_liebenbergite()
                    elif self.mineral == "Olivine":
                        data = Nesosilicates(data_type=True).create_olivine()
                    elif self.mineral == "Pyrope":
                        data = Nesosilicates(data_type=True).create_pyrope()
                    elif self.mineral == "Almandine":
                        data = Nesosilicates(data_type=True).create_almandine()
                    elif self.mineral == "Grossular":
                        data = Nesosilicates(data_type=True).create_grossular()
                    elif self.mineral == "Andradite":
                        data = Nesosilicates(data_type=True).create_andradite()
                    elif self.mineral == "Uvarovite":
                        data = Nesosilicates(data_type=True).create_uvarovite()
                    elif self.mineral == "Aluminium Garnet":
                        data = Nesosilicates(data_type=True).create_aluminium_garnet()
                    elif self.mineral == "Calcium Garnet":
                        data = Nesosilicates(data_type=True).create_calcium_garnet()
                    # Sorosilicates
                    elif self.mineral == "Epidote":
                        data = Sorosilicates(data_type=True).create_epidote()
                    elif self.mineral == "Zoisite":
                        data = Sorosilicates(data_type=True).create_zoisite()
                    elif self.mineral == "Gehlenite":
                        data = Sorosilicates(data_type=True).create_gehlenite()
                    # Inosilicates
                    elif self.mineral == "Enstatite":
                        data = Inosilicates(data_type=True).create_enstatite()
                    elif self.mineral == "Ferrosilite":
                        data = Inosilicates(data_type=True).create_ferrosilite()
                    elif self.mineral == "Donpeacorite":
                        data = Inosilicates(data_type=True).create_donpeacorite()
                    elif self.mineral == "Orthopyroxene":
                        data = Inosilicates(data_type=True).create_orthopyroxene()
                    elif self.mineral == "Diopside":
                        data = Inosilicates(data_type=True).create_diopside()
                    elif self.mineral == "Jadeite":
                        data = Inosilicates(data_type=True).create_jadeite()
                    elif self.mineral == "Aegirine":
                        data = Inosilicates(data_type=True).create_aegirine()
                    elif self.mineral == "Spodumene":
                        data = Inosilicates(data_type=True).create_spodumene()
                    elif self.mineral == "Wollastonite":
                        data = Inosilicates(data_type=True).create_wollastonite()
                    elif self.mineral == "Tremolite":
                        data = Inosilicates(data_type=True).create_tremolite()
                    elif self.mineral == "Actinolite":
                        data = Inosilicates(data_type=True).create_actinolite()
                    elif self.mineral == "Glaucophane":
                        data = Inosilicates(data_type=True).create_glaucophane()
                    elif self.mineral == "Augite":
                        data = Inosilicates(data_type=True).create_augite()
                    elif self.mineral == "Riebeckite":
                        data = Inosilicates(data_type=True).create_riebeckite()
                    elif self.mineral == "Arfvedsonite":
                        data = Inosilicates(data_type=True).create_arfvedsonite()
                    elif self.mineral == "Calcium Amphibole":
                        data = Inosilicates(data_type=True).create_calcium_amphibole()
                    elif self.mineral == "Sodium Amphibole":
                        data = Inosilicates(data_type=True).create_sodium_amphibole()
                    elif self.mineral == "Mg-Fe Pyroxene":
                        data = Inosilicates(data_type=True).create_mg_fe_pyroxene()
                    elif self.mineral == "Calcium Pyroxene":
                        data = Inosilicates(data_type=True).create_calium_pyroxene()
                    # Phosphates
                    elif self.mineral == "Apatite-F":
                        data = Phosphates(data_type=True).create_aptite_f()
                    elif self.mineral == "Apatite-Cl":
                        data = Phosphates(data_type=True).create_aptite_cl()
                    elif self.mineral == "Apatite-OH":
                        data = Phosphates(data_type=True).create_aptite_oh()
                    elif self.mineral == "Apatite":
                        data = Phosphates(data_type=True).create_aptite()
                    #
                    self.color_mineral = "#7C9097"
                    #
                    data_all.append(data)
                #
                if self.var_dict == False:
                    elements = np.array(data_all[0][-1])[:, 0]
                    self.element_list = []
                    for element in elements:
                        self.element_list.append(DP(dataset=data_all).extract_element_amounts(type="mineral", element=element)*100)
                    #
                    self.rho_b = DP(dataset=data_all).extract_densities(type="mineral", keyword="bulk")
                    self.molar_mass = DP(dataset=data_all).extract_molar_mass()
                    self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="bulk")
                    self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="shear")
                    self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="poisson")
                    self.vP = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vP")
                    self.vS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vS")
                    self.vPvS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vPvS")
                    self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="mineral")
                    self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="mineral")
                    #
                    if self.mineral == "Quartz":
                        self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Si")*100
                    elif self.mineral == "Alkalifeldspar":
                        self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="K")*100
                    elif self.mineral == "Plagioclase":
                        self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Ca")*100
                else:
                    self.molar_mass = DP(dataset=data_all).extract_data(keyword="M")
                    self.volume = DP(dataset=data_all).extract_data(keyword="V")
                    self.rho_b = DP(dataset=data_all).extract_data(keyword="rho")
                    self.vP = DP(dataset=data_all).extract_data(keyword="vP")
                    self.vS = DP(dataset=data_all).extract_data(keyword="vS")
                    self.vPvS = DP(dataset=data_all).extract_data(keyword="vP/vS")
                    self.bulk_mod = DP(dataset=data_all).extract_data(keyword="K")
                    self.shear_mod = DP(dataset=data_all).extract_data(keyword="G")
                    self.youngs_mod = DP(dataset=data_all).extract_data(keyword="E")
                    self.poisson = DP(dataset=data_all).extract_data(keyword="nu")
                    self.gamma_ray = DP(dataset=data_all).extract_data(keyword="GR")
                    self.photoelectricity = DP(dataset=data_all).extract_data(keyword="PE")
                    self.chemistry = DP(dataset=data_all).extract_data(keyword="chemistry")
                    elements = np.array(list(data_all[0]["chemistry"].keys()))
                    self.element_list = {}
                    for element in elements:
                        self.element_list[element] = []
                        for index, chem_dict in enumerate(self.chemistry, start=0):
                            self.element_list[element].append(abs(chem_dict[element]*100))
                    # if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
                    #     self.w_element = self.element_list["Fe"]
                    # elif self.mineral in ["Quartz"]:
                    #     self.w_element = self.element_list["Si"]
                    # elif self.mineral in ["Galena"]:
                    #     self.w_element = self.element_list["Pb"]
                    # elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
                    #     self.w_element = self.element_list["Cu"]
                    # elif self.mineral in ["Cassiterite"]:
                    #     self.w_element = self.element_list["Sn"]
                    # elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
                    #     self.w_element = self.element_list["Ca"]
                    # elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
                    #     self.w_element = self.element_list["Mg"]
                    # elif self.mineral in ["Halite"]:
                    #     self.w_element = self.element_list["Na"]
                    # elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
                    #     self.w_element = self.element_list["Cr"]
                    # elif self.mineral in ["Ilmenite", "Rutile"]:
                    #     self.w_element = self.element_list["Ti"]
                    # elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
                    #     self.w_element = self.element_list["K"]
                    # elif self.mineral in ["Illite", "Corundum"]:
                    #     self.w_element = self.element_list["Al"]
                    # elif self.mineral in ["Pyrolusite"]:
                    #     self.w_element = self.element_list["Mn"]
                #
                self.list_elements = list(self.chemistry[0].keys())
                opt_list_chem = ["No Selection"]
                opt_list_chem.extend(self.list_elements)
                self.opt_chem = SE(parent=self.parent_mineral, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                                   fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element",
                                                                  opt_list=opt_list_chem, active_bg=self.color_acc_02,
                                                                  command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
                self.gui_elements.append(self.opt_chem)
                #
                self.var_opt_chem.set("Select Element")
                #
                self.results = [self.molar_mass, self.rho_b, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod,
                                self.poisson, self.gamma_ray, self.photoelectricity]
                #
                self.entr_list_min = []
                self.entr_list_max = []
                self.entr_list_mean = []
                self.entr_list_std = []
                for i in range(10+len(elements)):
                    self.entr_list_min.append(tk.IntVar())
                    self.entr_list_max.append(tk.IntVar())
                    self.entr_list_mean.append(tk.IntVar())
                    self.entr_list_std.append(tk.IntVar())
                #
                i = 0
                for element in elements:
                    lbl_w = SE(parent=self.parent_mineral, row_id=25+i, column_id=3, bg=self.color_bg,
                       fg="black").create_label(text=str(element), relief=tk.RAISED)
                    self.lbl_w.append(lbl_w)
                    self.lbl_chem.append(lbl_w)
                    if self.var_dict == False:
                        self.results.append(self.element_list[i])
                    else:
                        self.results.append(self.element_list[element])
                    i += 1
                ## Entry Table
                for i in range(10+len(elements)):
                    if i < 10:
                        entr_min = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                           fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]), 3))
                        entr_max = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                           fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]), 3))
                        entr_mean = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                           fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                     var_entr_set=round(np.mean(self.results[i]), 3))
                        entr_std = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                           fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                     var_entr_set=round(np.std(self.results[i], ddof=1), 3))
                    elif i >= 10:
                        entr_min = SE(parent=self.parent_mineral, row_id=15+i, column_id=4, n_rows=1, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                                     var_entr_set=round(np.min(self.results[i]), 2))
                        entr_max = SE(parent=self.parent_mineral, row_id=15+i, column_id=5, n_rows=1, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                                     var_entr_set=round(np.max(self.results[i]), 2))
                        entr_mean = SE(parent=self.parent_mineral, row_id=15+i, column_id=6, n_rows=1, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                                      var_entr_set=round(np.mean(self.results[i]), 2))
                        entr_std = SE(parent=self.parent_mineral, row_id=15+i, column_id=7, n_rows=1, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                                     var_entr_set=round(np.std(self.results[i], ddof=1), 2))
                        #
                        self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                #
                self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
                self.data_plot_scatter = [[self.vP/1000, self.vS/1000, self.vP/self.vS], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
                self.create_3x3_histo(parent=self.parent_mineral, data=self.data_plot, row_id=2, column_id=9, n_rows=45, n_columns=9,
                                      color=self.color_mineral, labels=self.labels)
                #
                self.var_rb.set(0)
            else:
                print("Nothing to do!")
        elif var_rb.get() == 3:
            #
            try:
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except:
                pass
            #
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            #
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            try:
                plt.close("all")
            except:
                pass
            #
            self.var_actual_trace = var_rb.get()
            #
            if self.mineral == "Quartz":
                traces = ["Ti", "Ge", "Sn", "Al", "Fe", "Ga", "As", "B", "H", "Li", "Na", "Ag", "K", "Mg", "Cu", "Be",
                          "Mn"]
                n_trace = rd.randint(1, len(traces))
                selection_trace = rd.sample(traces, n_trace)
            elif self.mineral == "Apatite-F":
                traces = ["Ti", "Zr", "Hf", "Th", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Cr", "As",
                          "Cl", "H", "Rb", "Mn", "Co", "Sr", "Ba", "Pb"]
                n_trace = rd.randint(1, len(traces))
                selection_trace = rd.sample(traces, n_trace)
            data_all = []
            for i in range(self.var_entr.get()):
                # Oxides
                if self.mineral == "Quartz":
                    data = Oxides(impurity=selection_trace, data_type=True).create_quartz()
                elif self.mineral == "Magnetite":
                    data = Oxides(impurity="pure", data_type=True).create_magnetite()
                elif self.mineral == "Goethite":
                    data = Oxides(impurity="pure", data_type=True).create_goethite()
                elif self.mineral == "Cassiterite":
                    data = Oxides(impurity="pure", data_type=True).create_cassiterite()
                elif self.mineral == "Chromite":
                    data = Oxides(impurity="pure", data_type=True).create_chromite()
                elif self.mineral == "Magnesiochromite":
                    data = Oxides(impurity="pure", data_type=True).create_magnesiochromite()
                elif self.mineral == "Zincochromite":
                    data = Oxides(impurity="pure", data_type=True).create_zincochromite()
                elif self.mineral == "Trevorite":
                    data = Oxides(impurity="pure", data_type=True).create_trevorite()
                elif self.mineral == "Hematite":
                    data = Oxides(impurity="pure").create_hematite(dict=True)
                elif self.mineral == "Pyrolusite":
                    data = Oxides(impurity="pure", data_type=True).create_pyrolusite()
                elif self.mineral == "Rutile":
                    data = Oxides(impurity="pure", data_type=True).create_rutile()
                elif self.mineral == "Ilmenite":
                    data = Oxides(impurity="pure", data_type=True).create_ilmenite()
                elif self.mineral == "Litharge":
                    data = Oxides(impurity="pure", data_type=True).create_litharge()
                elif self.mineral == "Massicot":
                    data = Oxides(impurity="pure", data_type=True).create_massicot()
                elif self.mineral == "Minium":
                    data = Oxides(impurity="pure", data_type=True).create_minium()
                elif self.mineral == "Plattnerite":
                    data = Oxides(impurity="pure", data_type=True).create_plattnerite()
                elif self.mineral == "Scrutinyite":
                    data = Oxides(impurity="pure", data_type=True).create_scrutinyite()
                elif self.mineral == "Zincite":
                    data = Oxides(impurity="pure", data_type=True).create_zincite()
                elif self.mineral == "Uraninite":
                    data = Oxides(impurity="pure", data_type=True).create_uraninite()
                elif self.mineral == "Corundum":
                    data = Oxides(impurity="pure", data_type=True).create_corundum()
                elif self.mineral == "Aluminium Spinels":
                    data = Oxides(impurity="pure", data_type=True).create_aluminium_spinel()
                elif self.mineral == "Chromium Spinels":
                    data = Oxides(impurity="pure", data_type=True).create_chromium_spinel()
                elif self.mineral == "Iron Spinels":
                    data = Oxides(impurity="pure", data_type=True).create_iron_spinel()
                elif self.mineral == "Cuprospinel":
                    data = Oxides(impurity="pure", data_type=True).create_cuprospinel()
                elif self.mineral == "Jacobsite":
                    data = Oxides(impurity="pure", data_type=True).create_jacobsite()
                elif self.mineral == "Magnesioferrite":
                    data = Oxides(impurity="pure", data_type=True).create_magnesioferrite()
                elif self.mineral == "Franklinite":
                    data = Oxides(impurity="pure", data_type=True).create_franklinite()
                elif self.mineral == "Ulvöspinel":
                    data = Oxides(impurity="pure", data_type=True).create_ulvoespinel()
                elif self.mineral == "Pyrite":
                    data = Sulfides(impurity="pure", data_type=True).create_pyrite()
                elif self.mineral == "Chalcopyrite":
                    data = Sulfides(impurity="pure", data_type=True).create_chalcopyrite()
                elif self.mineral == "Galena":
                    data = Sulfides(impurity="pure", data_type=True).create_galena()
                elif self.mineral == "Acanthite":
                    data = Sulfides(impurity="pure", data_type=True).create_acanthite()
                elif self.mineral == "Chalcocite":
                    data = Sulfides(impurity="pure", data_type=True).create_chalcocite()
                elif self.mineral == "Bornite":
                    data = Sulfides(impurity="pure", data_type=True).create_bornite()
                elif self.mineral == "Sphalerite":
                    data = Sulfides(impurity="pure", data_type=True).create_sphalerite()
                elif self.mineral == "Pyrrhotite":
                    data = Sulfides(impurity="pure", data_type=True).create_pyrrhotite()
                elif self.mineral == "Millerite":
                    data = Sulfides(impurity="pure", data_type=True).create_millerite()
                elif self.mineral == "Pentlandite":
                    data = Sulfides(impurity="pure", data_type=True).create_pentlandite()
                elif self.mineral == "Covellite":
                    data = Sulfides(impurity="pure", data_type=True).create_covellite()
                elif self.mineral == "Cinnabar":
                    data = Sulfides(impurity="pure", data_type=True).create_cinnabar()
                elif self.mineral == "Stibnite":
                    data = Sulfides(impurity="pure", data_type=True).create_stibnite()
                elif self.mineral == "Molybdenite":
                    data = Sulfides(impurity="pure", data_type=True).create_molybdenite()
                elif self.mineral == "Fahlore":
                    data = Sulfides(impurity="pure", data_type=True).create_fahlore()
                elif self.mineral == "Realgar":
                    data = Sulfides(impurity="pure", data_type=True).create_realgar()
                elif self.mineral == "Orpiment":
                    data = Sulfides(impurity="pure", data_type=True).create_orpiment()
                elif self.mineral == "Marcasite":
                    data = Sulfides(impurity="pure", data_type=True).create_marcasite()
                elif self.mineral == "Calcite":
                    data = Carbonates(impurity="pure", data_type=True).create_calcite()
                elif self.mineral == "Dolomite":
                    data = Carbonates(impurity="pure", data_type=True).create_dolomite()
                elif self.mineral == "Magnesite":
                    data = Carbonates(impurity="pure", data_type=True).create_magnesite()
                elif self.mineral == "Siderite":
                    data = Carbonates(impurity="pure", data_type=True).create_siderite()
                elif self.mineral == "Rhodochrosite":
                    data = Carbonates(impurity="pure", data_type=True).create_rhodochrosite()
                elif self.mineral == "Aragonite":
                    data = Carbonates(impurity="pure", data_type=True).create_aragonite()
                elif self.mineral == "Cerussite":
                    data = Carbonates(impurity="pure", data_type=True).create_cerussite()
                elif self.mineral == "Ankerite":
                    data = Carbonates(impurity="pure", data_type=True).create_ankerite()
                elif self.mineral == "Azurite":
                    data = Carbonates(impurity="pure", data_type=True).create_azurite()
                elif self.mineral == "Malachite":
                    data = Carbonates(impurity="pure", data_type=True).create_malachite()
                # Halides
                elif self.mineral == "Halite":
                    data = Halides(impurity="pure", dict=True).create_halite()
                elif self.mineral == "Fluorite":
                    data = Halides(impurity="pure", dict=True).create_fluorite()
                elif self.mineral == "Sylvite":
                    data = Halides(impurity="pure", dict=True).create_sylvite()
                elif self.mineral == "Carnallite":
                    data = Halides(impurity="pure", dict=True).create_carnallite()
                # Phyllosilicates
                elif self.mineral == "Illite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_illite()
                elif self.mineral == "Kaolinite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()
                elif self.mineral == "Montmorillonite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_montmorillonite()
                elif self.mineral == "Chamosite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_chamosite()
                elif self.mineral == "Clinochlore":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_clinochlore()
                elif self.mineral == "Pennantite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_pennantite()
                elif self.mineral == "Nimite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_nimite()
                elif self.mineral == "Chlorite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_chlorite()
                elif self.mineral == "Vermiculite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_vermiculite()
                elif self.mineral == "Annite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_annite()
                elif self.mineral == "Phlogopite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_phlogopite()
                elif self.mineral == "Eastonite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_eastonite()
                elif self.mineral == "Siderophyllite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_siderophyllite()
                elif self.mineral == "Biotite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_biotite()
                elif self.mineral == "Muscovite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_muscovite()
                elif self.mineral == "Glauconite":
                    data = Phyllosilicates(impurity="pure", data_type=True).create_glauconite()
                # Tectosilicates
                elif self.mineral == "Alkalifeldspar":
                    data = Tectosilicates(impurity="pure", data_type=True).create_alkalifeldspar()
                elif self.mineral == "Plagioclase":
                    data = Tectosilicates(impurity="pure", data_type=True).create_plagioclase()
                elif self.mineral == "Scapolite":
                    data = Tectosilicates(impurity="pure", data_type=True).create_scapolite()
                elif self.mineral == "Barite":
                    data = Sulfates(impurity="pure", data_type=True).create_barite()
                elif self.mineral == "Anhydrite":
                    data = Sulfates(impurity="pure", data_type=True).create_anhydrite()
                elif self.mineral == "Gypsum":
                    data = Sulfates(impurity="pure", data_type=True).create_gypsum()
                elif self.mineral == "Celestite":
                    data = Sulfates(impurity="pure", data_type=True).create_celestine()
                elif self.mineral == "Anglesite":
                    data = Sulfates(impurity="pure", data_type=True).create_anglesite()
                elif self.mineral == "Hanksite":
                    data = Sulfates(impurity="pure", data_type=True).create_hanksite()
                elif self.mineral == "Jarosite":
                    data = Sulfates(impurity="pure", data_type=True).create_jarosite()
                elif self.mineral == "Alunite":
                    data = Sulfates(impurity="pure", data_type=True).create_alunite()
                elif self.mineral == "Chalcanthite":
                    data = Sulfates(impurity="pure", data_type=True).create_chalcanthite()
                elif self.mineral == "Kieserite":
                    data = Sulfates(impurity="pure", data_type=True).create_kieserite()
                elif self.mineral == "Scheelite":
                    data = Sulfates(impurity="pure", data_type=True).create_scheelite()
                elif self.mineral == "Hexahydrite":
                    data = Sulfates(impurity="pure", data_type=True).create_hexahydrite()
                elif self.mineral == "Kainite":
                    data = Sulfates(impurity="pure", data_type=True).create_kainite()
                # Nesosilicates
                elif self.mineral == "Zircon":
                    data = Nesosilicates(data_type=True).create_zircon()
                elif self.mineral == "Thorite":
                    data = Nesosilicates(data_type=True).create_thorite()
                elif self.mineral == "Titanite":
                        data = Nesosilicates(data_type=True).create_titanite()
                elif self.mineral == "Andalusite":
                    data = Nesosilicates(data_type=True).create_andalusite()
                elif self.mineral == "Kyanite":
                    data = Nesosilicates(data_type=True).create_kyanite()
                elif self.mineral == "Sillimanite":
                    data = Nesosilicates(data_type=True).create_sillimanite()
                elif self.mineral == "Topaz":
                    data = Nesosilicates(data_type=True).create_topaz()
                elif self.mineral == "Staurolite":
                    data = Nesosilicates(data_type=True).create_staurolite()
                elif self.mineral == "Forsterite":
                    data = Nesosilicates(data_type=True).create_forsterite()
                elif self.mineral == "Fayalite":
                    data = Nesosilicates(data_type=True).create_fayalite()
                elif self.mineral == "Tephroite":
                    data = Nesosilicates(data_type=True).create_tephroite()
                elif self.mineral == "Calcio-Olivine":
                    data = Nesosilicates(data_type=True).create_calcio_olivine()
                elif self.mineral == "Liebenbergite":
                    data = Nesosilicates(data_type=True).create_liebenbergite()
                elif self.mineral == "Olivine":
                    data = Nesosilicates(data_type=True).create_olivine()
                elif self.mineral == "Pyrope":
                    data = Nesosilicates(data_type=True).create_pyrope()
                elif self.mineral == "Almandine":
                    data = Nesosilicates(data_type=True).create_almandine()
                elif self.mineral == "Grossular":
                    data = Nesosilicates(data_type=True).create_grossular()
                elif self.mineral == "Andradite":
                    data = Nesosilicates(data_type=True).create_andradite()
                elif self.mineral == "Uvarovite":
                    data = Nesosilicates(data_type=True).create_uvarovite()
                elif self.mineral == "Aluminium Garnet":
                    data = Nesosilicates(data_type=True).create_aluminium_garnet()
                elif self.mineral == "Calcium Garnet":
                    data = Nesosilicates(data_type=True).create_calcium_garnet()
                # Sorosilicates
                elif self.mineral == "Epidote":
                    data = Sorosilicates(data_type=True).create_epidote()
                elif self.mineral == "Zoisite":
                    data = Sorosilicates(data_type=True).create_zoisite()
                elif self.mineral == "Gehlenite":
                    data = Sorosilicates(data_type=True).create_gehlenite()
                # Inosilicates
                elif self.mineral == "Enstatite":
                    data = Inosilicates(data_type=True).create_enstatite()
                elif self.mineral == "Ferrosilite":
                    data = Inosilicates(data_type=True).create_ferrosilite()
                elif self.mineral == "Donpeacorite":
                    data = Inosilicates(data_type=True).create_donpeacorite()
                elif self.mineral == "Orthopyroxene":
                    data = Inosilicates(data_type=True).create_orthopyroxene()
                elif self.mineral == "Diopside":
                    data = Inosilicates(data_type=True).create_diopside()
                elif self.mineral == "Jadeite":
                    data = Inosilicates(data_type=True).create_jadeite()
                elif self.mineral == "Aegirine":
                    data = Inosilicates(data_type=True).create_aegirine()
                elif self.mineral == "Spodumene":
                    data = Inosilicates(data_type=True).create_spodumene()
                elif self.mineral == "Wollastonite":
                    data = Inosilicates(data_type=True).create_wollastonite()
                elif self.mineral == "Tremolite":
                    data = Inosilicates(data_type=True).create_tremolite()
                elif self.mineral == "Actinolite":
                    data = Inosilicates(data_type=True).create_actinolite()
                elif self.mineral == "Glaucophane":
                    data = Inosilicates(data_type=True).create_glaucophane()
                elif self.mineral == "Augite":
                    data = Inosilicates(data_type=True).create_augite()
                elif self.mineral == "Riebeckite":
                    data = Inosilicates(data_type=True).create_riebeckite()
                elif self.mineral == "Arfvedsonite":
                    data = Inosilicates(data_type=True).create_arfvedsonite()
                elif self.mineral == "Calcium Amphibole":
                    data = Inosilicates(data_type=True).create_calcium_amphibole()
                elif self.mineral == "Sodium Amphibole":
                    data = Inosilicates(data_type=True).create_sodium_amphibole()
                elif self.mineral == "Mg-Fe Pyroxene":
                    data = Inosilicates(data_type=True).create_mg_fe_pyroxene()
                elif self.mineral == "Calcium Pyroxene":
                    data = Inosilicates(data_type=True).create_calium_pyroxene()
                # Phosphates
                elif self.mineral == "Apatite-F":
                    data = Phosphates(impurity=selection_trace, data_type=True).create_aptite_f()
                elif self.mineral == "Apatite-Cl":
                    data = Phosphates(data_type=True).create_aptite_cl()
                elif self.mineral == "Apatite-OH":
                    data = Phosphates(data_type=True).create_aptite_oh()
                elif self.mineral == "Apatite":
                    data = Phosphates(data_type=True).create_aptite()
                #
                self.color_mineral = "#7C9097"
                #
                data_all.append(data)
            #
            if self.var_dict == False:
                elements = np.array(data_all[0][-1])[:, 0]
                self.element_list = []
                for element in elements:
                    self.element_list.append(DP(dataset=data_all).extract_element_amounts(type="mineral", element=element)*100)
                #
                self.rho_b = DP(dataset=data_all).extract_densities(type="mineral", keyword="bulk")
                self.molar_mass = DP(dataset=data_all).extract_molar_mass()
                self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="bulk")
                self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="shear")
                self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="mineral", keyword="poisson")
                self.vP = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vP")
                self.vS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vS")
                self.vPvS = DP(dataset=data_all).extract_seismic_velocities(type="mineral", keyword="vPvS")
                self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="mineral")
                self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="mineral")
                #
                if self.mineral == "Quartz":
                    self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Si")*100
                elif self.mineral == "Alkalifeldspar":
                    self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="K")*100
                elif self.mineral == "Plagioclase":
                    self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Ca")*100
            else:
                self.molar_mass = DP(dataset=data_all).extract_data(keyword="M")
                self.volume = DP(dataset=data_all).extract_data(keyword="V")
                self.rho_b = DP(dataset=data_all).extract_data(keyword="rho")
                self.vP = DP(dataset=data_all).extract_data(keyword="vP")
                self.vS = DP(dataset=data_all).extract_data(keyword="vS")
                self.vPvS = DP(dataset=data_all).extract_data(keyword="vP/vS")
                self.bulk_mod = DP(dataset=data_all).extract_data(keyword="K")
                self.shear_mod = DP(dataset=data_all).extract_data(keyword="G")
                self.youngs_mod = DP(dataset=data_all).extract_data(keyword="E")
                self.poisson = DP(dataset=data_all).extract_data(keyword="nu")
                self.gamma_ray = DP(dataset=data_all).extract_data(keyword="GR")
                self.photoelectricity = DP(dataset=data_all).extract_data(keyword="PE")
                self.chemistry = DP(dataset=data_all).extract_data(keyword="chemistry")
                elements = np.array(list(data_all[0]["chemistry"].keys()))
                self.element_list = {}
                for element in elements:
                    self.element_list[element] = []
                    for index, chem_dict in enumerate(self.chemistry, start=0):
                        self.element_list[element].append(abs(chem_dict[element]*100))
                # if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
                #     self.w_element = self.element_list["Fe"]
                # elif self.mineral in ["Quartz"]:
                #     self.w_element = self.element_list["Si"]
                # elif self.mineral in ["Galena"]:
                #     self.w_element = self.element_list["Pb"]
                # elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
                #     self.w_element = self.element_list["Cu"]
                # elif self.mineral in ["Cassiterite"]:
                #     self.w_element = self.element_list["Sn"]
                # elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
                #     self.w_element = self.element_list["Ca"]
                # elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
                #     self.w_element = self.element_list["Mg"]
                # elif self.mineral in ["Halite"]:
                #     self.w_element = self.element_list["Na"]
                # elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
                #     self.w_element = self.element_list["Cr"]
                # elif self.mineral in ["Ilmenite", "Rutile"]:
                #     self.w_element = self.element_list["Ti"]
                # elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
                #     self.w_element = self.element_list["K"]
                # elif self.mineral in ["Illite", "Corundum"]:
                #     self.w_element = self.element_list["Al"]
                # elif self.mineral in ["Pyrolusite"]:
                #     self.w_element = self.element_list["Mn"]
            #
            self.list_elements = list(self.chemistry[0].keys())
            opt_list_chem = ["No Selection"]
            opt_list_chem.extend(self.list_elements)
            self.opt_chem = SE(parent=self.parent_mineral, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                               fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element",
                                                              opt_list=opt_list_chem, active_bg=self.color_acc_02,
                                                              command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
            self.gui_elements.append(self.opt_chem)
            #
            self.var_opt_chem.set("Select Element")
            #
            self.results = [self.molar_mass, self.rho_b, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod,
                            self.poisson, self.gamma_ray, self.photoelectricity]
            #
            self.entr_list_min = []
            self.entr_list_max = []
            self.entr_list_mean = []
            self.entr_list_std = []
            for i in range(10+len(elements)):
                self.entr_list_min.append(tk.IntVar())
                self.entr_list_max.append(tk.IntVar())
                self.entr_list_mean.append(tk.IntVar())
                self.entr_list_std.append(tk.IntVar())
            #
            i = 0
            for element in elements:
                lbl_w = SE(parent=self.parent_mineral, row_id=25+i, column_id=3, bg=self.color_bg,
                           fg="black").create_label(text=str(element), relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl_w)
                if self.var_dict == False:
                    self.results.append(self.element_list[i])
                else:
                    self.results.append(self.element_list[element])
                i += 1
            ## Entry Table
            for i in range(10+len(elements)):
                if i < 10:
                    entr_min = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                                 var_entr_set=round(np.min(self.results[i]), 3))
                    entr_max = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                                 var_entr_set=round(np.max(self.results[i]), 3))
                    entr_mean = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                                  var_entr_set=round(np.mean(self.results[i]), 3))
                    entr_std = SE(parent=self.parent_mineral, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                                 var_entr_set=round(np.std(self.results[i], ddof=1), 3))
                    self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                elif i >= 10:
                    entr_min = SE(parent=self.parent_mineral, row_id=15+i, column_id=4, n_rows=1, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                                 var_entr_set=round(np.min(self.results[i]), 2))
                    entr_max = SE(parent=self.parent_mineral, row_id=15+i, column_id=5, n_rows=1, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                                 var_entr_set=round(np.max(self.results[i]), 2))
                    entr_mean = SE(parent=self.parent_mineral, row_id=15+i, column_id=6, n_rows=1, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                                  var_entr_set=round(np.mean(self.results[i]), 2))
                    entr_std = SE(parent=self.parent_mineral, row_id=15+i, column_id=7, n_rows=1, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                                 var_entr_set=round(np.std(self.results[i], ddof=1), 2))
                    self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
            #
            self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000],
                              [self.bulk_mod, self.shear_mod, self.poisson],
                              [self.molar_mass, self.gamma_ray, self.photoelectricity]]
            self.data_plot_scatter = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                      [self.bulk_mod, self.shear_mod, self.poisson],
                                      [self.molar_mass, self.gamma_ray, self.photoelectricity]]
            self.create_3x3_histo(parent=self.parent_mineral, data=self.data_plot, row_id=2, column_id=9, n_rows=45, n_columns=9,
                                  color=self.color_mineral, labels=self.labels)
            #
            self.var_rb.set(0)
    #
    def select_opt(self, var_opt):
        try:
            self.fig.clf()
            self.ax.cla()
            self.canvas.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas:
                self.canvas.destroy()
        except AttributeError:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        if var_opt == "No Selection":
            self.var_rb.set(1)
            #
            data_x = np.array(self.rho_b)/1000
            xlabel = "Density $\\varrho$ g/ccm"
            self.create_3x3_scatter(parent=self.parent_mineral, data_x=data_x, data=self.data_plot_scatter, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_mineral, labels=self.labels_scatter,
                                    xlabel=xlabel)
            #
        else:
            self.var_rb.set(1)
            #
            data_c = []
            for item in self.chemistry:
                data_c.append(item[var_opt])
            self.data_c = np.array(data_c)*100
            element = var_opt
            data_x = np.array(self.rho_b)/1000
            xlabel = "Density $\\varrho$ g/ccm"
            self.create_3x3_scatter(parent=self.parent_mineral, data_x=data_x, data=self.data_plot_scatter, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_mineral, labels=self.labels_scatter,
                                    xlabel=xlabel)
            #
    #
    def __call__(self):
        return self.lbl_w, self.entr_w, self.exp_data, self.filename
#
class Rocks:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, rock, lbl_w, entr_w, gui_elements, exp_data, filename):
        #
        try:
            for lbl in lbl_w["physics"]:
                lbl.grid_forget()
            for entr in entr_w["physics"]:
                entr.grid_forget()
            lbl_w["physics"].clear()
            entr_w["physics"].clear()
        except:
            pass
        #
        try:
            for lbl in lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in entr_w["chemistry"]:
                entr.grid_forget()
            lbl_w["chemistry"].clear()
            entr_w["chemistry"].clear()
        except:
            pass
        #
        try:
            for gui_elmnt in gui_elements:
                gui_elmnt.grid_forget()
            gui_elmnt.clear()
        except:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        try:
            exp_data.clear()
            filename.clear()
        except:
            pass
        #
        self.parent_rock = parent
        self.color_bg = color_bg
        self.color_fg = color_fg
        self.color_acc_01 = color_acc[0]
        self.color_acc_02 = color_acc[1]
        self.var_entr = tk.IntVar()
        self.var_rockname = tk.StringVar()
        var_entr_start = 100
        self.var_phi0 = tk.IntVar()
        var_phi0_start = 5
        self.var_phi1 = tk.IntVar()
        var_phi1_start = 30
        self.var_rb = tk.IntVar()
        self.var_rb.set(0)
        self.var_rb_geochem = tk.IntVar()
        self.var_rb_geochem.set(2)
        self.rock = rock
        self.lbl_w = lbl_w
        self.entr_w = entr_w
        self.lbl_chem = []
        self.entr_chem = []
        self.gui_elements = gui_elements
        self.exp_data = exp_data
        self.filename = filename
        self.filename.extend([self.rock, var_entr_start, self.var_rockname])
        self.porosities = [var_phi0_start, var_phi1_start]
        self.opt_list_chem = ["No Selection"]
        self.list_oxides = {
            "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992], "Fe": ["Fe2O3", 1.4297],
            "H": ["H2O", 8.9360], "K": ["K2O", 1.2046], "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912],
            "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725], "S": ["SO3", 2.4972], "Si": ["SiO2", 2.1392],
            "Ti": ["TiO2", 1.6681], "Zn": ["ZnO", 1.2448]}
        #
        ## Labels
        #
        lbl_stat = SE(parent=self.parent_rock, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
           fg="black").create_label(text="Statistics - "+str(self.rock), relief=tk.RAISED)
        lbl_plt = SE(parent=self.parent_rock, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
           fg="black").create_label(text="Plots - "+str(self.rock), relief=tk.RAISED)
        #
        self.gui_elements.extend([lbl_stat, lbl_plt])
        #
        lbl_01 = SE(parent=self.parent_rock, row_id=4, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Density\n (g/ccm)", relief=tk.RAISED)
        lbl_02 = SE(parent=self.parent_rock, row_id=6, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="P-wave velocity\n (km/s)", relief=tk.RAISED)
        lbl_03 = SE(parent=self.parent_rock, row_id=8, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="S-wave velocity\n (km/s)", relief=tk.RAISED)
        lbl_04 = SE(parent=self.parent_rock, row_id=10, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Velocity ratio\n (1)", relief=tk.RAISED)
        lbl_05 = SE(parent=self.parent_rock, row_id=12, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Bulk modulus\n (GPa)", relief=tk.RAISED)
        lbl_06 = SE(parent=self.parent_rock, row_id=14, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Shear modulus\n (GPa)", relief=tk.RAISED)
        lbl_07 = SE(parent=self.parent_rock, row_id=16, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Poisson's ratio\n (1)", relief=tk.RAISED)
        lbl_08 = SE(parent=self.parent_rock, row_id=18, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Porosity\n (%)", relief=tk.RAISED)
        lbl_09 = SE(parent=self.parent_rock, row_id=20, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Gamma ray\n (API)", relief=tk.RAISED)
        lbl_10 = SE(parent=self.parent_rock, row_id=22, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Photoelectricity\n (barns/electron)", relief=tk.RAISED)
        #
        self.lbl_w["physics"].extend([lbl_01, lbl_02, lbl_03, lbl_04, lbl_05, lbl_06, lbl_07, lbl_08, lbl_09, lbl_10])
        #
        if self.rock in ["Sandstone", "Marl", "Conglomerate"]:
            var_phi0_start = 5
            var_phi1_start = 30
        elif self.rock in [
            "Shale", "Kupferschiefer", "Granulite", "Greenschist", "Greenschist (basaltic)", "Greenschist (ultramafic)",
            "Greenschist (pelitic)", "Compact Hematite Ore", "Banded Iron Formation", "Amphibolite (ortho)",
            "Rhyolite (Streckeisen)", "Trachyte (Streckeisen)", "Latite (Streckeisen)", "Andesite (Streckeisen)",
            "Basalt (Streckeisen)", "Dacite (Streckeisen)", "Granite (Streckeisen)", "Granodiorite (Streckeisen)",
            "Tonalite (Streckeisen)", "Gabbro (Streckeisen)", "Diorite (Streckeisen)", "Monzonite (Streckeisen)",
            "Syenite (Streckeisen)", "Granitoid (Streckeisen)", "Quarzolite (Streckeisen)", "Mudstone",
            "Foid-bearing Syenite (Streckeisen)", "Foid-bearing Monzonite (Streckeisen)",
            "Foid-bearing Monzodiorite (Streckeisen)", "Foid-bearing Monzogabbro (Streckeisen)",
            "Foid Monzosyenite (Streckeisen)", "Foid Monzodiorite (Streckeisen)",
            "Foid Monzogabbro (Streckeisen)", "Foidolite (Streckeisen)"]:
            var_phi0_start = 0
            var_phi1_start = 10
        elif self.rock in ["Limestone", "Limestone (old)", "Dolomite Rock", "Pyroclastic Rock"]:
            var_phi0_start = 0
            var_phi1_start = 40
        elif self.rock in ["Rock Salt", "Anhydrite (Rock)", "Potash", "Felsic Rock", "Intermediate Rock", "Granite",
                           "Gabbro", "Syenite", "Diorite", "Granodiorite", "Tonalite", "Monzonite", "Quartzolite",
                           "Qz-rich Granitoid"]:
            var_phi0_start = 0
            var_phi1_start = 2.5
        else:
            var_phi0_start = 0
            var_phi1_start = 20
        #
        lbl_11 = SE(parent=self.parent_rock, row_id=29, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_label(text="Number of samples", relief=tk.RAISED)
        entr_01 = SE(parent=self.parent_rock, row_id=29, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
                     fg=color_fg).create_entry(var_entr=self.var_entr, var_entr_set=var_entr_start,
                                               command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        lbl_12 = SE(parent=self.parent_rock, row_id=30, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_label(text="Minimum Porosity (%)", relief=tk.RAISED)
        lbl_13 = SE(parent=self.parent_rock, row_id=31, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_label(text="Maximum Porosity (%)", relief=tk.RAISED)
        lbl_14 = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                    fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        #
        entr_02 = SE(parent=self.parent_rock, row_id=30, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
           fg=color_fg).create_entry(var_entr=self.var_phi0, var_entr_set=var_phi0_start,
                                     command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        entr_03 = SE(parent=self.parent_rock, row_id=31, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
           fg=color_fg).create_entry(var_entr=self.var_phi1, var_entr_set=var_phi1_start,
                                     command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        #
        self.lbl_w["physics"].extend([lbl_11, lbl_12, lbl_13, lbl_14])
        self.entr_w["physics"].extend([entr_01, entr_02, entr_03])
        #
        rb_01 = SE(parent=self.parent_rock, row_id=32, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, value_rb=0, text="Histogram",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_rock, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, value_rb=1, text="Scatter",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_03 = SE(parent=self.parent_rock, row_id=34, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, value_rb=2,
                                          text="Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        rb_04 = SE(parent=self.parent_rock, row_id=35, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, value_rb=3,
                                          text="Minerals", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        #
        self.gui_elements.extend([rb_01, rb_02, rb_03, rb_04])
        #
        if self.rock == "Custom":
            self.create_custom_rock()
        else:
            #
            data_all = []
            if self.rock not in ["Conglomerate", "Sandstone", "Anhydrite (Rock)", "Limestone"]:
                for i in range(var_entr_start):
                    if self.rock == "Sandstone":
                        data = Sandstone(fluid="water", actualThickness=0).create_sandstone(
                            number=1, porosity=rd.uniform(self.var_phi0.get() / 100, self.var_phi1.get() / 100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Conglomerate":
                        data = Sandstone(fluid="water", actualThickness=0).create_conglomerate(
                            number=1, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Shale":
                        data = shale(fluid="water").create_simple_shale(
                            dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Mudstone":
                        data = Sandstone(fluid="water", actualThickness=0).create_mudstone(
                            number=1, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Limestone (old)":
                        data = limestone(fluid="water", actualThickness=0).create_simple_limestone(
                            dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1,
                            bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Limestone":
                        data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone_alternative(
                            number=1, porosity=[self.var_phi0.get() / 100, self.var_phi1.get() / 100])
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1,
                            bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Dolomite Rock":
                        data = dolomite(fluid="water", actualThickness=0).create_simple_dolomite(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Marl":
                        data = SedimentaryRocks(fluid="water", actualThickness=0).create_marl(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    #
                    elif self.rock == "Felsic Rock":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_felsic()
                    elif self.rock == "Intermediate Rock":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_intermediate()
                    elif self.rock == "Granite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granite()
                    elif self.rock == "Gabbro":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_gabbro()
                    elif self.rock == "Diorite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_diorite()
                    elif self.rock == "Granodiorite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granodiorite()
                    elif self.rock == "Monzonite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_monzonite()
                    elif self.rock == "Quartzolite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzolite()
                    elif self.rock == "Qz-rich Granitoid":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzrich_granitoid()
                    elif self.rock == "Syenite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_syenite()
                    elif self.rock == "Tonalite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_tonalite()
                    # Plutonic Rocks
                    elif self.rock == "Granite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_granite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Granodiorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_granodiorite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Tonalite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_tonalite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Gabbro (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_gabbro_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Diorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_diorite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Monzonite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_monzonite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Syenite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_syenite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Granitoid (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_granitoid_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Quarzolite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_quarzolite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid-bearing Syenite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_syenite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid-bearing Monzonite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_monzonite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid-bearing Monzodiorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_monzodiorite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid-bearing Monzogabbro (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_monzogabbro_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid Monzosyenite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foid_monzosyenite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid Monzodiorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foid_monzodiorite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid Monzogabbro (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foid_monzogabbro_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foidolite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidolite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    #
                    # Volcanic Rocks
                    elif self.rock == "Rhyolite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_rhyolite_generalized(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Dacite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_dacite(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Trachyte (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_trachyte_generalized(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Latite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_latite_generalized(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Andesite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_andesite(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Basalt (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_basalt(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    # Pyroclastic Rocks
                    elif self.rock == "Pyroclastic Rock":
                        data = Pyroclastic(fluid="water", actualThickness=0).create_pyroclastic_rock(
                            porosity=rd.uniform(self.var_phi0.get() / 100, self.var_phi1.get() / 100))
                    #
                    ## Evaporite Rocks
                    elif self.rock == "Rock Salt":
                        data = Evaporites(fluid="water", actualThickness=0).create_simple_rocksalt(
                            dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Anhydrite (Rock)":
                        data = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(number=1,
                            porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                    elif self.rock == "Potash":
                        data = Evaporites(fluid="water", actualThickness=0).create_potash(number=1,
                            porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                    ## Ore Rocks
                    elif self.rock == "Kupferschiefer":
                        data = Ores(
                            fluid="water", actualThickness=0,
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100),
                            data_type=True).create_kupferschiefer()
                    elif self.rock == "Compact Hematite Ore":
                        data = Ores(
                            fluid="water", actualThickness=0, data_type=True).create_compact_hematite_ore(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Banded Iron Formation":
                        data = Ores(
                            fluid="water", actualThickness=0, data_type=True).create_banded_iron_formation(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    ## Metamorphic Rocks
                    elif self.rock == "Granulite":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_granulite(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist (basaltic)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist_basaltic(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist (ultramafic)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist_ultramafic(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist (pelitic)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist_pelitic(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Amphibolite (ortho)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_amphibolite_ortho(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    #
                    data_all.append(data)
                self.exp_data.extend(data_all)
            else:
                if self.rock == "Conglomerate":
                    data = Sandstone(fluid="water", actualThickness=0).create_conglomerate(
                        number=var_entr_start, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                elif self.rock == "Sandstone":
                    data = Sandstone(fluid="water", actualThickness=0).create_sandstone(
                        number=var_entr_start, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                elif self.rock == "Anhydrite (Rock)":
                    data = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                        number=var_entr_start, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                elif self.rock == "Limestone":
                    data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone_alternative(
                        number=var_entr_start, porosity=[self.var_phi0.get() / 100, self.var_phi1.get() / 100])
                #
                rb_oxides = SE(
                    parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_radiobutton(
                    var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                    color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                    self.change_radiobutton(var_rb))
                self.gui_elements.append(rb_oxides)
                for i in range(var_entr_start):
                    mineral_data = {}
                    for key, value in data["mineralogy"].items():
                        mineral_data[key] = value[i]
                    element_data = {}
                    for key, value in data["chemistry"].items():
                        element_data[key] = value[i]
                    data_set = {
                        "rock": data["rock"], "mineralogy": mineral_data, "chemistry": element_data,
                        "phi": data["phi"][i], "fluid": data["fluid"], "rho_s": data["rho_s"][i], "rho": data["rho"][i],
                        "vP": data["vP"][i], "vS": data["vS"][i], "vP/vS": data["vP/vS"][i], "K": data["K"][i],
                        "G": data["G"][i], "E": data["E"][i], "nu": data["nu"][i], "GR": data["GR"][i],
                        "PE": data["PE"][i]}
                    data_all.append(data_set)
                self.exp_data.extend(data_all)
            #
            if self.rock not in ["Conglomerate"]:
                if isinstance(data_all[0], dict):
                    self.rho = DP(dataset=data_all).extract_data(keyword="rho")
                    self.vP = DP(dataset=data_all).extract_data(keyword="vP")
                    self.vS = DP(dataset=data_all).extract_data(keyword="vS")
                    self.vPvS = DP(dataset=data_all).extract_data(keyword="vP/vS")
                    self.bulk_mod = DP(dataset=data_all).extract_data(keyword="K")
                    self.shear_mod = DP(dataset=data_all).extract_data(keyword="G")
                    self.youngs_mod = DP(dataset=data_all).extract_data(keyword="E")
                    self.poisson = DP(dataset=data_all).extract_data(keyword="nu")
                    self.phi = DP(dataset=data_all).extract_data(keyword="phi")
                    self.gamma_ray = DP(dataset=data_all).extract_data(keyword="GR")
                    self.photoelectricity = DP(dataset=data_all).extract_data(keyword="PE")
                    self.chemistry = DP(dataset=data_all).extract_data(keyword="chemistry")
                    self.mineralogy = DP(dataset=data_all).extract_data(keyword="mineralogy")
                else:
                    self.rho = DP(dataset=data_all).extract_densities(type="random", keyword="bulk")
                    self.vP = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vP")
                    self.vS = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vS")
                    self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="bulk")
                    self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="shear")
                    self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="poisson")
                    self.phi = DP(dataset=data_all).extract_porosity(type="random")
                    self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="random")
                    self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="random")
            else:
                self.rho = np.array(data["rho"])
                self.vP = np.array(data["vP"])
                self.vS = np.array(data["vS"])
                self.vPvS = np.array(data["vP/vS"])
                self.bulk_mod = np.array(data["K"])
                self.shear_mod = np.array(data["G"])
                self.youngs_mod = np.array(data["E"])
                self.poisson = np.array(data["nu"])
                self.phi = np.array(data["phi"])
                self.gamma_ray = np.array(data["GR"])
                self.photoelectricity = np.array(data["PE"])
                self.chemistry = data["chemistry"]
                self.mineralogy = data["mineralogy"]
            #
            try:
                self.list_elements = list(self.chemistry[0].keys())
                self.list_minerals = list(self.mineralogy[0].keys())
            except:
                self.list_elements = list(self.chemistry.keys())
                self.list_minerals = list(self.mineralogy.keys())
            self.elements = {}
            self.minerals = {}
            for element in self.list_elements:
                self.elements[element] = []
                try:
                    for chemistry_data in self.chemistry:
                        self.elements[element].append(abs(chemistry_data[element]*100))
                except:
                    self.elements[element].extend(np.array(self.chemistry[element])*100)
            for mineral in self.list_minerals:
                self.minerals[mineral] = []
                try:
                    for mineralogy_data in self.mineralogy:
                        self.minerals[mineral].append(abs(mineralogy_data[mineral]*100))
                except:
                    self.minerals[mineral].extend(np.array(self.mineralogy[mineral])*100)
            #
            self.results = [self.rho, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod, self.poisson, self.phi,
                            self.gamma_ray, self.photoelectricity]
            #
            self.var_opt_chem = tk.StringVar()
            self.opt_list_chem = ["No Selection"]
            self.opt_list_chem.extend(self.list_elements)
            self.opt_list_chem.extend(self.list_minerals)
            self.opt_chem = SE(parent=self.parent_rock, row_id=36, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                               fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element/Mineral",
                                                              opt_list=self.opt_list_chem, active_bg=self.color_acc_02,
                                                              command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
            self.gui_elements.append(self.opt_chem)
            #
            self.entr_list_min = []
            self.entr_list_max = []
            self.entr_list_mean = []
            self.entr_list_std = []
            if len(self.list_elements) >= len(self.list_minerals):
                for i in range(10 + len(self.list_elements)):
                    self.entr_list_min.append(tk.IntVar())
                    self.entr_list_max.append(tk.IntVar())
                    self.entr_list_mean.append(tk.IntVar())
                    self.entr_list_std.append(tk.IntVar())
            else:
                for i in range(10 + len(self.list_minerals)):
                    self.entr_list_min.append(tk.IntVar())
                    self.entr_list_max.append(tk.IntVar())
                    self.entr_list_mean.append(tk.IntVar())
                    self.entr_list_std.append(tk.IntVar())
            #
            ## Entry Table
            for i in range(10):
                if i == 7:
                    entr_min = SE(parent=self.parent_rock, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]*100), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]*100), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results[i]*100), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results[i]*100, ddof=1), 3))
                else:
                    entr_min = SE(parent=self.parent_rock, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results[i]), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results[i], ddof=1), 3))
            for index, element in enumerate(self.list_elements, start=10):
                if element not in ["U"]:
                    entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(self.elements[element]), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(self.elements[element]), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(self.elements[element]), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(self.elements[element], ddof=1), 3))
                else:
                    ppm_amounts = np.array(self.elements[element])*10000
                    entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(ppm_amounts), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(ppm_amounts), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(ppm_amounts), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                #
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
            #
            for index, element in enumerate(self.list_elements, start=0):
                if element not in ["U"]:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
            #
            self.color_rock = "#7C9097"
            #
            self.labels = [["Densitiy $\\varrho$ (g/ccm)", "Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)"],
                           ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                           ["Porosity $\\phi$ (%)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
            self.labels_scatter_rho = [["Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)", "Velocity Ratio $v_P/v_S$ (1)"],
                           ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                           ["Porosity $\\phi$ (%)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
            self.labels_scatter_phi = [["Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)", "Velocity Ratio $v_P/v_S$ (1)"],
                           ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                           ["Densitiy $\\varrho$ (g/ccm)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
            self.data_plot = [[self.rho/1000, self.vP/1000, self.vS/1000],
                              [self.bulk_mod, self.shear_mod, self.poisson],
                              [self.phi*100, self.gamma_ray, self.photoelectricity]]
            self.data_plot_scatter_rho = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                          [self.bulk_mod, self.shear_mod, self.poisson],
                                          [self.phi*100, self.gamma_ray, self.photoelectricity]]
            self.data_plot_scatter_phi = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                          [self.bulk_mod, self.shear_mod, self.poisson],
                                          [self.rho/1000, self.gamma_ray, self.photoelectricity]]
            self.create_3x3_histo(parent=self.parent_rock, data=self.data_plot, row_id=2, column_id=9, n_rows=45,
                                  n_columns=9, color=self.color_rock, labels=self.labels)
    #
    def create_3x3_histo(self, parent, data, row_id, column_id, n_rows, n_columns, color, labels):
        #
        self.canvas = None
        self.fig, self.ax = plt.subplots(ncols=3, nrows=3, figsize=(9, 9), facecolor="#E9ECED")
        #
        for i in range(3):
            for j in range(3):
                self.ax[i][j].axvline(x=np.mean(data[i][j]), color="#E76F51", linewidth=3, linestyle="dashed")
                self.ax[i][j].hist(data[i][j], bins=15, color=color, edgecolor="black")
                #
                self.ax[i][j].set_xlabel(labels[i][j], fontsize="small")
                self.ax[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize="small")
                self.ax[i][j].grid(True)
                self.ax[i][j].set_axisbelow(True)
        self.fig.tight_layout()
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                               sticky="nesw")
    #
    def create_3x3_scatter(self, parent, data_x, data, row_id, column_id, n_rows, n_columns, color, labels, xlabel):
        #
        self.canvas = None
        self.fig, self.ax = plt.subplots(ncols=3, nrows=3, figsize=(9, 9), facecolor="#E9ECED")
        colors_ukraine = ["#0057b8", "#99C9FF", "#FFF099", "#ffd700"]
        cm_ukraine = LinearSegmentedColormap.from_list("Ukraine", colors_ukraine)
        #
        for i in range(3):
            for j in range(3):
                if self.var_opt_chem.get() in ["No Selection", "Select Element/Mineral"]:
                    self.ax[i][j].scatter(data_x, data[i][j], color=color, edgecolor="black", alpha=0.5)
                else:
                    if self.var_opt_chem.get() in self.list_elements:
                        plot = self.ax[i][j].scatter(data_x, data[i][j], c=self.data_c, cmap="viridis",
                                                     edgecolor="black", alpha=1)
                    else:
                        plot = self.ax[i][j].scatter(data_x, data[i][j], c=self.data_c, cmap="viridis",
                                                     edgecolor="black", alpha=1)
                #
                self.ax[i][j].set_xlabel(xlabel, fontsize="small")
                self.ax[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize="small")
                self.ax[i][j].grid(True)
                self.ax[i][j].set_axisbelow(True)
        self.fig.tight_layout()
        if self.var_opt_chem.get() not in ["No Selection", "Select Element/Mineral"]:
            self.fig.subplots_adjust(bottom=0.125)
            cbar_ax = self.fig.add_axes([0.15, 0.05, 0.8, 0.01])
            cbar = self.fig.colorbar(plot, format="%.3f", cax=cbar_ax, orientation="horizontal",
                                     ticks=np.linspace(min(self.data_c), max(self.data_c), 10, endpoint=True))
            cbar.set_label(self.var_opt_chem.get()+" (%)", rotation=0)
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def create_plot(self, parent, data, row_id, column_id, n_rows, n_columns, xlabel, color):
        #
        self.canvas_histo = None
        self.fig_histo = Figure(facecolor="#E9ECED")
        self.ax_histo = self.fig_histo.add_subplot()
        #
        self.ax_histo.axvline(x=np.mean(data), color="#E76F51", linewidth=3, linestyle="dashed")
        self.ax_histo.hist(data, bins=15, color=color, edgecolor="black")
        self.ax_histo.grid(color="grey", linestyle="dashed", which="both")
        self.ax_histo.set_axisbelow(True)
        self.ax_histo.set_xlabel(xlabel, fontsize="small")
        self.ax_histo.set_ylabel("Frequency", labelpad=0.5, fontsize="small")
        self.fig_histo.subplots_adjust(bottom=0.15, left=0.18)
        #
        self.canvas_histo = FigureCanvasTkAgg(self.fig_histo, master=parent)
        self.canvas_histo.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                               sticky="nesw")
    #
    def create_scatter_plot(self, parent, data_x, data_y, row_id, column_id, n_rows, n_columns, xlabel, ylabel, color):
        #
        self.canvas = None
        self.fig = Figure(facecolor="#E9ECED")
        self.ax = self.fig.add_subplot()
        #
        if self.var_opt_chem.get() in ["No Selection", "Select Element/Mineral"]:
            self.ax.scatter(data_x, data_y, color=color, edgecolor="black", alpha=0.5)
        else:
            if self.var_opt_chem.get() in self.list_elements:
                plot = self.ax.scatter(data_x, data_y, c=self.elements[self.var_opt_chem.get()], cmap="viridis",
                                       edgecolor="black", alpha=1)
            elif self.var_opt_chem.get() in self.list_minerals:
                plot = self.ax.scatter(data_x, data_y, c=self.minerals[self.var_opt_chem.get()], cmap="viridis",
                                       edgecolor="black", alpha=1)
            cbar = self.fig.colorbar(plot, format="%.1f")
            cbar.set_label(self.var_opt_chem.get()+" (%)", rotation=90)
        self.ax.grid(color="grey", linestyle="dashed", which="both")
        self.ax.set_axisbelow(True)
        self.ax.set_xlim(min(data_x), max(data_x))
        self.ax.set_xticks(np.around(np.linspace(round(0.995*min(data_x), 4), 1.005*round(max(data_x), 4), 5), 2))
        self.ax.set_xlabel(xlabel, fontsize="small")
        self.ax.set_ylabel(ylabel, labelpad=0.5, fontsize="small")
        self.fig.subplots_adjust(bottom=0.15, left=0.22)
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def enter_samples(self, var_entr, event):
        try:
            for lbl in self.lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in self.entr_w["chemistry"]:
                entr.grid_forget()
            self.lbl_w["chemistry"].clear()
            self.entr_w["chemistry"].clear()
        except:
            pass
        try:
            self.fig.clf()
            self.ax.cla()
            self.canvas.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas:
                self.canvas.destroy()
        except AttributeError:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        try:
            self.exp_data.clear()
        except:
            pass
        #
        self.filename[1] = var_entr.get()
        #
        self.var_rb_geochem.set(2)
        lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
               fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        self.lbl_w["chemistry"].append(lbl)
        for index, element in enumerate(self.list_elements, start=0):
            if element not in ["U"]:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element), relief=tk.RAISED)
            else:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
        #
        try:
            data_all = []
            if self.rock not in ["Conglomerate", "Sandstone", "Anhydrite (Rock)", "Limestone"]:
                for i in range(var_entr.get()):
                    if self.rock == "Sandstone":
                        data = Sandstone(fluid="water", actualThickness=0).create_sandstone(
                            number=1, porosity=rd.uniform(self.var_phi0.get() / 100, self.var_phi1.get() / 100))
                    elif self.rock == "Conglomerate":
                        data = Sandstone(fluid="water", actualThickness=0).create_conglomerate(
                            number=1, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                    elif self.rock == "Shale":
                        data = shale(fluid="water").create_simple_shale(
                            dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Mudstone":
                        data = Sandstone(fluid="water", actualThickness=0).create_mudstone(
                            number=1, porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                    elif self.rock == "Limestone":
                        data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone_alternative(
                            number=1, porosity=[self.var_phi0.get() / 100, self.var_phi1.get() / 100])
                    elif self.rock == "Limestone (old)":
                        data = limestone(fluid="water", actualThickness=0).create_simple_limestone(
                            dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Dolomite Rock":
                        data = dolomite(fluid="water", actualThickness=0).create_simple_dolomite(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Marl":
                        data = SedimentaryRocks(fluid="water", actualThickness=0).create_marl(
                            number=1, porosity=rd.uniform(self.var_phi0.get() / 100, self.var_phi1.get() / 100))
                    #
                    elif self.rock == "Felsic Rock":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_felsic()
                    elif self.rock == "Intermediate Rock":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_intermediate()
                    elif self.rock == "Granite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granite()
                    elif self.rock == "Gabbro":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_gabbro()
                    elif self.rock == "Diorite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_diorite()
                    elif self.rock == "Granodiorite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granodiorite()
                    elif self.rock == "Monzonite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_monzonite()
                    elif self.rock == "Quartzolite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzolite()
                    elif self.rock == "Qz-rich Granitoid":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzrich_granitoid()
                    elif self.rock == "Syenite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_syenite()
                    elif self.rock == "Tonalite":
                        data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_tonalite()
                    # Plutonic Rocks
                    elif self.rock == "Granite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get() / 100, self.var_phi1.get() / 100)).create_granite_streckeisen()
                    elif self.rock == "Granodiorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_granodiorite_streckeisen()
                    elif self.rock == "Tonalite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_tonalite_streckeisen()
                    elif self.rock == "Gabbro (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_gabbro_streckeisen()
                    elif self.rock == "Diorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_diorite_streckeisen()
                    elif self.rock == "Monzonite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_monzonite_streckeisen()
                    elif self.rock == "Syenite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_syenite_streckeisen()
                    elif self.rock == "Granitoid (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_granitoid_streckeisen()
                    elif self.rock == "Quarzolite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100)).create_quarzolite_streckeisen()
                    elif self.rock == "Foid-bearing Syenite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_syenite_streckeisen()
                    elif self.rock == "Foid-bearing Monzonite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_monzonite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid-bearing Monzodiorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_monzodiorite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid-bearing Monzogabbro (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidbearing_monzogabbro_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid Monzosyenite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foid_monzosyenite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid Monzodiorite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foid_monzodiorite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foid Monzogabbro (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foid_monzogabbro_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    elif self.rock == "Foidolite (Streckeisen)":
                        data = Plutonic(
                            fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(
                                self.var_phi0.get()/100,
                                self.var_phi1.get()/100)).create_foidolite_streckeisen()
                        rb_oxides = SE(
                            parent=self.parent_rock, row_id=34, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_01,
                            fg="black").create_radiobutton(
                            var_rb=self.var_rb_geochem, value_rb=4, text="Oxides",
                            color_bg=self.color_acc_01, command=lambda var_rb=self.var_rb_geochem:
                            self.change_radiobutton(var_rb))
                        self.gui_elements.append(rb_oxides)
                    #
                    # Volcanic Rocks
                    elif self.rock == "Rhyolite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_rhyolite_generalized(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Dacite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_dacite(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Trachyte (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_trachyte_generalized(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Latite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_latite_generalized(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Andesite (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_andesite(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Basalt (Streckeisen)":
                        data = Volcanic(fluid="water", actualThickness=0).create_basalt(
                            porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    # Pyroclastic Rocks
                    elif self.rock == "Pyroclastic Rock":
                        data = Pyroclastic(fluid="water", actualThickness=0).create_pyroclastic_rock(
                            porosity=rd.uniform(self.var_phi0.get() / 100, self.var_phi1.get() / 100))
                    #
                    ## Evaporite Rocks
                    elif self.rock == "Rock Salt":
                        data = Evaporites(fluid="water", actualThickness=0).create_simple_rocksalt(
                            dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Anhydrite (Rock)":
                        data = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                            porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                    elif self.rock == "Potash":
                        data = Evaporites(fluid="water", actualThickness=0).create_potash(number=1,
                            porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                    ## Ore Rocks
                    elif self.rock == "Kupferschiefer":
                        data = Ores(
                            fluid="water", actualThickness=0, porosity=rd.uniform(
                                self.var_phi0.get()/100, self.var_phi1.get()/100),
                            data_type=True).create_kupferschiefer()
                    elif self.rock == "Compact Hematite Ore":
                        data = Ores(
                            fluid="water", actualThickness=0).create_compact_hematite_ore(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Banded Iron Formation":
                        data = Ores(
                            fluid="water", actualThickness=0, data_type=True).create_banded_iron_formation(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    ## Metamorphic Rocks
                    elif self.rock == "Granulite":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_granulite(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist (basaltic)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist_basaltic(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist (ultramafic)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist_ultramafic(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Greenschist (pelitic)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_greenschist_pelitic(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    elif self.rock == "Amphibolite (ortho)":
                        data = MetamorphicRocks(fluid="water", actualThickness=0).create_amphibolite_ortho(
                            number=1, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                    #
                    data_all.append(data)
                self.exp_data.extend(data_all)
            else:
                if self.rock == "Conglomerate":
                    data = Sandstone(fluid="water", actualThickness=0).create_conglomerate(
                        number=var_entr.get(), porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                elif self.rock == "Sandstone":
                    data = Sandstone(fluid="water", actualThickness=0).create_sandstone(
                        number=var_entr.get(), porosity=[self.var_phi0.get()/100, self.var_phi1.get()/100])
                elif self.rock == "Anhydrite (Rock)":
                    data = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(
                        number=var_entr.get(), porosity=[self.var_phi0.get() / 100, self.var_phi1.get() / 100])
                elif self.rock == "Limestone":
                    data = CarbonateRocks(fluid="water", actualThickness=0).create_limestone_alternative(
                        number=var_entr.get(), porosity=[self.var_phi0.get() / 100, self.var_phi1.get() / 100])
                for i in range(var_entr.get()):
                    mineral_data = {}
                    for key, value in data["mineralogy"].items():
                        mineral_data[key] = value[i]
                    element_data = {}
                    for key, value in data["chemistry"].items():
                        element_data[key] = value[i]
                    data_set = {
                        "rock": data["rock"], "mineralogy": mineral_data, "chemistry": element_data,
                        "phi": data["phi"][i], "fluid": data["fluid"], "rho_s": data["rho_s"][i], "rho": data["rho"][i],
                        "vP": data["vP"][i], "vS": data["vS"][i], "vP/vS": data["vP/vS"][i], "K": data["K"][i],
                        "G": data["G"][i], "E": data["E"][i], "nu": data["nu"][i], "GR": data["GR"][i],
                        "PE": data["PE"][i]}
                    data_all.append(data_set)
                self.exp_data.extend(data_all)
            #
            if self.rock not in ["Conglomerate"]:
                if isinstance(data_all[0], dict):
                    self.rho = DP(dataset=data_all).extract_data(keyword="rho")
                    self.vP = DP(dataset=data_all).extract_data(keyword="vP")
                    self.vS = DP(dataset=data_all).extract_data(keyword="vS")
                    self.vPvS = DP(dataset=data_all).extract_data(keyword="vP/vS")
                    self.bulk_mod = DP(dataset=data_all).extract_data(keyword="K")
                    self.shear_mod = DP(dataset=data_all).extract_data(keyword="G")
                    self.youngs_mod = DP(dataset=data_all).extract_data(keyword="E")
                    self.poisson = DP(dataset=data_all).extract_data(keyword="nu")
                    self.phi = DP(dataset=data_all).extract_data(keyword="phi")
                    self.gamma_ray = DP(dataset=data_all).extract_data(keyword="GR")
                    self.photoelectricity = DP(dataset=data_all).extract_data(keyword="PE")
                    self.chemistry = DP(dataset=data_all).extract_data(keyword="chemistry")
                    self.mineralogy = DP(dataset=data_all).extract_data(keyword="mineralogy")
                else:
                    self.rho = DP(dataset=data_all).extract_densities(type="random", keyword="bulk")
                    self.vP = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vP")
                    self.vS = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vS")
                    self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="bulk")
                    self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="shear")
                    self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="poisson")
                    self.phi = DP(dataset=data_all).extract_porosity(type="random")
                    self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="random")
                    self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="random")
            else:
                self.rho = np.array(data["rho"])
                self.vP = np.array(data["vP"])
                self.vS = np.array(data["vS"])
                self.vPvS = np.array(data["vP/vS"])
                self.bulk_mod = np.array(data["K"])
                self.shear_mod = np.array(data["G"])
                self.youngs_mod = np.array(data["E"])
                self.poisson = np.array(data["nu"])
                self.phi = np.array(data["phi"])
                self.gamma_ray = np.array(data["GR"])
                self.photoelectricity = np.array(data["PE"])
                self.chemistry = data["chemistry"]
                self.mineralogy = data["mineralogy"]
            #
            try:
                self.list_elements = list(self.chemistry[0].keys())
                self.list_minerals = list(self.mineralogy[0].keys())
            except:
                self.list_elements = list(self.chemistry.keys())
                self.list_minerals = list(self.mineralogy.keys())
            self.elements = {}
            self.minerals = {}
            for element in self.list_elements:
                self.elements[element] = []
                try:
                    for chemistry_data in self.chemistry:
                        self.elements[element].append(abs(chemistry_data[element]*100))
                except:
                    self.elements[element].extend(np.array(self.chemistry[element])*100)
            for mineral in self.list_minerals:
                self.minerals[mineral] = []
                try:
                    for mineralogy_data in self.mineralogy:
                        self.minerals[mineral].append(abs(mineralogy_data[mineral]*100))
                except:
                    self.minerals[mineral].extend(np.array(self.mineralogy[mineral])*100)
            #
        except:
            self.generate_custom_rock_data()
        #
        self.results = [self.rho, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod, self.poisson, self.phi,
                        self.gamma_ray, self.photoelectricity]
        #
        self.entr_list_min = []
        self.entr_list_max = []
        self.entr_list_mean = []
        self.entr_list_std = []
        if len(self.list_elements) >= len(self.list_minerals):
            for i in range(10+len(self.list_elements)):
                self.entr_list_min.append(tk.IntVar())
                self.entr_list_max.append(tk.IntVar())
                self.entr_list_mean.append(tk.IntVar())
                self.entr_list_std.append(tk.IntVar())
        else:
            for i in range(10+len(self.list_minerals)):
                self.entr_list_min.append(tk.IntVar())
                self.entr_list_max.append(tk.IntVar())
                self.entr_list_mean.append(tk.IntVar())
                self.entr_list_std.append(tk.IntVar())
        #
        ## Entry Table
        for i in range(10):
            if i == 7:
                entr_min = SE(parent=self.parent_rock, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]*100), 3))
                entr_max = SE(parent=self.parent_rock, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]*100), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results[i]*100), 3))
                entr_std = SE(parent=self.parent_rock, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results[i]*100, ddof=1), 3))
            else:
                entr_min = SE(parent=self.parent_rock, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]), 3))
                entr_max = SE(parent=self.parent_rock, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results[i]), 3))
                entr_std = SE(parent=self.parent_rock, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results[i], ddof=1), 3))
        for index, element in enumerate(self.list_elements, start=10):
            if element not in ["U"]:
                entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(self.elements[element]), 3))
                entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(self.elements[element]), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(self.elements[element]), 3))
                entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(self.elements[element], ddof=1), 3))
            else:
                ppm_amounts = np.array(self.elements[element])*10000
                entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(ppm_amounts), 3))
                entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(ppm_amounts), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(ppm_amounts), 3))
                entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
            #
            self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        self.color_rock = "#7C9097"
        #
        self.data_plot = [[self.rho/1000, self.vP/1000, self.vS/1000],
                          [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.phi*100, self.gamma_ray, self.photoelectricity]]
        self.data_plot_scatter_rho = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                      [self.bulk_mod, self.shear_mod, self.poisson],
                                      [self.phi*100, self.gamma_ray, self.photoelectricity]]
        self.data_plot_scatter_phi = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                      [self.bulk_mod, self.shear_mod, self.poisson],
                                      [self.rho/1000, self.gamma_ray, self.photoelectricity]]
        self.create_3x3_histo(parent=self.parent_rock, data=self.data_plot, row_id=2, column_id=9, n_rows=45,
                              n_columns=9, color=self.color_rock, labels=self.labels)
        #
        self.var_rb.set(0)
    #
    def change_radiobutton(self, var_rb):
        if var_rb.get() == 0:
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass

            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            self.create_3x3_histo(parent=self.parent_rock, data=self.data_plot, row_id=2, column_id=9, n_rows=45,
                                  n_columns=9, color=self.color_rock, labels=self.labels)
            #
        elif var_rb.get() == 1:
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass

            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            self.var_prop = tk.IntVar()
            self.var_prop.set(8)
            rb_rho = SE(parent=self.parent_rock, row_id=37, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, value_rb=8,
                                                       text="Density", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            rb_phi = SE(parent=self.parent_rock, row_id=38, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, value_rb=9,
                                                       text="Porosity", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            self.rb_prop = [rb_rho, rb_phi]
            self.gui_elements.extend([rb_rho, rb_phi])
            #
            data_x_rho = np.array(self.rho)/1000
            xlabel="Densitiy $\\varrho$ (g/ccm)"
            #
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x_rho, data=self.data_plot_scatter_rho, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels_scatter_rho,
                                    xlabel=xlabel)
            #
        elif var_rb.get() == 2:
            #
            try:
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except:
                pass
            #
            lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                     fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
            for index, element in enumerate(self.list_elements, start=0):
                if element not in ["U"]:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
            #
            for index, element in enumerate(self.list_elements, start=10):
                if element not in ["U"]:
                    entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(self.elements[element]), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(self.elements[element]), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(self.elements[element]), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(self.elements[element], ddof=1), 3))
                else:
                    ppm_amounts = np.array(self.elements[element])*10000
                    entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(ppm_amounts), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(ppm_amounts), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(ppm_amounts), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                #
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        elif var_rb.get() == 3:
            #
            try:
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except:
                pass
            #
            lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                   fg="black").create_label(text="Mineralogical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
            for index, mineral in enumerate(self.list_minerals, start=0):
                if mineral not in ["Urn"]:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(mineral), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(mineral)+" (ppm)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
            #
            for index, mineral in enumerate(self.list_minerals, start=10):
                if mineral not in ["Urn"]:
                    entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(self.minerals[mineral]), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(self.minerals[mineral]), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(self.minerals[mineral]), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(self.minerals[mineral], ddof=1), 3))
                else:
                    ppm_amounts = np.array(self.minerals[mineral])*10000
                    entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(ppm_amounts), 3))
                    entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(ppm_amounts), 3))
                    entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(ppm_amounts), 3))
                    entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                #
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        elif var_rb.get() == 8:
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            data_x_rho = np.array(self.rho)/1000
            xlabel = "Densitiy $\\varrho$ (g/ccm)"
            #
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x_rho, data=self.data_plot_scatter_rho, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels_scatter_rho,
                                    xlabel=xlabel)
            #
        elif var_rb.get() == 9:
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            data_x_phi = np.array(self.phi)*100
            xlabel = "Porosity $\\phi$ (%)"
            #
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x_phi, data=self.data_plot_scatter_phi, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels_scatter_phi,
                                    xlabel=xlabel)
        #
        elif var_rb.get() == 4:
            #
            try:
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except:
                pass
            #
            if self.rock in ["Sandstone"]:
                w_Chl = round(np.mean(self.minerals["Chl"]), 4)
                w_Hem = round(np.mean(self.minerals["Hem"]), 4)
                self.list_oxides = {
                    "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992],
                    "Fe": ["FeO", w_Chl, 1.2865, "Fe2O3", w_Hem, 1.4297], "H": ["H2O", 8.9360], "K": ["K2O", 1.2046],
                    "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912], "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725],
                    "S": ["SO3", 2.4972], "Si": ["SiO2", 2.1392], "Ti": ["TiO2", 1.6681], "Zn": ["ZnO", 1.2448]}
            elif self.rock in ["Conglomerate"]:
                w_Bt = round(np.mean(self.minerals["Bt"]), 4)
                w_Hem = round(np.mean(self.minerals["Hem"]), 4)
                self.list_oxides = {
                    "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992],
                    "Fe": ["FeO", w_Bt, 1.2865, "Fe2O3", w_Hem, 1.4297], "H": ["H2O", 8.9360], "K": ["K2O", 1.2046],
                    "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912], "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725],
                    "S": ["SO3", 2.4972], "Si": ["SiO2", 2.1392], "Ti": ["TiO2", 1.6681], "Zn": ["ZnO", 1.2448]}
            elif self.rock in ["Limestone"]:
                w_Ilt = round(np.mean(self.minerals["Ilt"]), 4)
                self.list_oxides = {
                    "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992],
                    "Fe": ["FeO", 1.2865], "H": ["H2O", 8.9360], "K": ["K2O", 1.2046],
                    "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912], "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725],
                    "Si": ["SiO2", 2.1392]}
            elif self.rock in ["Mudstone"]:
                self.list_oxides = {
                    "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992],
                    "Fe": ["FeO", 1.2865], "H": ["H2O", 8.9360], "K": ["K2O", 1.2046],
                    "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912], "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725],
                    "Si": ["SiO2", 2.1392]}
            elif self.rock in ["Shale"]:
                self.list_oxides = {
                    "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992],
                    "Fe": ["FeO", 1.2865], "H": ["H2O", 8.9360], "K": ["K2O", 1.2046],
                    "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912], "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725],
                    "Si": ["SiO2", 2.1392], "U": ["UO2", 1.1344], "N": ["NO", 2.1422], "F": ["F", 1.0],
                    "S": ["SO2", 1.9981]}
            elif self.rock in ["Granite (Streckeisen)", "Tonalite (Streckeisen)", "Diorite (Streckeisen)",
                               "Gabbro (Streckeisen)", "Syenite (Streckeisen)", "Quarzolite (Streckeisen)",
                               "Granodiorite (Streckeisen)", "Granitoid (Streckeisen)", "Monzonite (Streckeisen)",
                               "Rhyolite (Streckeisen)", "Dacite (Streckeisen)", "Latite (Streckeisen)",
                               "Basalt (Streckeisen)", "Andesite (Streckeisen)", "Trachyte (Streckeisen)",
                               "Foid-bearing Syenite (Streckeisen)", "Foid-bearing Monzonite (Streckeisen)",
                               "Foid-bearing Monzodiorite (Streckeisen)", "Foid-bearing Monzogabbro (Streckeisen)",
                               "Foid Monzosyenite (Streckeisen)", "Foid Monzodiorite (Streckeisen)",
                               "Foid Monzogabbro (Streckeisen)", "Foidolite (Streckeisen)"]:
                self.list_oxides = {
                    "Al": ["Al2O3", 1.8895], "C": ["CO2", 3.6644], "Ca": ["CaO", 1.3992], "Fe": ["FeO", 1.2865],
                    "H": ["H2O", 8.9360], "K": ["K2O", 1.2046], "Mg": ["MgO", 1.6582], "Mn": ["MnO", 1.2912],
                    "Na": ["Na2O", 1.3480], "Ni": ["NiO", 1.2725], "S": ["SO3", 2.4972], "Si": ["SiO2", 2.1392],
                    "Ti": ["TiO2", 1.6681], "Zn": ["ZnO", 1.2448]}
            #
            lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                     fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
            self.list_elements.insert(0, self.list_elements.pop(self.list_elements.index("O")))
            self.additional_line = 0
            for index, element in enumerate(self.list_elements, start=1):
                if element not in self.list_oxides:
                    pass
                else:
                    if len(self.list_oxides[element]) == 2:
                        lbl = SE(parent=self.parent_rock, row_id=23+index + self.additional_line, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(self.list_oxides[element][0]), relief=tk.RAISED)
                    elif len(self.list_oxides[element]) == 6:
                        lbl = SE(parent=self.parent_rock, row_id=23+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(self.list_oxides[element][0]), relief=tk.RAISED)
                        self.lbl_w["chemistry"].append(lbl)
                        lbl = SE(parent=self.parent_rock, row_id=24+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(self.list_oxides[element][3]), relief=tk.RAISED)
                        self.additional_line += 1
                self.lbl_w["chemistry"].append(lbl)
            #
            self.additional_line = 0
            for index, element in enumerate(self.list_elements, start=10):
                if element in self.list_oxides:
                    if len(self.list_oxides[element]) == 2:
                        if element != "U":
                            entr_min = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=4,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_min[index+self.additional_line],
                                var_entr_set=round(np.min(self.elements[element])*self.list_oxides[element][1], 3))
                            entr_max = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=5,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_max[index+self.additional_line],
                                var_entr_set=round(np.max(self.elements[element])*self.list_oxides[element][1], 3))
                            entr_mean = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=6,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_mean[index+self.additional_line],
                                var_entr_set=round(np.mean(self.elements[element])*self.list_oxides[element][1], 3))
                            entr_std = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=7,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_std[index+self.additional_line],
                                var_entr_set=round(
                                    np.std(self.elements[element], ddof=1)*self.list_oxides[element][1], 3))
                        else:
                            entr_min = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=4,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_min[index+self.additional_line],
                                var_entr_set=round(np.min(self.elements[element])*self.list_oxides[element][1], 6))
                            entr_max = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=5,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_max[index+self.additional_line],
                                var_entr_set=round(np.max(self.elements[element])*self.list_oxides[element][1], 6))
                            entr_mean = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=6,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_mean[index+self.additional_line],
                                var_entr_set=round(np.mean(self.elements[element])*self.list_oxides[element][1], 6))
                            entr_std = SE(
                                parent=self.parent_rock, row_id=14+index + self.additional_line, column_id=7,
                                bg=self.color_bg, fg=self.color_fg).create_entry(
                                var_entr=self.entr_list_std[index+self.additional_line],
                                var_entr_set=round(
                                    np.std(self.elements[element], ddof=1)*self.list_oxides[element][1], 6))
                    elif len(self.list_oxides[element]) == 6:
                        entr_min = SE(parent=self.parent_rock, row_id=14+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(np.array(self.elements[element])*np.array(self.list_oxides[element][1]))*self.list_oxides[element][2], 3))
                        entr_max = SE(parent=self.parent_rock, row_id=14+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(np.array(self.elements[element])*np.array(self.list_oxides[element][1]))*self.list_oxides[element][2], 3))
                        entr_mean = SE(parent=self.parent_rock, row_id=14+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(np.array(self.elements[element])*np.array(self.list_oxides[element][1]))*self.list_oxides[element][2], 3))
                        entr_std = SE(parent=self.parent_rock, row_id=14+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(np.array(self.elements[element])*np.array(self.list_oxides[element][1]), ddof=1)*self.list_oxides[element][2], 3))
                        #
                        self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                        self.entr_list_min.append(tk.IntVar())
                        self.entr_list_max.append(tk.IntVar())
                        self.entr_list_mean.append(tk.IntVar())
                        self.entr_list_std.append(tk.IntVar())
                        #
                        entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index+1],
                                                                 var_entr_set=round(np.min(np.array(self.elements[element])*np.array(self.list_oxides[element][4]))*self.list_oxides[element][5], 3))
                        entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index+1],
                                                                 var_entr_set=round(np.max(np.array(self.elements[element])*np.array(self.list_oxides[element][4]))*self.list_oxides[element][5], 3))
                        entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index+1],
                                                                  var_entr_set=round(np.mean(np.array(self.elements[element])*np.array(self.list_oxides[element][4]))*self.list_oxides[element][5], 3))
                        entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index+1],
                                                                 var_entr_set=round(np.std(np.array(self.elements[element])*np.array(self.list_oxides[element][4]), ddof=1)*self.list_oxides[element][5], 3))
                        self.additional_line += 1
                    #
                    self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
    #
    def select_opt(self, var_opt):
        try:
            self.fig.clf()
            self.ax.cla()
            self.canvas.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas:
                self.canvas.destroy()
        except AttributeError:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        self.var_prop = tk.IntVar()
        self.var_prop.set(8)
        rb_rho = SE(parent=self.parent_rock, row_id=37, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                    fg="black").create_radiobutton(var_rb=self.var_prop, value_rb=8,
                                                   text="Density", color_bg=self.color_acc_01,
                                                   command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
        rb_phi = SE(parent=self.parent_rock, row_id=38, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                    fg="black").create_radiobutton(var_rb=self.var_prop, value_rb=9,
                                                   text="Porosity", color_bg=self.color_acc_01,
                                                   command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
        self.rb_prop = [rb_rho, rb_phi]
        self.gui_elements.extend(self.rb_prop)
        #
        if var_opt == "No Selection":
            self.var_rb.set(1)
            #
            data_x = np.array(self.rho)/1000
            xlabel = "Density $\\varrho$ g/ccm"
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x, data=self.data_plot_scatter_rho, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels_scatter_rho,
                                    xlabel=xlabel)
            #
        else:
            self.var_rb.set(1)
            #
            data_c = []
            # Mineralogy
            if var_opt in list(self.elements.keys()):
                try:
                    for item in self.chemistry:
                        data_c.append(item[var_opt])
                except:
                    data_c = np.array(self.elements[var_opt])/100
            if var_opt in list(self.minerals.keys()):
                try:
                    for item in self.mineralogy:
                        data_c.append(item[var_opt])
                except:
                    data_c = np.array(self.minerals[var_opt])/100
            self.data_c = np.array(data_c)*100
            data_x = np.array(self.rho)/1000
            xlabel = "Density $\\varrho$ g/ccm"
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x, data=self.data_plot_scatter_rho, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels_scatter_rho,
                                    xlabel=xlabel)
    #
    #def __call__(self):
    #    return self.lbl_w, self.entr_w
    #
    def create_custom_rock(self):
        #
        try:
            for lbl in self.lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in self.entr_w["chemistry"]:
                entr.grid_forget()
            self.lbl_w["chemistry"].clear()
            self.entr_w["chemistry"].clear()
            for lbl in self.lbl_w["physics"]:
                lbl.grid_forget()
            for entr in self.entr_w["physics"]:
                entr.grid_forget()
            self.lbl_w["physics"].clear()
            self.entr_w["physics"].clear()
        except:
            pass
        #
        try:
            self.fig.clf()
            self.ax.cla()
            self.canvas.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas:
                self.canvas.destroy()
        except AttributeError:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        try:
            self.exp_data.clear()
        except:
            pass
        #
        ## Variables
        self.custom_mineralogy = {}
        self.custom_mineralogy["mineralogy"] = {}
        self.custom_porosities = {}
        #
        ## Buttons
        btn_defmin = SE(parent=self.parent_rock, row_id=40, column_id=0, n_rows=2, n_columns=1, bg=self.color_acc_01,
                        fg="black").create_button(text="Define Mineralogy", command=lambda var_btn="Define Mineralogy": self.press_button(var_btn))
        btn_update = SE(parent=self.parent_rock, row_id=40, column_id=1, n_rows=2, n_columns=1, bg=self.color_acc_01,
                        fg="black").create_button(text="Generate Data", command=lambda var_btn="Generate Data": self.press_button(var_btn))
        #
        self.gui_elements.extend([btn_defmin, btn_update])
        #
        lbl_11 = SE(parent=self.parent_rock, row_id=29, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_label(text="Number of samples", relief=tk.RAISED)
        entr_01 = SE(parent=self.parent_rock, row_id=29, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
                     fg="black").create_entry(var_entr=self.var_entr, var_entr_set=self.filename[1],
                                              command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        lbl_12 = SE(parent=self.parent_rock, row_id=30, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_label(text="Minimum Porosity (%)", relief=tk.RAISED)
        lbl_13 = SE(parent=self.parent_rock, row_id=31, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
                    fg="black").create_label(text="Maximum Porosity (%)", relief=tk.RAISED)
        #
        entr_02 = SE(parent=self.parent_rock, row_id=30, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
                     fg="black").create_entry(var_entr=self.var_phi0, var_entr_set=self.porosities[0],
                                              command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        entr_03 = SE(parent=self.parent_rock, row_id=31, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
                     fg="black").create_entry(var_entr=self.var_phi1, var_entr_set=self.porosities[1],
                                              command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        lbl_name = SE(parent=self.parent_rock, row_id=39, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
            fg="black").create_label(text="Rock name", relief=tk.RAISED)
        entr_name = SE(parent=self.parent_rock, row_id=39, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
                     fg="black").create_entry(var_entr=self.filename[2], var_entr_set="Custom",
                                              command=lambda event, var_entr=self.filename[2]: self.enter_samples(var_entr, event))
        entr_name.bind("<Return>", lambda event, var_name=self.filename[2]: self.change_rockname(var_name, event))
        #
        self.lbl_w["physics"].extend([lbl_11, lbl_12, lbl_13, lbl_name])
        self.entr_w["physics"].extend([entr_01, entr_02, entr_03, entr_name])
        #
        ## Labels
        #
        lbl_stat = SE(parent=self.parent_rock, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
           fg="black").create_label(text="Statistics - "+str(self.rock), relief=tk.RAISED)
        lbl_plt = SE(parent=self.parent_rock, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
           fg="black").create_label(text="Plots - "+str(self.rock), relief=tk.RAISED)
        #
        self.gui_elements.extend([lbl_stat, lbl_plt])
        #
        lbl_01 = SE(parent=self.parent_rock, row_id=4, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Density\n (kg/m3)", relief=tk.RAISED)
        lbl_02 = SE(parent=self.parent_rock, row_id=6, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="P-wave velocity\n (m/s)", relief=tk.RAISED)
        lbl_03 = SE(parent=self.parent_rock, row_id=8, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="S-wave velocity\n (m/s)", relief=tk.RAISED)
        lbl_04 = SE(parent=self.parent_rock, row_id=10, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Velocity ratio\n (1)", relief=tk.RAISED)
        lbl_05 = SE(parent=self.parent_rock, row_id=12, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Bulk modulus\n (GPa)", relief=tk.RAISED)
        lbl_06 = SE(parent=self.parent_rock, row_id=14, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Shear modulus\n (GPa)", relief=tk.RAISED)
        lbl_07 = SE(parent=self.parent_rock, row_id=16, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Poisson's ratio\n (1)", relief=tk.RAISED)
        lbl_08 = SE(parent=self.parent_rock, row_id=18, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Porosity\n (%)", relief=tk.RAISED)
        lbl_09 = SE(parent=self.parent_rock, row_id=20, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Gamma ray\n (API)", relief=tk.RAISED)
        lbl_10 = SE(parent=self.parent_rock, row_id=22, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Photoelectricity\n (barns/electron)", relief=tk.RAISED)
        #
        self.lbl_w["physics"].extend([lbl_01, lbl_02, lbl_03, lbl_04, lbl_05, lbl_06, lbl_07, lbl_08, lbl_09, lbl_10])
        #
        lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                 fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        self.lbl_w["chemistry"].append(lbl)
        #
    def change_rockname(self, var_name, event):
        var_name.set(var_name.get())
        self.filename[0] = var_name.get()
        #
        lbl_stat = SE(parent=self.parent_rock, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
           fg="black").create_label(text="Statistics - "+str(var_name.get()), relief=tk.RAISED)
        lbl_plt = SE(parent=self.parent_rock, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
           fg="black").create_label(text="Plots - "+str(var_name.get()), relief=tk.RAISED)
        #
        self.gui_elements.extend([lbl_stat, lbl_plt])
        #
    #
    def press_button(self, var_btn):
        if var_btn == "Define Mineralogy":
            if len(self.custom_mineralogy["mineralogy"]) > 0:
                self.custom_mineralogy["mineralogy"].clear()
            ## CONSTANTS
            self.color_menu = "#264653"
            self.color_border = "#7C9097"
            self.color_bg = "#E9ECED"
            self.color_fg_dark = "black"
            self.color_fg_light = "white"
            self.color_accent_01 = "#E76F51"    # Orange
            self.color_accent_02 = "#F0A794"    # Orange light
            self.color_accent_03 = "#E9C46A"    # Yellow
            self.color_accent_04 = "#F3DEAD"    # Yellow light
            #
            self.window_custom_mineralogy = tk.Toplevel(self.parent_rock)
            self.window_custom_mineralogy.title("GebPy")
            self.window_custom_mineralogy.geometry("1400x960")
            self.window_custom_mineralogy.resizable(False, False)
            self.window_custom_mineralogy["bg"] = self.color_menu
            #

            self.var_custom_mineralogy = {}
            self.var_custom_mineralogy["checkbox"] = {}
            self.gui_custom_rock = {}
            list_mineral_classes = ["Oxides", "Carbonates", "Sulfates", "Sulfides", "Tectosilicates", "Phyllosilicates",
                                    "Nesosilicates", "Inosilicates", "Sorosilicates", "Organics"]
            for mineral_class in list_mineral_classes:
                self.gui_custom_rock[mineral_class] = {}
                self.gui_custom_rock[mineral_class]["start line"] = 0
                self.gui_custom_rock[mineral_class]["end line"] = 0
                self.gui_custom_rock[mineral_class]["labels"] = {}
                self.gui_custom_rock[mineral_class]["entries"] = {}
                self.gui_custom_rock[mineral_class]["checkbuttons"] = {}
            #
            ## LABELS
            # Oxides
            start_oxides = 0
            self.gui_custom_rock["Oxides"]["start line"] = start_oxides
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_oxides, column_id=0, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Oxides", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+1, column_id=0, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+1, column_id=1, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+1, column_id=2, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+1, column_id=3, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Oxides"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_oxides = ["Quartz", "Magnetite", "Hematite", "Corundum", "Ilmenite", "Rutile", "Pyrolusite",
                            "Cassiterite", "Coltan", "Chromite"]
            list_oxides.sort()
            for index, oxide in enumerate(list_oxides, start=2):
                self.var_custom_mineralogy["checkbox"][oxide] = tk.IntVar()
                self.entr_w["custom"][oxide] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+index,
                         column_id=0, bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=oxide, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+index,
                        column_id=1, bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][oxide],
                                                                                                     command=lambda var_cb=self.var_custom_mineralogy["checkbox"][oxide], name=oxide: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+index,
                          column_id=2, bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][oxide][0],
                                                                                                    var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_oxides+index,
                          column_id=3, bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][oxide][1],
                                                                                                    var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Oxides"]["labels"][oxide] = lbl
                self.gui_custom_rock[oxide] = {}
                self.gui_custom_rock[oxide]["checkbuttons"] = cb
                self.gui_custom_rock[oxide]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Oxides"]["end line"] = start_oxides + index
                #
                self.gui_custom_rock[oxide]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][oxide][0], name=oxide, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[oxide]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][oxide][1], name=oxide, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Tectosilicates
            start_tectosilicates = self.gui_custom_rock["Oxides"]["end line"] + 1
            self.gui_custom_rock["Tectosilicates"]["start line"] = start_tectosilicates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates, column_id=0, n_columns=4,
                                    bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Tectosilicates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+1, column_id=0, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+1, column_id=1, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+1, column_id=2, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+1, column_id=3, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Tectosilicates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min,
                                                                          lbl_max]
            list_tectosilicates = ["Alkali Feldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline"]
            list_tectosilicates.sort()
            for index, tectosilicate in enumerate(list_tectosilicates, start=2):
                self.var_custom_mineralogy["checkbox"][tectosilicate] = tk.IntVar()
                self.entr_w["custom"][tectosilicate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+index, column_id=0, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=tectosilicate, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+index, column_id=1, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][tectosilicate],
                                                                command=lambda var_cb=self.var_custom_mineralogy["checkbox"][tectosilicate], name=tectosilicate: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+index, column_id=2, bg=self.color_accent_02,
                          fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][tectosilicate][0],
                                                              var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_tectosilicates+index, column_id=3, bg=self.color_accent_02,
                          fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][tectosilicate][1],
                                                              var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Tectosilicates"]["labels"][tectosilicate] = lbl
                self.gui_custom_rock[tectosilicate] = {}
                self.gui_custom_rock[tectosilicate]["checkbuttons"] = cb
                self.gui_custom_rock[tectosilicate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Tectosilicates"]["end line"] = start_tectosilicates + index
                #
                self.gui_custom_rock[tectosilicate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][tectosilicate][0], name=tectosilicate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[tectosilicate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][tectosilicate][1], name=tectosilicate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Phyllosilicates
            start_phyllosilicates = self.gui_custom_rock["Tectosilicates"]["end line"] + 1
            self.gui_custom_rock["Phyllosilicates"]["start line"] = start_phyllosilicates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates, column_id=0, n_columns=4,
                           bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Phyllosilicates",
                                                                                        relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+1, column_id=0,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+1, column_id=1,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+1, column_id=2,
                         bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+1, column_id=3,
                         bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Phyllosilicates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min,
                                                                           lbl_max]
            list_phyllosilicates = ["Illite", "Kaolinite", "Montmorillonite", "Chlorite", "Vermiculite", "Biotite",
                                    "Muscovite", "Glauconite", "Nontronite", "Saponite"]
            list_phyllosilicates.sort()
            for index, phyllosilicate in enumerate(list_phyllosilicates, start=2):
                self.var_custom_mineralogy["checkbox"][phyllosilicate] = tk.IntVar()
                self.entr_w["custom"][phyllosilicate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+index, column_id=0,
                         bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=phyllosilicate, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+index, column_id=1,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][phyllosilicate],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][phyllosilicate], name=phyllosilicate: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+index, column_id=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][phyllosilicate][0],
                                                                                           var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_phyllosilicates+index, column_id=3,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][phyllosilicate][1],
                                                                                           var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Phyllosilicates"]["labels"][phyllosilicate] = lbl
                self.gui_custom_rock[phyllosilicate] = {}
                self.gui_custom_rock[phyllosilicate]["checkbuttons"] = cb
                self.gui_custom_rock[phyllosilicate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Phyllosilicates"]["end line"] = start_phyllosilicates + index
                #
                self.gui_custom_rock[phyllosilicate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][phyllosilicate][0], name=phyllosilicate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[phyllosilicate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][phyllosilicate][1], name=phyllosilicate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Nesosilicates
            start_nesosilicates = self.gui_custom_rock["Phyllosilicates"]["end line"] + 1
            self.gui_custom_rock["Nesosilicates"]["start line"] = start_nesosilicates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates, column_id=0, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Nesosilicates",
                                                                                         relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+1, column_id=0,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+1, column_id=1,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+1, column_id=2,
                         bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+1, column_id=3,
                         bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Nesosilicates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_nesosilicates = ["Olivine", "Ca-Garnet", "Al-Garnet", "Zircon", "Thorite", "Andalusite", "Kyanite", "Sillimanite",
                                  "Topaz", "Staurolite"]
            list_nesosilicates.sort()
            for index, nesosilicate in enumerate(list_nesosilicates, start=2):
                self.var_custom_mineralogy["checkbox"][nesosilicate] = tk.IntVar()
                self.entr_w["custom"][nesosilicate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+index, column_id=0,
                         bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=nesosilicate,
                                                                                      relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+index, column_id=1,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][nesosilicate],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][nesosilicate], name=nesosilicate: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+index, column_id=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][nesosilicate][0],
                                                                                           var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_nesosilicates+index, column_id=3,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][nesosilicate][1],
                                                                                           var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Nesosilicates"]["labels"][nesosilicate] = lbl
                self.gui_custom_rock[nesosilicate] = {}
                self.gui_custom_rock[nesosilicate]["checkbuttons"] = cb
                self.gui_custom_rock[nesosilicate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Nesosilicates"]["end line"] = start_nesosilicates + index
                #
                self.gui_custom_rock[nesosilicate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][nesosilicate][0], name=nesosilicate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[nesosilicate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][nesosilicate][1], name=nesosilicate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Inosilicates
            start_inosilicates = start_oxides
            self.gui_custom_rock["Inosilicates"]["start line"] = start_inosilicates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates, column_id=4, n_columns=4,
                           bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Inosilicates",
                                                                                        relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+1, column_id=4, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+1, column_id=5, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+1, column_id=6, bg=self.color_accent_01,
                         fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+1, column_id=7, bg=self.color_accent_01,
                         fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Inosilicates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_inosilicates = ["Enstatite", "Diopside", "Augite", "Jadeite", "Aegirine", "Actinolite", "Arfvedsonite",
                                 "Ca-Amphibole", "Ca-Pyroxene", "Ferrosilite", "Glaucophane", "Mg-Fe-Pyroxene",
                                 "Riebeckite", "Na-Amphibole", "Spodumene", "Tremolite", "Wollastonite", "Donpeacorite",
                                 "Orthopyroxene"]
            list_inosilicates.sort()
            for index, inosilicate in enumerate(list_inosilicates, start=2):
                self.var_custom_mineralogy["checkbox"][inosilicate] = tk.IntVar()
                self.entr_w["custom"][inosilicate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+index, column_id=4, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=inosilicate, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+index, column_id=5, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][inosilicate],
                                                                command=lambda var_cb=self.var_custom_mineralogy["checkbox"][inosilicate], name=inosilicate: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+index, column_id=6, bg=self.color_accent_02,
                          fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][inosilicate][0],
                                                              var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_inosilicates+index, column_id=7, bg=self.color_accent_02,
                          fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][inosilicate][1],
                                                              var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Inosilicates"]["labels"][inosilicate] = lbl
                self.gui_custom_rock[inosilicate] = {}
                self.gui_custom_rock[inosilicate]["checkbuttons"] = cb
                self.gui_custom_rock[inosilicate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Inosilicates"]["end line"] = start_inosilicates + index
                #
                self.gui_custom_rock[inosilicate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][inosilicate][0], name=inosilicate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[inosilicate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][inosilicate][1], name=inosilicate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Sorosilicates
            start_sorosilicates = self.gui_custom_rock["Inosilicates"]["end line"] + 1
            self.gui_custom_rock["Sorosilicates"]["start line"] = start_sorosilicates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates, column_id=4, n_columns=4,
                           bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Sorosilicates",
                                                                                        relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+1, column_id=4, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+1, column_id=5, bg=self.color_accent_01,
                          fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+1, column_id=6, bg=self.color_accent_01,
                         fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+1, column_id=7, bg=self.color_accent_01,
                         fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Sorosilicates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_sorosilicates = ["Epidote", "Zoisite", "Gehlenite"]
            list_sorosilicates.sort()
            for index, sorosilicate in enumerate(list_sorosilicates, start=2):
                self.var_custom_mineralogy["checkbox"][sorosilicate] = tk.IntVar()
                self.entr_w["custom"][sorosilicate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+index, column_id=4, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=sorosilicate, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+index, column_id=5,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][sorosilicate],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][sorosilicate], name=sorosilicate: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+index, column_id=6,
                          bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][sorosilicate][0],
                                                                                       var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_sorosilicates+index, column_id=7,
                          bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][sorosilicate][1],
                                                                                       var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Sorosilicates"]["labels"][sorosilicate] = lbl
                self.gui_custom_rock[sorosilicate] = {}
                self.gui_custom_rock[sorosilicate]["checkbuttons"] = cb
                self.gui_custom_rock[sorosilicate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Sorosilicates"]["end line"] = start_sorosilicates + index
                #
                self.gui_custom_rock[sorosilicate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][sorosilicate][0], name=sorosilicate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[sorosilicate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][sorosilicate][1], name=sorosilicate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Carbonates
            start_carbonates = self.gui_custom_rock["Sorosilicates"]["end line"] + 1
            self.gui_custom_rock["Carbonates"]["start line"] = start_carbonates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates, column_id=4, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Carbonates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+1, column_id=4, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+1, column_id=5, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+1, column_id=6, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+1, column_id=7, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Carbonates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_carbonates = ["Calcite", "Dolomite", "Magnesite", "Rhodochrosite", "Siderite", "Aragonite",
                               "Cerussite", "Ankerite", "Azurite", "Malachite"]
            list_carbonates.sort()
            for index, carbonate in enumerate(list_carbonates, start=2):
                self.var_custom_mineralogy["checkbox"][carbonate] = tk.IntVar()
                self.entr_w["custom"][carbonate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+index, column_id=4,
                         bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=carbonate, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+index, column_id=5,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][carbonate],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][carbonate], name=carbonate: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+index, column_id=6,
                          bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][carbonate][0],
                                                                                       var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_carbonates+index, column_id=7,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][carbonate][1],
                                                                                           var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Carbonates"]["labels"][carbonate] = lbl
                self.gui_custom_rock[carbonate] = {}
                self.gui_custom_rock[carbonate]["checkbuttons"] = cb
                self.gui_custom_rock[carbonate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Carbonates"]["end line"] = start_carbonates + index
                #
                self.gui_custom_rock[carbonate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][carbonate][0], name=carbonate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[carbonate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][carbonate][1], name=carbonate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Sulfides
            start_sulfides = start_oxides
            self.gui_custom_rock["Sulfides"]["start line"] = start_sulfides
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides, column_id=8, n_columns=4,
                           bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Sulfides", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+1, column_id=8,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+1, column_id=9,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+1, column_id=10,
                         bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+1, column_id=11,
                         bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Sulfides"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_sulfides = ["Pyrite", "Chalcopyrite", "Bornite", "Covellite", "Molybdenite", "Sphalerite", "Galena",
                             "Fahlore"]
            list_sulfides.sort()
            for index, sulfide in enumerate(list_sulfides, start=2):
                self.var_custom_mineralogy["checkbox"][sulfide] = tk.IntVar()
                self.entr_w["custom"][sulfide] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+index, column_id=8,
                         bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=sulfide, relief=tk.RAISED)
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+index, column_id=9,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][sulfide],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][sulfide], name=sulfide: self.marked_checkbox(var_cb, name))
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+index, column_id=10,
                          bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][sulfide][0],
                                                                                       var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_sulfides+index, column_id=11,
                          bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][sulfide][1],
                                                                                       var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Sulfides"]["labels"][sulfide] = lbl
                self.gui_custom_rock[sulfide] = {}
                self.gui_custom_rock[sulfide]["checkbuttons"] = cb
                self.gui_custom_rock[sulfide]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Sulfides"]["end line"] = start_sulfides + index
                #
                self.gui_custom_rock[sulfide]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][sulfide][0], name=sulfide, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[sulfide]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][sulfide][1], name=sulfide, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Sulfates
            start_sulfates = self.gui_custom_rock["Sulfides"]["end line"] + 1
            self.gui_custom_rock["Sulfates"]["start line"] = start_sulfates
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates, column_id=8, n_columns=4,
                              bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Sulfates",
                                                                                           relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+1, column_id=8,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+1, column_id=9,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+1, column_id=10,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+1, column_id=11,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Sulfates"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_sulfates = ["Barite", "Anhydrite", "Gypsum", "Alunite", "Jarosite", "Anglesite", "Hanksite",
                             "Celestine", "Kieserite", "Chalcanthite", "Scheelite", "Kainite"]
            list_sulfates.sort()
            for index, sulfate in enumerate(list_sulfates, start=2):
                self.var_custom_mineralogy["checkbox"][sulfate] = tk.IntVar()
                self.entr_w["custom"][sulfate] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+index, column_id=8,
                         bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=sulfate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][sulfate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+index, column_id=9,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][sulfate],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][sulfate], name=sulfate: self.marked_checkbox(var_cb, name))
                self.entr_w["custom"][sulfate] = [tk.StringVar(), tk.StringVar()]
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+index, column_id=10,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][sulfate][0],
                                                                                           var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_sulfates+index, column_id=11,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][sulfate][1],
                                                                                           var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Sulfates"]["labels"][sulfate] = lbl
                self.gui_custom_rock[sulfate] = {}
                self.gui_custom_rock[sulfate]["checkbuttons"] = cb
                self.gui_custom_rock[sulfate]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Sulfates"]["end line"] = start_sulfates + index
                #
                self.gui_custom_rock[sulfate]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][sulfate][0], name=sulfate, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[sulfate]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][sulfate][1], name=sulfate, pos=1: self.set_entry(var_entr, name, pos, event))
            #
            # Organics
            start_organics = self.gui_custom_rock["Sulfates"]["end line"] + 1
            self.gui_custom_rock["Organics"]["start line"] = start_organics
            lbl_title = SE(parent=self.window_custom_mineralogy, row_id=start_organics, column_id=8, n_columns=4,
                              bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Organics",
                                                                                           relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=start_organics+1, column_id=8,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_part = SE(parent=self.window_custom_mineralogy, row_id=start_organics+1, column_id=9,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_min = SE(parent=self.window_custom_mineralogy, row_id=start_organics+1, column_id=10,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_max = SE(parent=self.window_custom_mineralogy, row_id=start_organics+1, column_id=11,
                          bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            self.gui_custom_rock["Organics"]["labels"]["header"] = [lbl_title, lbl_name, lbl_part, lbl_min, lbl_max]
            list_organics = ["Organic Matter"]
            list_organics.sort()
            for index, organic in enumerate(list_organics, start=2):
                self.var_custom_mineralogy["checkbox"][organic] = tk.IntVar()
                self.entr_w["custom"][organic] = [tk.StringVar(), tk.StringVar()]
                #
                lbl = SE(parent=self.window_custom_mineralogy, row_id=start_organics+index, column_id=8,
                         bg=self.color_accent_02, fg=self.color_fg_dark).create_label(text=organic, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][organic] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=start_organics+index, column_id=9,
                        bg=self.color_accent_02, fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][organic],
                                                                                        command=lambda var_cb=self.var_custom_mineralogy["checkbox"][organic], name=organic: self.marked_checkbox(var_cb, name))
                self.entr_w["custom"][organic] = [tk.StringVar(), tk.StringVar()]
                entr_min = SE(parent=self.window_custom_mineralogy, row_id=start_organics+index, column_id=10,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][organic][0],
                                                                                           var_entr_set=0.0, width=5)
                entr_max = SE(parent=self.window_custom_mineralogy, row_id=start_organics+index, column_id=11,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_entry(var_entr=self.entr_w["custom"][organic][1],
                                                                                           var_entr_set=1.0, width=5)
                #
                self.gui_custom_rock["Organics"]["labels"][organic] = lbl
                self.gui_custom_rock[organic] = {}
                self.gui_custom_rock[organic]["checkbuttons"] = cb
                self.gui_custom_rock[organic]["entries"] = [entr_min, entr_max]
                self.gui_custom_rock["Organics"]["end line"] = start_organics + index
                #
                self.gui_custom_rock[organic]["entries"][0].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][organic][0], name=organic, pos=0: self.set_entry(var_entr, name, pos, event))
                self.gui_custom_rock[organic]["entries"][1].bind("<Return>", lambda event, var_entr=self.entr_w["custom"][organic][1], name=organic, pos=1: self.set_entry(var_entr, name, pos, event))
        #
        elif var_btn == "Generate Data":
            try:
                self.list_elements.clear()
                self.list_minerals.clear()
            except:
                pass
            #
            self.var_rb.set(0)
            self.generate_custom_rock_data()
    #
    def marked_checkbox(self, var_cb, name):
        if var_cb.get() == 1:
            try:
                self.custom_mineralogy["mineralogy"][name] = [self.gui_custom_rock["Sulfates"]["entries"][name][0].get(),
                                                              self.gui_custom_rock["Sulfates"]["entries"][name][1].get()]
            except:
                pass
            try:
                self.custom_mineralogy["mineralogy"][name] = [self.gui_custom_rock[name]["entries"][0].get(),
                                                              self.gui_custom_rock[name]["entries"][1].get()]
            except:
                pass
        else:
            del self.custom_mineralogy["mineralogy"][name]
    #
    def set_entry(self, var_entr, name, pos, event):
        var_entr.set(var_entr.get())
        cb_state = self.var_custom_mineralogy["checkbox"][name].get()
        if cb_state == 1:
            self.custom_mineralogy["mineralogy"][name][pos] = var_entr.get()
        else:
            self.custom_mineralogy["mineralogy"][name] = [self.gui_custom_rock[name]["entries"][0].get(),
                                                          self.gui_custom_rock[name]["entries"][1].get()]
            self.var_custom_mineralogy["checkbox"][name].set(1)
    #
    def generate_custom_rock_data(self):
        #
        try:
            for lbl in self.lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in self.entr_w["chemistry"]:
                entr.grid_forget()
            self.lbl_w["chemistry"].clear()
            self.entr_w["chemistry"].clear()
        except:
            pass
        #
        try:
            self.fig.clf()
            self.ax.cla()
            self.canvas.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas:
                self.canvas.destroy()
        except AttributeError:
            pass
        #
        try:
            plt.close("all")
        except:
            pass
        #
        water = fluids.Water.water("")
        n_samples = self.var_entr.get()
        phi_min = round(self.var_phi0.get()/100, 4)
        phi_max = round(self.var_phi1.get()/100, 4)
        self.phi = [round(rd.uniform(phi_min, phi_max), 4) for i in range(n_samples)]
        assemblage = list(self.custom_mineralogy["mineralogy"].keys())
        self.minerals = {}
        self.elements = {}
        elements_list = []
        mineral_list = []
        data_minerals = {}
        self.mineralogy_custom = []
        self.chemistry_custom = []
        #
        for mineral in assemblage:
            # Tectosilicates
            if mineral in ["Alkali Feldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline"]:
                data = Tectosilicates(data_type=True, mineral=mineral).get_data()
            # Phyllosilicates
            elif mineral in ["Illite", "Kaolinite", "Montmorillonite", "Chlorite", "Vermiculite", "Biotite",
                             "Muscovite", "Glauconite"]:
                data = Phyllosilicates(data_type=True, mineral=mineral).get_data()
            # Nesosilicates
            elif mineral in ["Zircon", "Thorite", "Topaz", "Pyrope", "Olivine", "Ca-Garnet", "Al-Garnet", "Uvarovite",
                             "Andratite", "Grossular", "Almandine", "Liebenbergite", "Ca-Olivine", "Tephroite",
                             "Fayalite", "Forsterite", "Staurolite", "Sillimanite", "Kyanite", "Andalusite",
                             "Spessartine"]:
                data = Nesosilicates(data_type=True, mineral=mineral).get_data()
            # Sorosilicates
            elif mineral in ["Epidote", "Zoisite", "Gehlenite"]:
                data = Sorosilicates(data_type=True, mineral=mineral).get_data()
            # Inosilicates
            elif mineral in ["Enstatite", "Diopside", "Augite", "Jadeite", "Aegirine", "Actinolite", "Arfvedsonite",
                             "Ca-Amphibole", "Ca-Pyroxene", "Ferrosilite", "Glaucophane", "Mg-Fe-Pyroxene",
                             "Riebeckite", "Na-Amphibole", "Spodumene", "Tremolite", "Wollastonite", "Donpeacorite",
                             "Orthopyroxene"]:
                data = Inosilicates(data_type=True, mineral=mineral).get_data()
            # Carbonates
            elif mineral in ["Calcite", "Dolomite", "Magnesite", "Rhodochrosite", "Siderite", "Aragonite", "Cerussite",
                             "Ankerite", "Azurite", "Malachite"]:
                data = Carbonates(data_type=True, mineral=mineral).get_data()
            # Sulfides
            elif mineral in ["Pyrite", "Chalcopyrite", "Bornite", "Covellite", "Molybdenite", "Sphalerite", "Galena",
                             "Fahlore"]:
                data = Sulfides(data_type=True, mineral=mineral).get_data()
            # Sulfates
            elif mineral in ["Alunite", "Anglesite", "Anhydrite", "Barite", "Celestine", "Chalcanthite", "Gypsum",
                             "Hanksite", "Jarosite", "Kieserite", "Scheelite", "Kainite"]:
                data = Sulfates(data_type=True, mineral=mineral).get_data()
            # Organics
            elif mineral in ["Organic Matter"]:
                data = Organics(data_type=True, mineral=mineral).get_data()
            # Oxides
            else:
                data = Oxides(data_type=True, mineral=mineral).get_data()
            self.minerals[data["mineral"]] = []
            mineral_list.append(data["mineral"])
            elements_list.extend(list(data["chemistry"].keys()))
            if data["state"] == "variable":
                # Tectosilicates
                if mineral in ["Alkali Feldspar", "Plagioclase", "Scapolite", "Danburite", "Nepheline"]:
                    dataset = Tectosilicates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Phyllosilicates
                elif mineral in ["Illite", "Kaolinite", "Montmorillonite", "Chlorite", "Vermiculite", "Biotite",
                                 "Muscovite", "Glauconite"]:
                    dataset = Phyllosilicates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Nesosilicates
                elif mineral in ["Zircon", "Thorite", "Topaz", "Pyrope", "Olivine", "Ca-Garnet", "Al-Garnet", "Uvarovite",
                                 "Andratite", "Grossular", "Almandine", "Liebenbergite", "Ca-Olivine", "Tephroite",
                                 "Fayalite", "Forsterite", "Staurolite", "Sillimanite", "Kyanite", "Andalusite",
                                 "Spessartine"]:
                    dataset = Nesosilicates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Sorosilicates
                elif mineral in ["Epidote", "Zoisite", "Gehlenite"]:
                    dataset = Sorosilicates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Inosilicates
                elif mineral in ["Enstatite", "Diopside", "Augite", "Jadeite", "Aegirine", "Actinolite", "Arfvedsonite",
                                 "Ca-Amphibole", "Ca-Pyroxene", "Ferrosilite", "Glaucophane", "Mg-Fe-Pyroxene",
                                 "Riebeckite", "Na-Amphibole", "Spodumene", "Tremolite", "Wollastonite", "Donpeacorite",
                                 "Orthopyroxene"]:
                    dataset = Inosilicates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Carbonates
                elif mineral in ["Calcite", "Dolomite", "Magnesite", "Rhodochrosite", "Siderite", "Aragonite",
                                 "Cerussite", "Ankerite", "Azurite", "Malachite"]:
                    dataset = Carbonates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Sulfides
                elif mineral in ["Pyrite", "Chalcopyrite", "Bornite", "Covellite", "Molybdenite", "Sphalerite",
                                 "Galena", "Fahlore"]:
                    dataset = Sulfides(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Sulfates
                elif mineral in ["Alunite", "Anglesite", "Anhydrite", "Barite", "Celestine", "Chalcanthite", "Gypsum",
                                 "Hanksite", "Jarosite", "Kieserite", "Scheelite", "Kainite"]:
                    dataset = Sulfates(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Organics
                elif mineral in ["Organic Matter"]:
                    dataset = Organics(data_type=True, mineral=mineral).get_data(number=n_samples)
                # Oxides
                else:
                    dataset = Oxides(data_type=True, mineral=mineral).get_data(number=n_samples)
                data_minerals[data["mineral"]] = dataset
            else:
                data_minerals[data["mineral"]] = data
        #
        elements_list = list(dict.fromkeys(elements_list))
        elements_list.sort()
        mineral_chemistry = {}
        for mineral in mineral_list:
            try:
                mineral_chemistry[mineral] = list(data_minerals[mineral]["chemistry"].keys())
            except:
                mineral_chemistry[mineral] = list(data_minerals[mineral][0]["chemistry"].keys())
        for element in elements_list:
            self.elements[element] = []
        n = 0
        self.mineral_helper = {}
        self.element_helper = {}
        for mineral in mineral_list:
            self.mineral_helper[mineral] = []
        for element in elements_list:
            self.element_helper[element] = []
        index = 0
        while len(self.element_helper[elements_list[0]]) < n_samples:
            condition_01 = False
            condition_02 = False
            while condition_01 == False and condition_02 == False:
                amount_total = 0
                amount_helper = {}
                for j, mineral in enumerate(assemblage):
                    amount_helper[mineral_list[j]] = 0
                    if j < len(assemblage)-1:
                        amount = round(rd.uniform(float(self.custom_mineralogy["mineralogy"][mineral][0]),
                                                  float(self.custom_mineralogy["mineralogy"][mineral][1])), 4)
                        amount_total += amount
                        amount_helper[mineral_list[j]] = amount
                    else:
                        amount = round(1 - amount_total, 4)
                        if amount >= 0:
                            amount_helper[mineral_list[j]] = amount
                        else:
                            continue
                w_minerals = list(amount_helper.values())
                if sum(w_minerals) == 1 and amount <= float(self.custom_mineralogy["mineralogy"][mineral][1]):
                    condition_01 = True
                    for mineral, amount in amount_helper.items():
                        self.mineral_helper[mineral].append(amount)
                    amount_total = 0
                    amount_helper = {}
                    for k, element in enumerate(elements_list):
                        amount_helper[element] = 0
                        for mineral in mineral_list:
                            if element in mineral_chemistry[mineral] and k < len(elements_list)-1:
                                try:
                                    value = round(self.mineral_helper[mineral][-1]*data_minerals[mineral]["chemistry"][element], 4)
                                except:
                                    value = round(self.mineral_helper[mineral][-1]*data_minerals[mineral][index]["chemistry"][element], 4)
                                amount_helper[element] += value
                                amount_total += value
                            elif element not in mineral_chemistry[mineral]:
                                amount_helper[element] += 0
                            else:
                                amount_helper[element] = round(1 - amount_total, 4)
                    w_elements = list(amount_helper.values())
                    if sum(w_elements) == 1:
                        if all(w >= 0 for w in w_elements) == True:
                            for element, amount in amount_helper.items():
                                self.element_helper[element].append(abs(amount))
                            index += 1
                            condition_02 = True
                        else:
                            continue
                    else:
                        if len(self.mineral_helper[mineral_list[0]]) == len(self.element_helper[elements_list[0]]):
                            pass
                        elif len(self.mineral_helper[mineral_list[0]]) > len(self.element_helper[elements_list[0]]):
                            for mineral, amount in self.mineral_helper.items():
                                self.mineral_helper[mineral].pop()
                        elif len(self.mineral_helper[mineral_list[0]]) < len(self.element_helper[elements_list[0]]):
                            for element, amount in self.element_helper.items():
                                self.element_helper[element].pop()
                        continue
                else:
                    if len(self.mineral_helper[mineral_list[0]]) == len(self.element_helper[elements_list[0]]):
                        pass
                    elif len(self.mineral_helper[mineral_list[0]]) > len(self.element_helper[elements_list[0]]):
                        for mineral, amount in self.mineral_helper.items():
                            self.mineral_helper[mineral].pop()
                    elif len(self.mineral_helper[mineral_list[0]]) < len(self.element_helper[elements_list[0]]):
                        for element, amount in self.element_helper.items():
                            self.element_helper[element].pop()
                    continue
        #
        # print(self.mineral_helper)
        # print(self.element_helper)
        # for mineral, value in self.mineral_helper.items():
        #     print(mineral, len(value))
        # for element, value in self.element_helper.items():
        #     print(element, len(value))
        #
        self.list_elements = elements_list
        self.list_minerals = mineral_list
        self.elements = {}
        self.minerals = {}
        for mineral, values in self.mineral_helper.items():
            self.minerals[mineral] = list(np.around(np.array(values)*100, 2))
        for element, values in self.element_helper.items():
            self.elements[element] = list(np.around(np.array(values)*100, 2))
        for index in range(n_samples):
            self.mineralogy_custom.append({})
            for mineral, value in self.mineral_helper.items():
                self.mineralogy_custom[index][mineral] = value[index]
            self.chemistry_custom.append({})
            for element, value in self.element_helper.items():
                self.chemistry_custom[index][element] = value[index]
        #
        self.rho = []
        self.bulk_mod = []
        self.shear_mod = []
        self.youngs_mod = []
        self.poisson = []
        self.gamma_ray = []
        self.photoelectricity = []
        self.vP = []
        self.vS = []
        self.vPvS = []
        for i in range(n_samples):
            rho_value = 0
            K_list = []
            G_list = []
            w = []
            gr_value = 0
            pe_value = 0
            poisson_value = 0
            for mineral in mineral_list:
                try:
                    rho_value += self.mineral_helper[mineral][i]*data_minerals[mineral]["rho"]
                    gr_value += self.mineral_helper[mineral][i]*data_minerals[mineral]["GR"]
                    pe_value += self.mineral_helper[mineral][i]*data_minerals[mineral]["PE"]
                    w.append(self.mineral_helper[mineral][i])
                    K_list.append(round(data_minerals[mineral]["K"]*self.mineral_helper[mineral][i], 3))
                    G_list.append(round(data_minerals[mineral]["G"]*self.mineral_helper[mineral][i], 3))
                    poisson_value += self.minerals[mineral][i]*data_minerals[mineral]["nu"]
                except:
                    rho_value += self.mineral_helper[mineral][i]*data_minerals[mineral][i]["rho"]
                    gr_value += self.mineral_helper[mineral][i]*data_minerals[mineral][i]["GR"]
                    pe_value += self.mineral_helper[mineral][i]*data_minerals[mineral][i]["PE"]
                    w.append(self.mineral_helper[mineral][i])
                    K_list.append(round(data_minerals[mineral][i]["K"]*self.mineral_helper[mineral][i], 3))
                    G_list.append(round(data_minerals[mineral][i]["G"]*self.mineral_helper[mineral][i], 3))
                    poisson_value += self.mineral_helper[mineral][i]*data_minerals[mineral][i]["nu"]
            self.gamma_ray.append(round(gr_value, 2))
            self.photoelectricity.append(round(pe_value, 2))
            K_geo = elast.calc_geometric_mean(self, w, K_list)
            G_geo = elast.calc_geometric_mean(self, w, G_list)
            self.bulk_mod.append(round(K_geo, 2))
            self.shear_mod.append(round(G_geo, 2))
            self.youngs_mod.append(round((9*K_geo*G_geo)/(3*K_geo - G_geo), 2))
            #self.poisson.append(round(poisson_value, 3))
            self.poisson.append(round((3*K_geo-2*G_geo)/(6*K_geo + 2*G_geo), 4))
            vP_solid = np.sqrt((K_geo*10**9+4/3*G_geo*10**9)/(rho_value))
            vS_solid = np.sqrt((G_geo*10**9)/(rho_value))
            rho = (1 - self.phi[i])*rho_value + self.phi[i]*water[2]/1000
            vP = (1 - self.phi[i])*vP_solid + self.phi[i]*water[4][0]
            vS = (1 - self.phi[i])*vS_solid
            self.rho.append(round(rho, 3))
            self.vP.append(round(vP, 2))
            self.vS.append(round(vS, 2))
            self.vPvS.append(round(vP/vS, 3))
        self.rho = np.array(self.rho)
        self.vP = np.array(self.vP)
        self.vS = np.array(self.vS)
        self.vPvS = np.array(self.vPvS)
        self.bulk_mod = np.array(self.bulk_mod)
        self.shear_mod = np.array(self.shear_mod)
        self.poisson = np.array(self.poisson)
        self.phi = np.array(self.phi)
        self.gamma_ray = np.array(self.gamma_ray)
        self.photoelectricity = np.array(self.photoelectricity)
        #
        for index in range(n_samples):
            self.exp_data.append({"rock": self.filename[2].get(), "phi": self.phi[index], "fluid": "water",
                                  "rho": self.rho[index], "vP": self.vP[index], "vS": self.vS[index],
                                  "vP/vS": self.vPvS[index], "K": self.bulk_mod[index], "G": self.shear_mod[index],
                                  "E": self.youngs_mod[index], "nu": self.poisson[index], "GR": self.gamma_ray[index],
                                  "PE": self.photoelectricity[index], "mineralogy": self.mineralogy_custom[index],
                                  "chemistry": self.chemistry_custom[index]})
        #
        self.results = [self.rho, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod, self.poisson, self.phi,
                        self.gamma_ray, self.photoelectricity]
        #
        self.var_opt_chem = tk.StringVar()
        try:
            if len(self.opt_list_chem) > 0:
                self.opt_list_chem.clear()
        except:
            pass
        self.opt_list_chem = ["No Selection"]
        self.opt_list_chem.extend(self.list_elements)
        self.opt_list_chem.extend(self.list_minerals)
        self.opt_chem = SE(parent=self.parent_rock, row_id=36, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                           fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element/Mineral",
                                                          opt_list=self.opt_list_chem, active_bg=self.color_acc_02,
                                                          command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
        self.gui_elements.append(self.opt_chem)
        #
        lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                 fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        self.lbl_w["chemistry"].append(lbl)
        #
        self.entr_list_min = []
        self.entr_list_max = []
        self.entr_list_mean = []
        self.entr_list_std = []
        for i in range(10+len(self.list_elements)):
            self.entr_list_min.append(tk.IntVar())
            self.entr_list_max.append(tk.IntVar())
            self.entr_list_mean.append(tk.IntVar())
            self.entr_list_std.append(tk.IntVar())
        #
        ## Entry Table
        for i in range(10):
            if i == 7:
                entr_min = SE(parent=self.parent_rock, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]*100), 3))
                entr_max = SE(parent=self.parent_rock, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]*100), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results[i]*100), 3))
                entr_std = SE(parent=self.parent_rock, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results[i]*100, ddof=1), 3))
            else:
                entr_min = SE(parent=self.parent_rock, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results[i]), 3))
                entr_max = SE(parent=self.parent_rock, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results[i]), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results[i]), 3))
                entr_std = SE(parent=self.parent_rock, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results[i], ddof=1), 3))
        for index, element in enumerate(self.list_elements, start=10):
            if element not in ["U"]:
                entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(self.elements[element]), 3))
                entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(self.elements[element]), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(self.elements[element]), 3))
                entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(self.elements[element], ddof=1), 3))
            else:
                ppm_amounts = np.array(self.elements[element])*10000
                entr_min = SE(parent=self.parent_rock, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(ppm_amounts), 3))
                entr_max = SE(parent=self.parent_rock, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(ppm_amounts), 3))
                entr_mean = SE(parent=self.parent_rock, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(ppm_amounts), 3))
                entr_std = SE(parent=self.parent_rock, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
            #
            self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        for index, element in enumerate(self.list_elements, start=0):
            if element not in ["U"]:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element), relief=tk.RAISED)
            else:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
        #
        self.color_rock = "#7C9097"
        #
        self.labels = [["Densitiy $\\varrho$ (g/ccm)", "Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)"],
                       ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                       ["Porosity $\\phi$ (%)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
        self.labels_scatter_rho = [["Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)", "Velocity Ratio $v_P/v_S$ (1)"],
                       ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                       ["Porosity $\\phi$ (%)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
        self.labels_scatter_phi = [["Seismic velocity $v_P$ (km/s)", "Seismic velocity $v_S$ (km/s)", "Velocity Ratio $v_P/v_S$ (1)"],
                       ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)", "Poisson's ratio $\\mu$ (1)"],
                       ["Densitiy $\\varrho$ (g/ccm)", "Gamma ray GR (API)", "Photoelectricity PE (barns/electron)"]]
        self.data_plot = [[self.rho/1000, self.vP/1000, self.vS/1000],
                          [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.phi*100, self.gamma_ray, self.photoelectricity]]
        self.data_plot_scatter_rho = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                      [self.bulk_mod, self.shear_mod, self.poisson],
                                      [self.phi*100, self.gamma_ray, self.photoelectricity]]
        self.data_plot_scatter_phi = [[self.vP/1000, self.vS/1000, self.vP/self.vS],
                                      [self.bulk_mod, self.shear_mod, self.poisson],
                                      [self.rho/1000, self.gamma_ray, self.photoelectricity]]
        self.create_3x3_histo(parent=self.parent_rock, data=self.data_plot, row_id=2, column_id=9, n_rows=45,
                              n_columns=9, color=self.color_rock, labels=self.labels)
    #
    def __call__(self):
        return self.lbl_w, self.entr_w, self.exp_data, self.filename
#
class Subsurface:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, subsurface, lbl_w, entr_w, gui_elements):
        #
        try:
            for lbl in lbl_w["physics"]:
                lbl.grid_forget()
            for entr in entr_w["physics"]:
                entr.grid_forget()
            lbl_w["physics"].clear()
            entr_w["physics"].clear()
        except:
            pass
        #
        try:
            for lbl in lbl_w["chemistry"]:
                lbl.grid_forget()
            for entr in entr_w["chemistry"]:
                entr.grid_forget()
            lbl_w["chemistry"].clear()
            entr_w["chemistry"].clear()
        except:
            pass
        #
        try:
            for gui_elmnt in gui_elements:
                gui_elmnt.grid_forget()
            gui_elmnt.clear()
        except:
            pass
        #
        self.parent_subsurface = parent
        self.color_bg = color_bg
        self.color_fg = color_fg
        self.color_acc_01 = color_acc[0]
        self.color_acc_02 = color_acc[1]
        self.var_entr = tk.IntVar()
        var_entr_start = 100
        self.var_phi0 = tk.IntVar()
        var_phi0_start = 5
        self.var_phi1 = tk.IntVar()
        var_phi1_start = 30
        self.var_rb = tk.IntVar()
        self.var_rb_geochem = tk.IntVar()
        var_rb_start = 0
        var_rb_geochem_start = 2
        self.subsurface = subsurface
        self.lbl_w = lbl_w
        self.entr_w = entr_w
        self.gui_elements = gui_elements
        #
        if self.subsurface == "random":
            self.create_random_sequences(thickness=1000, style="siliciclastic")
            #
            lbl_stat = SE(parent=self.parent_subsurface, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
               fg="black").create_label(text="Statistics - Random Subsurface (siliciclastic)", relief=tk.RAISED)
            lbl_plt = SE(parent=self.parent_subsurface, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
               fg="black").create_label(text="Statistics - Random Subsurface (siliciclastic)", relief=tk.RAISED)
            #
            self.gui_elements.extend([lbl_stat, lbl_plt])
            #
        ## Labels
        lbl_01 = SE(parent=self.parent_subsurface, row_id=4, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Density\n (g/ccm)", relief=tk.RAISED)
        lbl_02 = SE(parent=self.parent_subsurface, row_id=6, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="P-wave velocity\n (km/s)", relief=tk.RAISED)
        lbl_03 = SE(parent=self.parent_subsurface, row_id=8, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="S-wave velocity\n (km/s)", relief=tk.RAISED)
        lbl_04 = SE(parent=self.parent_subsurface, row_id=10, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Velocity ratio\n (1)", relief=tk.RAISED)
        lbl_05 = SE(parent=self.parent_subsurface, row_id=12, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Bulk modulus\n (GPa)", relief=tk.RAISED)
        lbl_06 = SE(parent=self.parent_subsurface, row_id=14, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Shear modulus\n (GPa)", relief=tk.RAISED)
        lbl_07 = SE(parent=self.parent_subsurface, row_id=16, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Poisson's ratio\n (1)", relief=tk.RAISED)
        lbl_08 = SE(parent=self.parent_subsurface, row_id=18, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Porosity\n (%)", relief=tk.RAISED)
        lbl_09 = SE(parent=self.parent_subsurface, row_id=20, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Gamma ray\n (API)", relief=tk.RAISED)
        lbl_10 = SE(parent=self.parent_subsurface, row_id=22, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Photoelectricity\n (barns/electron)", relief=tk.RAISED)
        #
        self.lbl_w["physics"].extend([lbl_01, lbl_02, lbl_03, lbl_04, lbl_05, lbl_06, lbl_07, lbl_08, lbl_09, lbl_10])
    #
    def create_real_world_sequences(self, name, data_units):
        self.name = name
        self.data_rocks = {}
        #
        lbl_stat = SE(parent=self.parent_subsurface, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
                      fg="black").create_label(text="Statistics - "+str(name), relief=tk.RAISED)
        lbl_plt = SE(parent=self.parent_subsurface, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
                     fg="black").create_label(text="Statistics - "+str(name), relief=tk.RAISED)
        #
        self.gui_elements.extend([lbl_stat, lbl_plt])
        #
        self.results_sorted = {}
        self.results_units = {}
        self.list_rocks_short = []
        properties = ["thickness", "rock", "rho", "vP", "vS", "vPvS", "phi", "K", "G", "Poisson", "GR", "PE", "Top",
                      "Bottom"]
        for prop in properties:
            self.results_sorted[prop] = []
        #
        #self.results_sorted["Bottom"].append(100.0)
        for unit in data_units:
            for key, value in unit.items():
                rock = value["rock"]
                break
            if rock not in self.list_rocks_short:
                self.list_rocks_short.append(rock)
        for rock in self.list_rocks_short:
            self.results_units[rock] = {}
            for property in properties:
                self.results_units[rock][property] = []
        n_units = []
        self.results_sorted["Top"].append(0.0)
        #
        self.data_chemistry = {}
        self.data_mineralogy = {}
        for rock in self.list_rocks_short:
            self.data_chemistry[rock] = {}
            self.data_mineralogy[rock] = {}
        #
        i = 0
        for unit in data_units:
            n_unit = len(unit)
            n_units.append(n_unit)
            for key, value in unit.items():
                if i == 0:
                    self.results_sorted["Bottom"].append(key)
                else:
                    self.results_sorted["Top"].append(key)
                    self.results_sorted["Bottom"].append(key)
                try:
                    self.results_sorted["rho"].append(value["rho"][0])
                    self.results_sorted["vP"].append(value["vP"][0])
                    self.results_sorted["vS"].append(value["vS"][0])
                    self.results_sorted["vPvS"].append(value["vP/vS"][0])
                    self.results_sorted["phi"].append(value["phi"][0])
                    self.results_sorted["K"].append(value["K"][0])
                    self.results_sorted["G"].append(value["G"][0])
                    self.results_sorted["Poisson"].append(value["nu"][0])
                    self.results_sorted["GR"].append(value["GR"][0])
                    self.results_sorted["PE"].append(value["PE"][0])
                    self.results_sorted["rock"].append(value["rock"])
                except:
                    self.results_sorted["rho"].append(value["rho"])
                    self.results_sorted["vP"].append(value["vP"])
                    self.results_sorted["vS"].append(value["vS"])
                    self.results_sorted["vPvS"].append(value["vP/vS"])
                    self.results_sorted["phi"].append(value["phi"])
                    self.results_sorted["K"].append(value["K"])
                    self.results_sorted["G"].append(value["G"])
                    self.results_sorted["Poisson"].append(value["nu"])
                    self.results_sorted["GR"].append(value["GR"])
                    self.results_sorted["PE"].append(value["PE"])
                    self.results_sorted["rock"].append(value["rock"])
                #
                if len(self.data_chemistry[value["rock"]]) == 0:
                    for element in value["chemistry"].keys():
                       self.data_chemistry[value["rock"]][element] = []
                #
                for element, value_element in value["chemistry"].items():
                    self.data_chemistry[value["rock"]][element].append(value_element)
                #
                if len(self.data_mineralogy[value["rock"]]) == 0:
                    for mineral in value["mineralogy"].keys():
                       self.data_mineralogy[value["rock"]][mineral] = []
                #
                for mineral, value_mineral in value["mineralogy"].items():
                    self.data_mineralogy[value["rock"]][mineral].append(value_mineral)
                #
                i += 1
        #
        self.results_sorted["Bottom"][-1]
        np.sort(np.array(self.results_sorted["Bottom"]))
        np.sort(np.array(self.results_sorted["Top"]))
        self.results_sorted["thickness"].append(self.results_sorted["Bottom"][0])
        for index, item in enumerate(self.results_sorted["Bottom"][1:]):
            self.results_sorted["thickness"].append(item - self.results_sorted["thickness"][-1])
        a = np.around([j - i for i, j in zip(self.results_sorted["Bottom"][:-1], self.results_sorted["Bottom"][1:])], 4)
        a = np.insert(a, 0, a[0])
        self.results_sorted["thickness"] = list(a)
        #
        start = 0
        end = n_units[0]
        for index, rock in enumerate(self.list_rocks_short):
            for property in properties:
                self.results_units[rock][property].extend(self.results_sorted[property][start:end])
            start += n_units[index]
            end += n_units[index]
        #
        # for index, value in enumerate(self.results_sorted["Top"]):
        #     print(value, self.results_sorted["Bottom"][index])
        #
        # print("Rocks:", self.list_rocks_short)
        # print("Rocks:", self.results_sorted["rock"], len(self.results_sorted["rock"]))
        # print("Top:", self.results_sorted["Top"], len(self.results_sorted["Top"]))
        # print("Bottom:", self.results_sorted["Bottom"], len(self.results_sorted["Bottom"]))
        # print("Thickness:", self.results_sorted["thickness"], len(self.results_sorted["thickness"]))
        #
        self.unit_sections = {}
        for rock in self.list_rocks_short:
            if rock == "Sandstone":
                self.unit_sections[rock] = {"Intervals": [], "Color": "tan"}
            elif rock in ["Shale", "Mudstone"]:
                self.unit_sections[rock] = {"Intervals": [], "Color": "olivedrab"}
            elif rock in ["Granite", "Gabbro", "Diorite"]:
                self.unit_sections[rock] = {"Intervals": [], "Color": "darkorange"}
            elif rock == "Kupferschiefer":
                self.unit_sections[rock] = {"Intervals": [], "Color": "gray"}
            elif rock in ["limestone", "Limestone"]:
                self.unit_sections[rock] = {"Intervals": [], "Color": "skyblue"}
            elif rock == "Anhydrite":
                self.unit_sections[rock] = {"Intervals": [], "Color": "orchid"}
            elif rock == "Dolomite":
                self.unit_sections[rock] = {"Intervals": [], "Color": "lightcyan"}
            elif rock == "Rock Salt":
                self.unit_sections[rock] = {"Intervals": [], "Color": "lavender"}
            elif rock == "Potash":
                self.unit_sections[rock] = {"Intervals": [], "Color": "yellowgreen"}
            elif rock == "Marl":
                self.unit_sections[rock] = {"Intervals": [], "Color": "moccasin"}
        for index, rock in enumerate(self.results_sorted["rock"]):
            if index > 1:
                top = round(self.results_sorted["Bottom"][index - 1], 4)
                bottom = round(self.results_sorted["Top"][index], 4)
            else:
                if index == 0:
                    top = round(self.results_sorted["Top"][index], 4)
                    bottom = round(self.results_sorted["Bottom"][index], 4)
                else:
                    top = round(self.results_sorted["Bottom"][index - 1], 4)
                    bottom = round(self.results_sorted["Top"][index], 4)
            self.unit_sections[rock]["Intervals"].append(
                [round(self.results_sorted["Top"][index], 4),
                 round(self.results_sorted["Top"][index] + self.results_sorted["thickness"][index], 4)])
            self.unit_sections[rock]["Intervals"].append([top, bottom])
        # print("rho:", self.results_sorted["rho"], len(self.results_sorted["rho"]))
        # print("vP:", self.results_sorted["vP"], len(self.results_sorted["vP"]))
        # print("vS:", self.results_sorted["vS"], len(self.results_sorted["vS"]))
        # print("vP/vS:", self.results_sorted["vPvS"], len(self.results_sorted["vPvS"]))
        # print("phi:", self.results_sorted["phi"], len(self.results_sorted["phi"]))
        # print("K:", self.results_sorted["K"], len(self.results_sorted["K"]))
        # print("G:", self.results_sorted["G"], len(self.results_sorted["G"]))
        # print("Poisson:", self.results_sorted["Poisson"], len(self.results_sorted["Poisson"]))
        # print("GR:", self.results_sorted["GR"], len(self.results_sorted["GR"]))
        # print("PE:", self.results_sorted["PE"], len(self.results_sorted["PE"]))
        # print(self.results_units)
        #
        self.var_rb_stat = tk.IntVar()
        self.var_rb_stat.set(0)
        self.var_rb_geochem = tk.IntVar()
        self.var_rb_geochem.set(3)
        self.var_rb_lith = tk.IntVar()
        self.var_rb_lith.set(5)
        self.current_rb_lith = tk.IntVar()
        #
        self.results_plot = {}
        for rock in self.list_rocks_short:
            self.results_plot[rock] = {}
            for prop in properties:
                self.results_plot[rock][prop] = []
        for index, value in enumerate(self.results_sorted["thickness"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["thickness"].append(value)
        for index, value in enumerate(self.results_sorted["Top"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["Top"].append(value)
        for index, value in enumerate(self.results_sorted["Bottom"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["Bottom"].append(value)
        for index, value in enumerate(self.results_sorted["rho"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["rho"].append(value)
        for index, value in enumerate(self.results_sorted["vP"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["vP"].append(value)
        for index, value in enumerate(self.results_sorted["vS"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["vS"].append(value)
        for index, value in enumerate(self.results_sorted["vPvS"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["vPvS"].append(value)
        for index, value in enumerate(self.results_sorted["phi"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["phi"].append(value)
        for index, value in enumerate(self.results_sorted["K"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["K"].append(value)
        for index, value in enumerate(self.results_sorted["G"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["G"].append(value)
        for index, value in enumerate(self.results_sorted["Poisson"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["Poisson"].append(value)
        for index, value in enumerate(self.results_sorted["GR"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["GR"].append(value)
        for index, value in enumerate(self.results_sorted["PE"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["PE"].append(value)
        #
        rb_01 = SE(parent=self.parent_subsurface, row_id=29, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                   fg="black").create_radiobutton(var_rb=self.var_rb_stat, value_rb=0,
                                                  text="Well Log Plot", color_bg=self.color_acc_01,
                                                  command=lambda var_rb=self.var_rb_stat: self.change_radiobutton(
                                                      var_rb))
        rb_02 = SE(parent=self.parent_subsurface, row_id=30, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                   fg="black").create_radiobutton(var_rb=self.var_rb_stat, value_rb=1,
                                                  text="Histogram Plot", color_bg=self.color_acc_01,
                                                  command=lambda var_rb=self.var_rb_stat: self.change_radiobutton(
                                                      var_rb))
        rb_03 = SE(parent=self.parent_subsurface, row_id=31, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                   fg="black").create_radiobutton(var_rb=self.var_rb_stat, value_rb=2,
                                                  text="Scatter Plot", color_bg=self.color_acc_01,
                                                  command=lambda var_rb=self.var_rb_stat: self.change_radiobutton(
                                                      var_rb))
        rb_04 = SE(parent=self.parent_subsurface, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                   fg="black").create_radiobutton(var_rb=self.var_rb_geochem, value_rb=3,
                                                  text="Elements", color_bg=self.color_acc_01,
                                                  command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(
                                                      var_rb))
        rb_05 = SE(parent=self.parent_subsurface, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                   fg="black").create_radiobutton(var_rb=self.var_rb_geochem, value_rb=4,
                                                  text="Minerals", color_bg=self.color_acc_01,
                                                  command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(
                                                      var_rb))
        self.gui_elements.extend([rb_01, rb_02, rb_03, rb_04, rb_05])
        for index, rock in enumerate(self.list_rocks_short, start=0):
            rb = SE(parent=self.parent_subsurface, row_id=36 + index, column_id=0, n_rows=1, n_columns=2,
                    bg=self.color_acc_01,
                    fg="black").create_radiobutton(var_rb=self.var_rb_lith, value_rb=index,
                                                   text=rock, color_bg=self.color_acc_01,
                                                   command=lambda var_rb=self.var_rb_lith: self.change_radiobutton(
                                                       var_rb))
            self.gui_elements.append(rb)
        #
        self.entr_list_min = []
        self.entr_list_max = []
        self.entr_list_mean = []
        self.entr_list_std = []
        #
        n_elements = 0
        for rock, element_data in self.data_chemistry.items():
            n_element_rock = len(element_data.keys())
            if n_element_rock > n_elements:
                n_elements = n_element_rock
        #
        for i in range(10 + n_elements):
            self.entr_list_min.append(tk.IntVar())
            self.entr_list_max.append(tk.IntVar())
            self.entr_list_mean.append(tk.IntVar())
            self.entr_list_std.append(tk.IntVar())
        #
        self.results_0 = self.results_units[self.list_rocks_short[0]]
        #
        ## Entry Table
        categories = ["rho", "vP", "vS", "vPvS", "K", "G", "Poisson", "phi", "GR", "PE"]
        for i, category in enumerate(categories):
            if i == 7:
                entr_min = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                             var_entr_set=round(np.min(self.results_0[category] * 100), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                             var_entr_set=round(np.max(self.results_0[category] * 100), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                              var_entr_set=round(np.mean(self.results_0[category] * 100), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                             var_entr_set=round(np.std(self.results_0[category] * 100, ddof=1),
                                                                                3))
            else:
                entr_min = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i],
                                                             var_entr_set=round(np.min(self.results_0[category]), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i],
                                                             var_entr_set=round(np.max(self.results_0[category]), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i],
                                                              var_entr_set=round(np.mean(self.results_0[category]), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=4 + 2 * i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i],
                                                             var_entr_set=round(np.std(self.results_0[category], ddof=1), 3))
            #
            self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                 fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        #
        self.lbl_w["chemistry"].append(lbl)
        #
        self.list_elements_0 = list(self.data_chemistry[self.list_rocks_short[0]].keys())
        self.elements_0 = self.data_chemistry[self.list_rocks_short[0]]
        #
        for index, element in enumerate(self.list_elements_0, start=0):
            if element not in ["U"]:
                lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element), relief=tk.RAISED)
            else:
                lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
            #
            self.lbl_w["chemistry"].append(lbl)
        #
        for index, element in enumerate(self.list_elements_0, start=10):
            if element not in ["U"]:
                entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(self.elements_0[element]), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(self.elements_0[element]), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(self.elements_0[element]), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(self.elements_0[element], ddof=1), 3))
            else:
                ppm_amounts = np.array(self.elements_0[element])*10000
                entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(ppm_amounts), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(ppm_amounts), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(ppm_amounts), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
            #
            self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        self.create_well_log_plot(
            parent=self.parent_subsurface, data_x=self.results_sorted["GR"], data_y=self.results_sorted["Top"],
            row_id=2, column_id=9, n_rows=45, n_columns=9, unit_sections=self.unit_sections)
        #
    def create_random_sequences(self, thickness, style, n_parts=20):
        #
        self.var_rb_stat = tk.IntVar()
        var_rb_stat_start = 0
        self.var_rb_geochem = tk.IntVar()
        var_rb_geochem_start = 3
        self.var_rb_lith = tk.IntVar()
        var_rb_lith_start = 5
        self.current_rb_lith = tk.IntVar()
        self.name = "Random"
        #
        results_subsurface = {}
        self.results_sorted = {}
        properties = ["thickness", "rock", "rho", "vP", "vS", "vPvS", "phi", "K", "G", "Poisson", "GR", "PE", "Top",
                      "Bottom"]
        for prop in properties:
            self.results_sorted[prop] = []
        #
        if style == "siliciclastic":
            rocks = ["Sandstone", "Shale"]
            basement = ["Granite", "Gabbro", "Diorite"]
            list_rocks = []
            self.list_rocks_short = []
            #
            n_units = rd.randint(10, 20)
            for i in range(n_units):
                if len(list_rocks) < n_units - 1:
                    magicnumber = rd.randint(0, len(rocks)-1)
                    list_rocks.append(rocks[magicnumber])
                    if rocks[magicnumber] not in self.list_rocks_short:
                        self.list_rocks_short.append(rocks[magicnumber])
                else:
                    magicnumber = rd.randint(0, len(basement)-1)
                    list_rocks.append(basement[magicnumber])
                    if basement[magicnumber] not in self.list_rocks_short:
                        self.list_rocks_short.append(basement[magicnumber])
            list_thickness = self.split_thickness(thickness=thickness, n_units=n_units)
            #
            for index, rock in enumerate(list_rocks, start=0):
                results_subsurface[index] = {}
                if rock == "Sandstone":
                    if index == 0:
                        thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                        for i, d in enumerate(thickness_unit, start=0):
                            if i == 0:
                                results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                results_subsurface[index][i]["thickness"] = d
                                self.results_sorted["Top"].append(0)
                                self.results_sorted["thickness"].append(d)
                                self.results_sorted["Bottom"].append(d)
                                self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                                self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                                self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                                self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                                self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                                self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                                self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                                self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                                self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                                self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                                self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                            else:
                                results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                   amounts=results_subsurface[index][i-1]["mineralogy"])
                                results_subsurface[index][i]["thickness"] = d
                                self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                                self.results_sorted["thickness"].append(d)
                                self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                                self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                                self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                                self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                                self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                                self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                                self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                                self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                                self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                                self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                                self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                                self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                    else:
                        thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                        for i, d in enumerate(thickness_unit, start=0):
                            if i == 0:
                                if list_rocks[index-1] == "Shale":
                                    magicnumber = rd.randint(0, 2)
                                    if magicnumber == 0:
                                        results_subsurface[index][i] = sandstone(fluid="gas", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                    elif magicnumber == 1:
                                        results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                    elif magicnumber == 2:
                                        results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                elif list_rocks[index-1] == "Sandstone" and results_subsurface[index-1][i]["fluid"] == "gas":
                                    results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                elif list_rocks[index-1] == "Sandstone" and results_subsurface[index-1][i]["fluid"] == "oil":
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                else:
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                #
                                results_subsurface[index][i]["thickness"] = d
                                self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                                self.results_sorted["thickness"].append(d)
                                self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                                self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                                self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                                self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                                self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                                self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                                self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                                self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                                self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                                self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                                self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                                self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                            else:
                                if list_rocks[index-1] == "Shale":
                                    magicnumber = rd.randint(0, 2)
                                    if magicnumber == 0:
                                        results_subsurface[index][i] = sandstone(fluid="gas", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                         amounts=results_subsurface[index][i-1]["mineralogy"])
                                    elif magicnumber == 1:
                                        results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                         amounts=results_subsurface[index][i-1]["mineralogy"])
                                    if magicnumber == 2:
                                        results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                           amounts=results_subsurface[index][i-1]["mineralogy"])
                                elif list_rocks[index-1] == "Sandstone" and results_subsurface[index-1][i]["fluid"] == "gas":
                                    results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                     amounts=results_subsurface[index][i-1]["mineralogy"])
                                elif list_rocks[index-1] == "Sandstone" and results_subsurface[index-1][i]["fluid"] == "oil":
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                       amounts=results_subsurface[index][i-1]["mineralogy"])
                                else:
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                       amounts=results_subsurface[index][i-1]["mineralogy"])
                                #
                                results_subsurface[index][i]["thickness"] = d
                                self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                                self.results_sorted["thickness"].append(d)
                                self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                                self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                                self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                                self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                                self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                                self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                                self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                                self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                                self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                                self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                                self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                                self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                #
                elif rock == "Shale":
                    if index == 0:
                        thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                        for i, d in enumerate(thickness_unit, start=0):
                            if i == 0:
                                results_subsurface[index][i] = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(0, 0.1))
                                self.results_sorted["Top"].append(0)
                                self.results_sorted["thickness"].append(d)
                                self.results_sorted["Bottom"].append(d)
                                self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                                self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                                self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                                self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                                self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                                self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                                self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                                self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                                self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                                self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                                self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                            else:
                                results_subsurface[index][i] = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(0, 0.1), amounts=results_subsurface[index][i-1]["mineralogy"])
                                self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                                self.results_sorted["thickness"].append(d)
                                self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                                self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                                self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                                self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                                self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                                self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                                self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                                self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                                self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                                self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                                self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                                self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                    else:
                        thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                        for i, d in enumerate(thickness_unit, start=0):
                            if i == 0:
                                results_subsurface[index][i] = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(0, 0.1))
                            else:
                                results_subsurface[index][i] = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(0, 0.1), amounts=results_subsurface[index][i-1]["mineralogy"])
                            #
                            self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                            self.results_sorted["thickness"].append(d)
                            self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                            self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                            self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                            self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                            self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                            self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                            self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                            self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                            self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                            self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                            self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                            self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                #
                elif rock == "Granite":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_granite()
                        else:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_granite(amounts=results_subsurface[index][i-1]["mineralogy"])
                        #
                        self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                        self.results_sorted["thickness"].append(d)
                        self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                        self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                        self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                        self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                        self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                        self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                        self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                        self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                        self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                        self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                        self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                        self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                #
                elif rock == "Gabbro":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_gabbro()
                        else:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_gabbro(amounts=results_subsurface[index][i-1]["mineralogy"])
                        #
                        self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                        self.results_sorted["thickness"].append(d)
                        self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                        self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                        self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                        self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                        self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                        self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                        self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                        self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                        self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                        self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                        self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                        self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                #
                elif rock == "Diorite":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=n_parts)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True,
                                                                    porosity=rd.uniform(0, 0.05)).create_simple_diorite()
                        else:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True,
                                                                    porosity=rd.uniform(0, 0.05)).create_simple_diorite(amounts=results_subsurface[index][i-1]["mineralogy"])
                        #
                        self.results_sorted["Top"].append(self.results_sorted["Bottom"][-1])
                        self.results_sorted["thickness"].append(d)
                        self.results_sorted["Bottom"].append(int(self.results_sorted["Bottom"][-1] + d))
                        self.results_sorted["rock"].append(results_subsurface[index][i]["rock"])
                        self.results_sorted["rho"].append(results_subsurface[index][i]["rho"])
                        self.results_sorted["vP"].append(results_subsurface[index][i]["vP"])
                        self.results_sorted["vS"].append(results_subsurface[index][i]["vS"])
                        self.results_sorted["vPvS"].append(results_subsurface[index][i]["vP/vS"])
                        self.results_sorted["phi"].append(results_subsurface[index][i]["phi"])
                        self.results_sorted["K"].append(results_subsurface[index][i]["K"])
                        self.results_sorted["G"].append(results_subsurface[index][i]["G"])
                        self.results_sorted["Poisson"].append(results_subsurface[index][i]["nu"])
                        self.results_sorted["GR"].append(results_subsurface[index][i]["GR"])
                        self.results_sorted["PE"].append(results_subsurface[index][i]["PE"])
                #
                results_subsurface[index][i]["thickness"] = list_thickness[index]
        #
        self.results_plot = {}
        for rock in self.list_rocks_short:
            self.results_plot[rock] = {}
            for prop in properties:
                self.results_plot[rock][prop] = []
        for index, value in enumerate(self.results_sorted["thickness"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["thickness"].append(value)
        for index, value in enumerate(self.results_sorted["Top"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["Top"].append(value)
        for index, value in enumerate(self.results_sorted["Bottom"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["Bottom"].append(value)
        for index, value in enumerate(self.results_sorted["rho"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["rho"].append(value)
        for index, value in enumerate(self.results_sorted["vP"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["vP"].append(value)
        for index, value in enumerate(self.results_sorted["vS"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["vS"].append(value)
        for index, value in enumerate(self.results_sorted["vPvS"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["vPvS"].append(value)
        for index, value in enumerate(self.results_sorted["phi"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["phi"].append(value)
        for index, value in enumerate(self.results_sorted["K"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["K"].append(value)
        for index, value in enumerate(self.results_sorted["G"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["G"].append(value)
        for index, value in enumerate(self.results_sorted["Poisson"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["Poisson"].append(value)
        for index, value in enumerate(self.results_sorted["GR"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["GR"].append(value)
        for index, value in enumerate(self.results_sorted["PE"], start=0):
            self.results_plot[self.results_sorted["rock"][index]]["PE"].append(value)
        #
        rb_01 = SE(parent=self.parent_subsurface, row_id=29, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_stat, var_rb_set=var_rb_stat_start, value_rb=0,
                                          text="Well Log Plot", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_stat: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_subsurface, row_id=30, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_stat, var_rb_set=var_rb_stat_start, value_rb=1,
                                          text="Histogram Plot", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_stat: self.change_radiobutton(var_rb))
        rb_03 = SE(parent=self.parent_subsurface, row_id=31, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_stat, var_rb_set=var_rb_stat_start, value_rb=2,
                                          text="Scatter Plot", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_stat: self.change_radiobutton(var_rb))
        rb_04 = SE(parent=self.parent_subsurface, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, var_rb_set=var_rb_geochem_start, value_rb=3,
                                          text="Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        rb_05 = SE(parent=self.parent_subsurface, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, var_rb_set=var_rb_geochem_start, value_rb=4,
                                          text="Minerals", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        self.gui_elements.extend([rb_01, rb_02, rb_03, rb_04, rb_05])
        for index, rock in enumerate(self.list_rocks_short, start=0):
            rb = SE(parent=self.parent_subsurface, row_id=36+index, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                    fg="black").create_radiobutton(var_rb=self.var_rb_lith, var_rb_set=var_rb_lith_start, value_rb=var_rb_lith_start+index,
                                                   text=rock, color_bg=self.color_acc_01,
                                                   command=lambda var_rb=self.var_rb_lith: self.change_radiobutton(var_rb))
            self.gui_elements.append(rb)
        #
        data_all = {}
        for rock in self.list_rocks_short:
            data_all[rock] = []
        for item in results_subsurface.values():
            for item_part in item.values():
                data_all[item_part["rock"]].append(item_part)
        #
        #self.extract_data(data=data_all)
        #
        if "Sandstone" in self.list_rocks_short:
            self.rho_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="rho")
            self.vP_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="vP")
            self.vS_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="vS")
            self.vPvS_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="vP/vS")
            self.bulk_mod_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="K")
            self.shear_mod_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="G")
            self.youngs_mod_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="E")
            self.poisson_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="nu")
            self.phi_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="phi")
            self.gamma_ray_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="GR")
            self.photoelectricity_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="PE")
            self.chemistry_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="chemistry")
            self.mineralogy_sst = DP(dataset=data_all["Sandstone"]).extract_data(keyword="mineralogy")
            #
            self.list_elements_sst = list(self.chemistry_sst[0].keys())
            self.list_minerals_sst = list(self.mineralogy_sst[0].keys())
            self.elements_sst = {}
            self.minerals_sst = {}
            for element in self.list_elements_sst:
                self.elements_sst[element] = []
                for chemistry_data in self.chemistry_sst:
                    self.elements_sst[element].append(abs(chemistry_data[element]*100))
            for mineral in self.list_minerals_sst:
                self.minerals_sst[mineral] = []
                for mineralogy_data in self.mineralogy_sst:
                    self.minerals_sst[mineral].append(abs(mineralogy_data[mineral]*100))
            #
            self.results_sst = [self.rho_sst, self.vP_sst, self.vS_sst, self.vPvS_sst, self.bulk_mod_sst,
                                self.shear_mod_sst, self.poisson_sst, self.phi_sst, self.gamma_ray_sst,
                                self.photoelectricity_sst]
        #
        if "Shale" in self.list_rocks_short:
            self.rho_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="rho")
            self.vP_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="vP")
            self.vS_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="vS")
            self.vPvS_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="vP/vS")
            self.bulk_mod_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="K")
            self.shear_mod_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="G")
            self.youngs_mod_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="E")
            self.poisson_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="nu")
            self.phi_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="phi")
            self.gamma_ray_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="GR")
            self.photoelectricity_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="PE")
            self.chemistry_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="chemistry")
            self.mineralogy_sh = DP(dataset=data_all["Shale"]).extract_data(keyword="mineralogy")
            #
            self.list_elements_sh = list(self.chemistry_sh[0].keys())
            self.list_minerals_sh = list(self.mineralogy_sh[0].keys())
            self.elements_sh = {}
            self.minerals_sh = {}
            for element in self.list_elements_sh:
                self.elements_sh[element] = []
                for chemistry_data in self.chemistry_sh:
                    self.elements_sh[element].append(abs(chemistry_data[element]*100))
            for mineral in self.list_minerals_sh:
                self.minerals_sh[mineral] = []
                for mineralogy_data in self.mineralogy_sh:
                    self.minerals_sh[mineral].append(abs(mineralogy_data[mineral]*100))
            #
            self.results_sh = [self.rho_sh, self.vP_sh, self.vS_sh, self.vPvS_sh, self.bulk_mod_sh, self.shear_mod_sh,
                               self.poisson_sh, self.phi_sh, self.gamma_ray_sh, self.photoelectricity_sh]
        #
        if "Granite" in self.list_rocks_short:
            self.rho_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="rho")
            self.vP_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="vP")
            self.vS_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="vS")
            self.vPvS_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="vP/vS")
            self.bulk_mod_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="K")
            self.shear_mod_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="G")
            self.youngs_mod_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="E")
            self.poisson_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="nu")
            self.phi_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="phi")
            self.gamma_ray_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="GR")
            self.photoelectricity_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="PE")
            self.chemistry_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="chemistry")
            self.mineralogy_granite = DP(dataset=data_all["Granite"]).extract_data(keyword="mineralogy")
            #
            self.list_elements_granite = list(self.chemistry_granite[0].keys())
            self.list_minerals_granite = list(self.mineralogy_granite[0].keys())
            self.elements_granite = {}
            self.minerals_granite = {}
            for element in self.list_elements_granite:
                self.elements_granite[element] = []
                for chemistry_data in self.chemistry_granite:
                    self.elements_granite[element].append(abs(chemistry_data[element]*100))
            for mineral in self.list_minerals_granite:
                self.minerals_granite[mineral] = []
                for mineralogy_data in self.mineralogy_granite:
                    self.minerals_granite[mineral].append(abs(mineralogy_data[mineral]*100))
            #
            self.results_granite = [self.rho_granite, self.vP_granite, self.vS_granite, self.vPvS_granite,
                                    self.bulk_mod_granite, self.shear_mod_granite, self.poisson_granite,
                                    self.phi_granite, self.gamma_ray_granite, self.photoelectricity_granite]
            #
        elif "Gabbro" in self.list_rocks_short:
            self.rho_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="rho")
            self.vP_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="vP")
            self.vS_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="vS")
            self.vPvS_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="vP/vS")
            self.bulk_mod_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="K")
            self.shear_mod_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="G")
            self.youngs_mod_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="E")
            self.poisson_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="nu")
            self.phi_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="phi")
            self.gamma_ray_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="GR")
            self.photoelectricity_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="PE")
            self.chemistry_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="chemistry")
            self.mineralogy_gabbro = DP(dataset=data_all["Gabbro"]).extract_data(keyword="mineralogy")
            #
            self.list_elements_gabbro = list(self.chemistry_gabbro[0].keys())
            self.list_minerals_gabbro = list(self.mineralogy_gabbro[0].keys())
            self.elements_gabbro = {}
            self.minerals_gabbro = {}
            for element in self.list_elements_gabbro:
                self.elements_gabbro[element] = []
                for chemistry_data in self.chemistry_gabbro:
                    self.elements_gabbro[element].append(abs(chemistry_data[element]*100))
            for mineral in self.list_minerals_gabbro:
                self.minerals_gabbro[mineral] = []
                for mineralogy_data in self.mineralogy_gabbro:
                    self.minerals_gabbro[mineral].append(abs(mineralogy_data[mineral]*100))
            #
            self.results_gabbro = [self.rho_gabbro, self.vP_gabbro, self.vS_gabbro, self.vPvS_gabbro,
                                   self.bulk_mod_gabbro, self.shear_mod_gabbro, self.poisson_gabbro, self.phi_gabbro,
                                   self.gamma_ray_gabbro, self.photoelectricity_gabbro]
            #
        elif "Diorite" in self.list_rocks_short:
            self.rho_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="rho")
            self.vP_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="vP")
            self.vS_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="vS")
            self.vPvS_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="vP/vS")
            self.bulk_mod_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="K")
            self.shear_mod_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="G")
            self.youngs_mod_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="E")
            self.poisson_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="nu")
            self.phi_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="phi")
            self.gamma_ray_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="GR")
            self.photoelectricity_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="PE")
            self.chemistry_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="chemistry")
            self.mineralogy_diorite = DP(dataset=data_all["Diorite"]).extract_data(keyword="mineralogy")
            #
            self.list_elements_diorite = list(self.chemistry_diorite[0].keys())
            self.list_minerals_diorite = list(self.mineralogy_diorite[0].keys())
            self.elements_diorite = {}
            self.minerals_diorite = {}
            for element in self.list_elements_diorite:
                self.elements_diorite[element] = []
                for chemistry_data in self.chemistry_diorite:
                    self.elements_diorite[element].append(abs(chemistry_data[element]*100))
            for mineral in self.list_minerals_diorite:
                self.minerals_diorite[mineral] = []
                for mineralogy_data in self.mineralogy_diorite:
                    self.minerals_diorite[mineral].append(abs(mineralogy_data[mineral]*100))
            #
            self.results_diorite = [self.rho_diorite, self.vP_diorite, self.vS_diorite, self.vPvS_diorite,
                                    self.bulk_mod_diorite, self.shear_mod_diorite, self.poisson_diorite,
                                    self.phi_diorite, self.gamma_ray_diorite, self.photoelectricity_diorite]
        #
        self.entr_list_min = []
        self.entr_list_max = []
        self.entr_list_mean = []
        self.entr_list_std = []
        #
        if "Sandstone" in self.list_rocks_short[0]:
            self.var_rb_lith.set(var_rb_lith_start)
            self.current_rb_lith.set(var_rb_lith_start)
            self.list_elements_0 = self.list_elements_sst
            self.list_minerals_0 = self.list_minerals_sst
            self.results_0 = self.results_sst
            self.elements_0 = self.elements_sst
            self.minerals_0 = self.minerals_sst
        elif "Shale" in self.list_rocks_short[0]:
            self.var_rb_lith.set(var_rb_lith_start)
            self.current_rb_lith.set(var_rb_lith_start)
            self.list_elements_0 = self.list_elements_sh
            self.list_minerals_0 = self.list_minerals_sh
            self.results_0 = self.results_sh
            self.elements_0 = self.elements_sh
            self.minerals_0 = self.minerals_sh
        #
        n_elements = 0
        if len(self.list_elements_sst) > n_elements:
            n_elements = len(self.list_elements_sst)
        if len(self.list_elements_sh) > n_elements:
            n_elements = len(self.list_elements_sh)
        try:
            if len(self.list_elements_granite) > n_elements:
                n_elements = len(self.list_elements_granite)
        except:
            pass
        try:
            if len(self.list_elements_gabbro) > n_elements:
                n_elements = len(self.list_elements_gabbro)
        except:
            pass
        try:
            if len(self.list_elements_diorite) > n_elements:
                n_elements = len(self.list_elements_diorite)
        except:
            pass
        #
        for i in range(10+n_elements):
            self.entr_list_min.append(tk.IntVar())
            self.entr_list_max.append(tk.IntVar())
            self.entr_list_mean.append(tk.IntVar())
            self.entr_list_std.append(tk.IntVar())
        #
        ## Entry Table
        for i in range(10):
            if i == 7:
                entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]*100), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]*100), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]*100), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i]*100, ddof=1), 3))
            else:
                entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i], ddof=1), 3))
            #
            self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
               fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        self.lbl_w["chemistry"].append(lbl)
        for index, element in enumerate(self.list_elements_0, start=0):
            if element not in ["U"]:
                lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element), relief=tk.RAISED)
            else:
                lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
        #
        for index, element in enumerate(self.list_elements_0, start=10):
            if element not in ["U"]:
                entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(self.elements_0[element]), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(self.elements_0[element]), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(self.elements_0[element]), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(self.elements_0[element], ddof=1), 3))
            else:
                ppm_amounts = np.array(self.elements_0[element])*10000
                entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                             var_entr_set=round(np.min(ppm_amounts), 3))
                entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                             var_entr_set=round(np.max(ppm_amounts), 3))
                entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                               fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                              var_entr_set=round(np.mean(ppm_amounts), 3))
                entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                              fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                             var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
            #
            self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        self.create_well_log_plot(parent=self.parent_subsurface, data_x=self.results_sorted["GR"],
                                  data_y=self.results_sorted["Top"], row_id=2, column_id=9, n_rows=45, n_columns=9)
        #
    def change_radiobutton(self, var_rb):
        if var_rb.get() == 0:   # Well Log Plot
            try:
                for rb in self.rb_prop:
                    rb.grid_forget()
                self.rb_prop.clear()
            except:
                pass
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            if self.name == "Random":
                self.create_well_log_plot(
                    parent=self.parent_subsurface, data_x=self.results_sorted["GR"], data_y=self.results_sorted["Top"],
                    row_id=2, column_id=9, n_rows=45, n_columns=9)
            elif self.name in ["Zechstein", "Muschelkalk"]:
                self.create_well_log_plot(
                    parent=self.parent_subsurface, data_x=self.results_sorted["GR"], data_y=self.results_sorted["Top"],
                    row_id=2, column_id=9, n_rows=45, n_columns=9, unit_sections=self.unit_sections)
            #
        elif var_rb.get() == 1: # Histogram Plot
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            if self.name == "Random":
                data_rho = [np.array(self.results_plot["Sandstone"]["phi"])*100,
                            np.array(self.results_plot["Shale"]["phi"])*100,
                            np.array(self.results_plot[self.list_rocks_short[-1]]["phi"])*100]
                data_vP = [self.results_plot["Sandstone"]["vP"], self.results_plot["Shale"]["vP"],
                           self.results_plot[self.list_rocks_short[-1]]["vP"]]
                data_vS = [self.results_plot["Sandstone"]["vS"], self.results_plot["Shale"]["vS"],
                           self.results_plot[self.list_rocks_short[-1]]["vS"]]
                data_vPvS = [self.results_plot["Sandstone"]["vPvS"], self.results_plot["Shale"]["vPvS"],
                             self.results_plot[self.list_rocks_short[-1]]["vPvS"]]
                data_K = [self.results_plot["Sandstone"]["K"], self.results_plot["Shale"]["K"],
                          self.results_plot[self.list_rocks_short[-1]]["K"]]
                data_G = [self.results_plot["Sandstone"]["G"], self.results_plot["Shale"]["G"],
                          self.results_plot[self.list_rocks_short[-1]]["G"]]
                data_nu = [self.results_plot["Sandstone"]["Poisson"], self.results_plot["Shale"]["Poisson"],
                           self.results_plot[self.list_rocks_short[-1]]["Poisson"]]
                data_phi = [np.array(self.results_plot["Sandstone"]["phi"])*100,
                            np.array(self.results_plot["Shale"]["phi"])*100,
                            np.array(self.results_plot[self.list_rocks_short[-1]]["phi"])*100]
                data_gr = [self.results_plot["Sandstone"]["GR"], self.results_plot["Shale"]["GR"],
                           self.results_plot[self.list_rocks_short[-1]]["GR"]]
                data_pe = [self.results_plot["Sandstone"]["PE"], self.results_plot["Shale"]["PE"],
                           self.results_plot[self.list_rocks_short[-1]]["PE"]]
                #
                colors = ["tan", "olivedrab", "darkorange"]
                litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
            elif self.name == "Zechstein":
                data_rho = [np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Rock Salt"]["phi"]) * 100,
                            np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Kupferschiefer"]["phi"]) * 100,
                            np.array(self.results_plot["Potash"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100]
                data_vP = [self.results_plot["Anhydrite"]["vP"], self.results_plot["Dolomite"]["vP"],
                           self.results_plot["Rock Salt"]["vP"], self.results_plot["Limestone"]["vP"],
                           self.results_plot["Kupferschiefer"]["vP"], self.results_plot["Potash"]["vP"],
                           self.results_plot["Mudstone"]["vP"]]
                data_vS = [self.results_plot["Anhydrite"]["vS"], self.results_plot["Dolomite"]["vS"],
                           self.results_plot["Rock Salt"]["vS"], self.results_plot["Limestone"]["vS"],
                           self.results_plot["Kupferschiefer"]["vS"], self.results_plot["Potash"]["vS"],
                           self.results_plot["Mudstone"]["vS"]]
                data_vPvS = [self.results_plot["Anhydrite"]["vPvS"], self.results_plot["Dolomite"]["vPvS"],
                           self.results_plot["Rock Salt"]["vPvS"], self.results_plot["Limestone"]["vPvS"],
                           self.results_plot["Kupferschiefer"]["vPvS"], self.results_plot["Potash"]["vPvS"],
                             self.results_plot["Mudstone"]["vPvS"]]
                data_K = [self.results_plot["Anhydrite"]["K"], self.results_plot["Dolomite"]["K"],
                           self.results_plot["Rock Salt"]["K"], self.results_plot["Limestone"]["K"],
                           self.results_plot["Kupferschiefer"]["K"], self.results_plot["Potash"]["K"],
                          self.results_plot["Mudstone"]["K"]]
                data_G = [self.results_plot["Anhydrite"]["G"], self.results_plot["Dolomite"]["G"],
                           self.results_plot["Rock Salt"]["G"], self.results_plot["Limestone"]["G"],
                           self.results_plot["Kupferschiefer"]["G"], self.results_plot["Potash"]["G"],
                          self.results_plot["Mudstone"]["G"]]
                data_nu = [self.results_plot["Anhydrite"]["Poisson"], self.results_plot["Dolomite"]["Poisson"],
                           self.results_plot["Rock Salt"]["Poisson"], self.results_plot["Limestone"]["Poisson"],
                           self.results_plot["Kupferschiefer"]["Poisson"], self.results_plot["Potash"]["Poisson"],
                           self.results_plot["Mudstone"]["Poisson"]]
                data_phi = [np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Rock Salt"]["phi"]) * 100,
                            np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Kupferschiefer"]["phi"]) * 100,
                            np.array(self.results_plot["Potash"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100]
                data_gr = [self.results_plot["Anhydrite"]["GR"], self.results_plot["Dolomite"]["GR"],
                           self.results_plot["Rock Salt"]["GR"], self.results_plot["Limestone"]["GR"],
                           self.results_plot["Kupferschiefer"]["GR"], self.results_plot["Potash"]["GR"],
                           self.results_plot["Mudstone"]["GR"]]
                data_pe = [self.results_plot["Anhydrite"]["PE"], self.results_plot["Dolomite"]["PE"],
                           self.results_plot["Rock Salt"]["PE"], self.results_plot["Limestone"]["PE"],
                           self.results_plot["Kupferschiefer"]["PE"], self.results_plot["Potash"]["PE"],
                           self.results_plot["Mudstone"]["PE"]]
                #
                colors = ["orchid", "lightcyan", "lavender", "skyblue", "gray", "yellowgreen", "olivedrab"]
                litho_list = ["Anhydrite", "Dolomite", "Rock Salt", "Limestone", "Kupferschiefer", "Potash", "Mudstone"]
                #
            elif self.name == "Muschelkalk":
                data_rho = [np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100,
                            np.array(self.results_plot["Marl"]["phi"]) * 100]
                data_vP = [self.results_plot["Limestone"]["vP"], self.results_plot["Dolomite"]["vP"],
                           self.results_plot["Anhydrite"]["vP"], self.results_plot["Mudstone"]["vP"],
                           self.results_plot["Marl"]["vP"]]
                data_vS = [self.results_plot["Limestone"]["vS"], self.results_plot["Dolomite"]["vS"],
                           self.results_plot["Anhydrite"]["vS"], self.results_plot["Mudstone"]["vS"],
                           self.results_plot["Marl"]["vS"]]
                data_vPvS = [self.results_plot["Limestone"]["vPvS"], self.results_plot["Dolomite"]["vPvS"],
                           self.results_plot["Anhydrite"]["vPvS"], self.results_plot["Mudstone"]["vPvS"],
                           self.results_plot["Marl"]["vPvS"]]
                data_K = [self.results_plot["Limestone"]["K"], self.results_plot["Dolomite"]["K"],
                           self.results_plot["Anhydrite"]["K"], self.results_plot["Mudstone"]["K"],
                           self.results_plot["Marl"]["K"]]
                data_G = [self.results_plot["Limestone"]["G"], self.results_plot["Dolomite"]["G"],
                           self.results_plot["Anhydrite"]["G"], self.results_plot["Mudstone"]["G"],
                           self.results_plot["Marl"]["G"]]
                data_nu = [self.results_plot["Limestone"]["Poisson"], self.results_plot["Dolomite"]["Poisson"],
                           self.results_plot["Anhydrite"]["Poisson"], self.results_plot["Mudstone"]["Poisson"],
                           self.results_plot["Marl"]["Poisson"]]
                data_phi = [np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100,
                            np.array(self.results_plot["Marl"]["phi"]) * 100]
                data_gr = [self.results_plot["Limestone"]["GR"], self.results_plot["Dolomite"]["GR"],
                           self.results_plot["Anhydrite"]["GR"], self.results_plot["Mudstone"]["GR"],
                           self.results_plot["Marl"]["GR"]]
                data_pe = [self.results_plot["Limestone"]["PE"], self.results_plot["Dolomite"]["PE"],
                           self.results_plot["Anhydrite"]["PE"], self.results_plot["Mudstone"]["PE"],
                           self.results_plot["Marl"]["PE"]]
                #
                colors = ["skyblue", "lightcyan", "orchid", "olivedrab", "moccasin"]
                litho_list = ["Limestone", "Dolomite", "Anhydrite", "Mudstone", "Marl"]
                #
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_rho, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Density $\\varrho$ (g/ccm)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            #
            self.create_3x3_histo(parent=self.parent_subsurface, data=data, row_id=2, column_id=9, n_rows=45,
                                  n_columns=9, colors=colors, labels=labels, lithos=litho_list)
            #
        elif var_rb.get() == 2: # Scatter Plot
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
                for rb in self.rb_prop:
                    rb.grid_forget()
                self.rb_prop.clear()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            self.var_prop = tk.IntVar()
            self.var_prop.set(8)
            rb_rho = SE(parent=self.parent_subsurface, row_id=40, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, value_rb=8,
                                                       text="Density", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            rb_phi = SE(parent=self.parent_subsurface, row_id=41, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, value_rb=9,
                                                       text="Porosity", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            self.rb_prop = [rb_rho, rb_phi]
            #
            if self.name == "Random":
                data_x_rho = [np.array(self.results_plot["Sandstone"]["rho"])/1000,
                              np.array(self.results_plot["Shale"]["rho"])/1000,
                              np.array(self.results_plot[self.list_rocks_short[-1]]["rho"])/1000]
                #
                data_vP = [self.results_plot["Sandstone"]["vP"], self.results_plot["Shale"]["vP"],
                           self.results_plot[self.list_rocks_short[-1]]["vP"]]
                data_vS = [self.results_plot["Sandstone"]["vS"], self.results_plot["Shale"]["vS"],
                           self.results_plot[self.list_rocks_short[-1]]["vS"]]
                data_vPvS = [self.results_plot["Sandstone"]["vPvS"], self.results_plot["Shale"]["vPvS"],
                             self.results_plot[self.list_rocks_short[-1]]["vPvS"]]
                data_K = [self.results_plot["Sandstone"]["K"], self.results_plot["Shale"]["K"],
                          self.results_plot[self.list_rocks_short[-1]]["K"]]
                data_G = [self.results_plot["Sandstone"]["G"], self.results_plot["Shale"]["G"],
                          self.results_plot[self.list_rocks_short[-1]]["G"]]
                data_nu = [self.results_plot["Sandstone"]["Poisson"], self.results_plot["Shale"]["Poisson"],
                           self.results_plot[self.list_rocks_short[-1]]["Poisson"]]
                data_phi = [np.array(self.results_plot["Sandstone"]["phi"])*100,
                            np.array(self.results_plot["Shale"]["phi"])*100,
                            np.array(self.results_plot[self.list_rocks_short[-1]]["phi"])*100]
                data_gr = [self.results_plot["Sandstone"]["GR"], self.results_plot["Shale"]["GR"],
                           self.results_plot[self.list_rocks_short[-1]]["GR"]]
                data_pe = [self.results_plot["Sandstone"]["PE"], self.results_plot["Shale"]["PE"],
                           self.results_plot[self.list_rocks_short[-1]]["PE"]]
                #
                colors = ["tan", "olivedrab", "darkorange"]
                litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
            elif self.name == "Zechstein":
                data_x_rho = [np.array(self.results_plot["Anhydrite"]["rho"]) / 1000,
                              np.array(self.results_plot["Dolomite"]["rho"]) / 1000,
                              np.array(self.results_plot["Rock Salt"]["rho"]) / 1000,
                              np.array(self.results_plot["Limestone"]["rho"]) / 1000,
                              np.array(self.results_plot["Kupferschiefer"]["rho"]) / 1000,
                              np.array(self.results_plot["Potash"]["rho"]) / 1000,
                              np.array(self.results_plot["Mudstone"]["rho"]) / 1000]
                #
                data_vP = [self.results_plot["Anhydrite"]["vP"], self.results_plot["Dolomite"]["vP"],
                           self.results_plot["Rock Salt"]["vP"], self.results_plot["Limestone"]["vP"],
                           self.results_plot["Kupferschiefer"]["vP"], self.results_plot["Potash"]["vP"],
                           self.results_plot["Mudstone"]["vP"]]
                data_vS = [self.results_plot["Anhydrite"]["vS"], self.results_plot["Dolomite"]["vS"],
                           self.results_plot["Rock Salt"]["vS"], self.results_plot["Limestone"]["vS"],
                           self.results_plot["Kupferschiefer"]["vS"], self.results_plot["Potash"]["vS"],
                           self.results_plot["Mudstone"]["vS"]]
                data_vPvS = [self.results_plot["Anhydrite"]["vPvS"], self.results_plot["Dolomite"]["vPvS"],
                             self.results_plot["Rock Salt"]["vPvS"], self.results_plot["Limestone"]["vPvS"],
                             self.results_plot["Kupferschiefer"]["vPvS"], self.results_plot["Potash"]["vPvS"],
                             self.results_plot["Mudstone"]["vPvS"]]
                data_K = [self.results_plot["Anhydrite"]["K"], self.results_plot["Dolomite"]["K"],
                          self.results_plot["Rock Salt"]["K"], self.results_plot["Limestone"]["K"],
                          self.results_plot["Kupferschiefer"]["K"], self.results_plot["Potash"]["K"],
                          self.results_plot["Mudstone"]["K"]]
                data_G = [self.results_plot["Anhydrite"]["G"], self.results_plot["Dolomite"]["G"],
                          self.results_plot["Rock Salt"]["G"], self.results_plot["Limestone"]["G"],
                          self.results_plot["Kupferschiefer"]["G"], self.results_plot["Potash"]["G"],
                          self.results_plot["Mudstone"]["G"]]
                data_nu = [self.results_plot["Anhydrite"]["Poisson"], self.results_plot["Dolomite"]["Poisson"],
                           self.results_plot["Rock Salt"]["Poisson"], self.results_plot["Limestone"]["Poisson"],
                           self.results_plot["Kupferschiefer"]["Poisson"], self.results_plot["Potash"]["Poisson"],
                           self.results_plot["Mudstone"]["Poisson"]]
                data_phi = [np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Rock Salt"]["phi"]) * 100,
                            np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Kupferschiefer"]["phi"]) * 100,
                            np.array(self.results_plot["Potash"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100]
                data_gr = [self.results_plot["Anhydrite"]["GR"], self.results_plot["Dolomite"]["GR"],
                           self.results_plot["Rock Salt"]["GR"], self.results_plot["Limestone"]["GR"],
                           self.results_plot["Kupferschiefer"]["GR"], self.results_plot["Potash"]["GR"],
                           self.results_plot["Mudstone"]["GR"]]
                data_pe = [self.results_plot["Anhydrite"]["PE"], self.results_plot["Dolomite"]["PE"],
                           self.results_plot["Rock Salt"]["PE"], self.results_plot["Limestone"]["PE"],
                           self.results_plot["Kupferschiefer"]["PE"], self.results_plot["Potash"]["PE"],
                           self.results_plot["Mudstone"]["PE"]]
                #
                colors = ["orchid", "lightcyan", "lavender", "skyblue", "gray", "yellowgreen", "olivedrab"]
                litho_list = ["Anhydrite", "Dolomite", "Rock Salt", "Limestone", "Kupferschiefer", "Potash", "Mudstone"]
                #
            elif self.name == "Muschelkalk":
                data_x_rho = [np.array(self.results_plot["Limestone"]["rho"]) / 1000,
                              np.array(self.results_plot["Dolomite"]["rho"]) / 1000,
                              np.array(self.results_plot["Anhydrite"]["rho"]) / 1000,
                              np.array(self.results_plot["Mudstone"]["rho"]) / 1000,
                              np.array(self.results_plot["Marl"]["rho"]) / 1000]
                #
                data_vP = [self.results_plot["Limestone"]["vP"], self.results_plot["Dolomite"]["vP"],
                           self.results_plot["Anhydrite"]["vP"], self.results_plot["Mudstone"]["vP"],
                           self.results_plot["Marl"]["vP"]]
                data_vS = [self.results_plot["Limestone"]["vS"], self.results_plot["Dolomite"]["vS"],
                           self.results_plot["Anhydrite"]["vS"], self.results_plot["Mudstone"]["vS"],
                           self.results_plot["Marl"]["vS"]]
                data_vPvS = [self.results_plot["Limestone"]["vPvS"], self.results_plot["Dolomite"]["vPvS"],
                           self.results_plot["Anhydrite"]["vPvS"], self.results_plot["Mudstone"]["vPvS"],
                           self.results_plot["Marl"]["vPvS"]]
                data_K = [self.results_plot["Limestone"]["K"], self.results_plot["Dolomite"]["K"],
                           self.results_plot["Anhydrite"]["K"], self.results_plot["Mudstone"]["K"],
                           self.results_plot["Marl"]["K"]]
                data_G = [self.results_plot["Limestone"]["G"], self.results_plot["Dolomite"]["G"],
                           self.results_plot["Anhydrite"]["G"], self.results_plot["Mudstone"]["G"],
                           self.results_plot["Marl"]["G"]]
                data_nu = [self.results_plot["Limestone"]["Poisson"], self.results_plot["Dolomite"]["Poisson"],
                           self.results_plot["Anhydrite"]["Poisson"], self.results_plot["Mudstone"]["Poisson"],
                           self.results_plot["Marl"]["Poisson"]]
                data_phi = [np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100,
                            np.array(self.results_plot["Marl"]["phi"]) * 100]
                data_gr = [self.results_plot["Limestone"]["GR"], self.results_plot["Dolomite"]["GR"],
                           self.results_plot["Anhydrite"]["GR"], self.results_plot["Mudstone"]["GR"],
                           self.results_plot["Marl"]["GR"]]
                data_pe = [self.results_plot["Limestone"]["PE"], self.results_plot["Dolomite"]["PE"],
                           self.results_plot["Anhydrite"]["PE"], self.results_plot["Mudstone"]["PE"],
                           self.results_plot["Marl"]["PE"]]
                #
                colors = ["skyblue", "lightcyan", "orchid", "olivedrab", "moccasin"]
                litho_list = ["Limestone", "Dolomite", "Anhydrite", "Mudstone", "Marl"]
                #
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_phi, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Porosity $\\phi$ (%)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            #
            self.create_3x3_scatter(parent=self.parent_subsurface, data_x=data_x_rho, data=data, row_id=2, column_id=9,
                                    n_rows=45, n_columns=9, colors=colors, labels=labels,
                                    xlabel="Densitiy $\\varrho$ (g/ccm)", lithos=litho_list)
        #
        elif var_rb.get() == 3: # Elements
            #
            try:
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except:
                pass
            #
            if self.name == "Random":
                if "Sandstone" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Sandstone":
                    self.list_elements_0 = self.list_elements_sst
                    self.elements_0 = self.elements_sst
                elif "Shale" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Shale":
                    self.list_elements_0 = self.list_elements_sh
                    self.elements_0 = self.elements_sh
            elif self.name == "Zechstein":
                if "Mudstone" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Mudstone":
                    self.list_elements_0 = list(self.data_chemistry["Mudstone"].keys())
                    self.elements_0 = self.data_chemistry["Mudstone"]
                elif "Anhydrite" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Anhydrite":
                    self.list_elements_0 = list(self.data_chemistry["Anhydrite"].keys())
                    self.elements_0 = self.data_chemistry["Anhydrite"]
                elif "Rock Salt" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Rock Salt":
                    self.list_elements_0 = list(self.data_chemistry["Rock Salt"].keys())
                    self.elements_0 = self.data_chemistry["Rock Salt"]
                elif "Potash" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Potash":
                    self.list_elements_0 = list(self.data_chemistry["Potash"].keys())
                    self.elements_0 = self.data_chemistry["Potash"]
                elif "Dolomite" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Dolomite":
                    self.list_elements_0 = list(self.data_chemistry["Dolomite"].keys())
                    self.elements_0 = self.data_chemistry["Dolomite"]
                elif "Limestone" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Limestone":
                    self.list_elements_0 = list(self.data_chemistry["Limestone"].keys())
                    self.elements_0 = self.data_chemistry["Limestone"]
                elif "Kupferschiefer" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Kupferschiefer":
                    self.list_elements_0 = list(self.data_chemistry["Kupferschiefer"].keys())
                    self.elements_0 = self.data_chemistry["Kupferschiefer"]
            #
            lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                   fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
            for index, element in enumerate(self.list_elements_0, start=0):
                if element not in ["U"]:
                    lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
            #
            for index, element in enumerate(self.list_elements_0, start=10):
                if element not in ["U"]:
                    entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(self.elements_0[element]), 3))
                    entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(self.elements_0[element]), 3))
                    entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(self.elements_0[element]), 3))
                    entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(self.elements_0[element], ddof=1), 3))
                else:
                    ppm_amounts = np.array(self.elements_0[element])*10000
                    entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(ppm_amounts), 3))
                    entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(ppm_amounts), 3))
                    entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(ppm_amounts), 3))
                    entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                #
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
        #
        elif var_rb.get() == 4: # Minerals
            #
            try:
                for lbl in self.lbl_w["chemistry"]:
                    lbl.grid_forget()
                for entr in self.entr_w["chemistry"]:
                    entr.grid_forget()
                self.lbl_w["chemistry"].clear()
                self.entr_w["chemistry"].clear()
            except:
                pass
            #
            lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                   fg="black").create_label(text="Mineralogical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w["chemistry"].append(lbl)
            for index, mineral in enumerate(self.list_minerals_0, start=0):
                if mineral not in ["Urn"]:
                    lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(mineral), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(mineral)+" (ppm)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
            #
            for index, mineral in enumerate(self.list_minerals_0, start=10):
                if mineral not in ["Urn"]:
                    entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(self.minerals_0[mineral]), 3))
                    entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(self.minerals_0[mineral]), 3))
                    entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(self.minerals_0[mineral]), 3))
                    entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(self.minerals_0[mineral], ddof=1), 3))
                else:
                    ppm_amounts = np.array(self.minerals_0[mineral])*10000
                    entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                 var_entr_set=round(np.min(ppm_amounts), 3))
                    entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                 var_entr_set=round(np.max(ppm_amounts), 3))
                    entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                   fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                  var_entr_set=round(np.mean(ppm_amounts), 3))
                    entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                  fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                 var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                #
                self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                #
        elif var_rb.get() == 5:
            if self.current_rb_lith.get() == var_rb.get():
                pass
            else:
                #
                try:
                    for lbl in self.lbl_w["chemistry"]:
                        lbl.grid_forget()
                    for entr in self.entr_w["chemistry"]:
                        entr.grid_forget()
                    self.lbl_w["chemistry"].clear()
                    self.entr_w["chemistry"].clear()
                except:
                    pass
                #
                self.current_rb_lith.set(var_rb.get())
                if self.list_rocks_short[0] == "Sandstone":
                    self.list_elements_0 = self.list_elements_sst
                    self.list_minerals_0 = self.list_minerals_sst
                    self.results_0 = self.results_sst
                    self.elements_0 = self.elements_sst
                    self.minerals_0 = self.minerals_sst
                elif self.list_rocks_short[0] == "Shale":
                    self.list_elements_0 = self.list_elements_sh
                    self.list_minerals_0 = self.list_minerals_sh
                    self.results_0 = self.results_sh
                    self.elements_0 = self.elements_sh
                    self.minerals_0 = self.minerals_sh
                #
                ## Entry Table
                for i in range(10):
                    if i == 7:
                        entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]*100), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]*100), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]*100), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i]*100, ddof=1), 3))
                    else:
                        entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i], ddof=1), 3))
                    #
                    self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
                #
                self.var_rb_geochem.set(3)
                lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                       fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
                for index, element in enumerate(self.list_elements_0, start=0):
                    if element not in ["U"]:
                        lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(element), relief=tk.RAISED)
                    else:
                        lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                    self.lbl_w["chemistry"].append(lbl)
                #
                for index, element in enumerate(self.list_elements_0, start=10):
                    if element not in ["U"]:
                        entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                     var_entr_set=round(np.min(self.elements_0[element]), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                     var_entr_set=round(np.max(self.elements_0[element]), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                      var_entr_set=round(np.mean(self.elements_0[element]), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                     var_entr_set=round(np.std(self.elements_0[element], ddof=1), 3))
                    else:
                        ppm_amounts = np.array(self.elements_0[element])*10000
                        entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                     var_entr_set=round(np.min(ppm_amounts), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                     var_entr_set=round(np.max(ppm_amounts), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                      var_entr_set=round(np.mean(ppm_amounts), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                     var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                    #
                    self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                    #
        elif var_rb.get() == 6:
            if self.current_rb_lith.get() == var_rb.get():
                pass
            else:
                #
                try:
                    for lbl in self.lbl_w["chemistry"]:
                        lbl.grid_forget()
                    for entr in self.entr_w["chemistry"]:
                        entr.grid_forget()
                    self.lbl_w["chemistry"].clear()
                    self.entr_w["chemistry"].clear()
                except:
                    pass
                #
                self.current_rb_lith.set(var_rb.get())
                if self.name == "Random":
                    if self.list_rocks_short[1] == "Sandstone":
                        self.list_elements_0 = self.list_elements_sst
                        self.list_minerals_0 = self.list_minerals_sst
                        self.results_0 = self.results_sst
                        self.elements_0 = self.elements_sst
                        self.minerals_0 = self.minerals_sst
                    elif self.list_rocks_short[1] == "Shale":
                        self.list_elements_0 = self.list_elements_sh
                        self.list_minerals_0 = self.list_minerals_sh
                        self.results_0 = self.results_sh
                        self.elements_0 = self.elements_sh
                        self.minerals_0 = self.minerals_sh
                elif self.name == "Zechstein":
                    if self.list_rocks_short[1] == "Mudstone":
                        self.list_elements_0 = list(self.data_chemistry["Mudstone"].keys())
                        self.list_minerals_0 = list(self.data_mineralogy["Mudstone"].keys())
                        self.results_0 = self.results_units["Mudstone"]
                        self.elements_0 = self.data_chemistry["Mudstone"]
                        self.minerals_0 = self.data_mineralogy["Mudstone"]
                    elif self.list_rocks_short[1] == "Anhydrite":
                        self.list_elements_0 = list(self.data_chemistry["Anhydrite"].keys())
                        self.list_minerals_0 = list(self.data_mineralogy["Anhydrite"].keys())
                        self.results_0 = self.results_units["Anhydrite"]
                        self.elements_0 = self.data_chemistry["Anhydrite"]
                        self.minerals_0 = self.data_mineralogy["Anhydrite"]
                #
                ## Entry Table
                for i in range(10):
                    if i == 7:
                        entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]*100), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]*100), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]*100), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i]*100, ddof=1), 3))
                    else:
                        entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i], ddof=1), 3))
                    #
                    self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
                #
                self.var_rb_geochem.set(3)
                lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                       fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
                for index, element in enumerate(self.list_elements_0, start=0):
                    if element not in ["U"]:
                        lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(element), relief=tk.RAISED)
                    else:
                        lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                    self.lbl_w["chemistry"].append(lbl)
                #
                for index, element in enumerate(self.list_elements_0, start=10):
                    if element not in ["U"]:
                        entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                     var_entr_set=round(np.min(self.elements_0[element]), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                     var_entr_set=round(np.max(self.elements_0[element]), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                      var_entr_set=round(np.mean(self.elements_0[element]), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                     var_entr_set=round(np.std(self.elements_0[element], ddof=1), 3))
                    else:
                        ppm_amounts = np.array(self.elements_0[element])*10000
                        entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                     var_entr_set=round(np.min(ppm_amounts), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                     var_entr_set=round(np.max(ppm_amounts), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                      var_entr_set=round(np.mean(ppm_amounts), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                     var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                    #
                    self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                    #
        elif var_rb.get() == 7:
            if self.current_rb_lith.get() == var_rb.get():
                pass
            else:
                #
                try:
                    for lbl in self.lbl_w["chemistry"]:
                        lbl.grid_forget()
                    for entr in self.entr_w["chemistry"]:
                        entr.grid_forget()
                    self.lbl_w["chemistry"].clear()
                    self.entr_w["chemistry"].clear()
                except:
                    pass
                #
                self.current_rb_lith.set(var_rb.get())
                if self.name == "Random":
                    if self.list_rocks_short[-1] == "Granite":
                        self.list_elements_0 = self.list_elements_granite
                        self.list_minerals_0 = self.list_minerals_granite
                        self.results_0 = self.results_granite
                        self.elements_0 = self.elements_granite
                        self.minerals_0 = self.minerals_granite
                    elif self.list_rocks_short[-1] == "Gabbro":
                        self.list_elements_0 = self.list_elements_gabbro
                        self.list_minerals_0 = self.list_minerals_gabbro
                        self.results_0 = self.results_gabbro
                        self.elements_0 = self.elements_gabbro
                        self.minerals_0 = self.minerals_gabbro
                    elif self.list_rocks_short[-1] == "Diorite":
                        self.list_elements_0 = self.list_elements_diorite
                        self.list_minerals_0 = self.list_minerals_diorite
                        self.results_0 = self.results_diorite
                        self.elements_0 = self.elements_diorite
                        self.minerals_0 = self.minerals_diorite
                #
                ## Entry Table
                for i in range(10):
                    if i == 7:
                        entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]*100), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]*100), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]*100), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i]*100, ddof=1), 3))
                    else:
                        entr_min = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=4, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[i], var_entr_set=round(np.min(self.results_0[i]), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=5, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[i], var_entr_set=round(np.max(self.results_0[i]), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=6, n_rows=2, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[i], var_entr_set=round(np.mean(self.results_0[i]), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=4+2*i, column_id=7, n_rows=2, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[i], var_entr_set=round(np.std(self.results_0[i], ddof=1), 3))
                    #
                    self.entr_w["physics"].extend([entr_min, entr_max, entr_mean, entr_std])
                #
                self.var_rb_geochem.set(3)
                lbl = SE(parent=self.parent_subsurface, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                       fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
                self.lbl_w["chemistry"].append(lbl)
                for index, element in enumerate(self.list_elements_0, start=0):
                    if element not in ["U"]:
                        lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(element), relief=tk.RAISED)
                    else:
                        lbl = SE(parent=self.parent_subsurface, row_id=25+index, column_id=3, bg=self.color_bg,
                                 fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                    self.lbl_w["chemistry"].append(lbl)
                #
                for index, element in enumerate(self.list_elements_0, start=10):
                    if element not in ["U"]:
                        entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                     var_entr_set=round(np.min(self.elements_0[element]), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                     var_entr_set=round(np.max(self.elements_0[element]), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                      var_entr_set=round(np.mean(self.elements_0[element]), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                     var_entr_set=round(np.std(self.elements_0[element], ddof=1), 3))
                    else:
                        ppm_amounts = np.array(self.elements_0[element])*10000
                        entr_min = SE(parent=self.parent_subsurface, row_id=15+index, column_id=4, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_min[index],
                                                                     var_entr_set=round(np.min(ppm_amounts), 3))
                        entr_max = SE(parent=self.parent_subsurface, row_id=15+index, column_id=5, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_max[index],
                                                                     var_entr_set=round(np.max(ppm_amounts), 3))
                        entr_mean = SE(parent=self.parent_subsurface, row_id=15+index, column_id=6, bg=self.color_bg,
                                       fg=self.color_fg).create_entry(var_entr=self.entr_list_mean[index],
                                                                      var_entr_set=round(np.mean(ppm_amounts), 3))
                        entr_std = SE(parent=self.parent_subsurface, row_id=15+index, column_id=7, bg=self.color_bg,
                                      fg=self.color_fg).create_entry(var_entr=self.entr_list_std[index],
                                                                     var_entr_set=round(np.std(ppm_amounts, ddof=1), 3))
                    #
                    self.entr_w["chemistry"].extend([entr_min, entr_max, entr_mean, entr_std])
                    #
        elif var_rb.get() == 8: # Density
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            if self.name == "Random":
                data_x_rho = [np.array(self.results_plot["Sandstone"]["rho"]) / 1000,
                              np.array(self.results_plot["Shale"]["rho"]) / 1000,
                              np.array(self.results_plot[self.list_rocks_short[-1]]["rho"]) / 1000]
                #
                data_vP = [self.results_plot["Sandstone"]["vP"], self.results_plot["Shale"]["vP"],
                           self.results_plot[self.list_rocks_short[-1]]["vP"]]
                data_vS = [self.results_plot["Sandstone"]["vS"], self.results_plot["Shale"]["vS"],
                           self.results_plot[self.list_rocks_short[-1]]["vS"]]
                data_vPvS = [self.results_plot["Sandstone"]["vPvS"], self.results_plot["Shale"]["vPvS"],
                             self.results_plot[self.list_rocks_short[-1]]["vPvS"]]
                data_K = [self.results_plot["Sandstone"]["K"], self.results_plot["Shale"]["K"],
                          self.results_plot[self.list_rocks_short[-1]]["K"]]
                data_G = [self.results_plot["Sandstone"]["G"], self.results_plot["Shale"]["G"],
                          self.results_plot[self.list_rocks_short[-1]]["G"]]
                data_nu = [self.results_plot["Sandstone"]["Poisson"], self.results_plot["Shale"]["Poisson"],
                           self.results_plot[self.list_rocks_short[-1]]["Poisson"]]
                data_phi = [np.array(self.results_plot["Sandstone"]["phi"]) * 100,
                            np.array(self.results_plot["Shale"]["phi"]) * 100,
                            np.array(self.results_plot[self.list_rocks_short[-1]]["phi"]) * 100]
                data_gr = [self.results_plot["Sandstone"]["GR"], self.results_plot["Shale"]["GR"],
                           self.results_plot[self.list_rocks_short[-1]]["GR"]]
                data_pe = [self.results_plot["Sandstone"]["PE"], self.results_plot["Shale"]["PE"],
                           self.results_plot[self.list_rocks_short[-1]]["PE"]]
                #
                colors = ["tan", "olivedrab", "darkorange"]
                litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
            elif self.name == "Zechstein":
                data_x_rho = [np.array(self.results_plot["Anhydrite"]["rho"]) / 1000,
                              np.array(self.results_plot["Dolomite"]["rho"]) / 1000,
                              np.array(self.results_plot["Rock Salt"]["rho"]) / 1000,
                              np.array(self.results_plot["Limestone"]["rho"]) / 1000,
                              np.array(self.results_plot["Kupferschiefer"]["rho"]) / 1000,
                              np.array(self.results_plot["Potash"]["rho"]) / 1000,
                              np.array(self.results_plot["Mudstone"]["rho"]) / 1000]
                #
                data_vP = [self.results_plot["Anhydrite"]["vP"], self.results_plot["Dolomite"]["vP"],
                           self.results_plot["Rock Salt"]["vP"], self.results_plot["Limestone"]["vP"],
                           self.results_plot["Kupferschiefer"]["vP"], self.results_plot["Potash"]["vP"],
                           self.results_plot["Mudstone"]["vP"]]
                data_vS = [self.results_plot["Anhydrite"]["vS"], self.results_plot["Dolomite"]["vS"],
                           self.results_plot["Rock Salt"]["vS"], self.results_plot["Limestone"]["vS"],
                           self.results_plot["Kupferschiefer"]["vS"], self.results_plot["Potash"]["vS"],
                           self.results_plot["Mudstone"]["vS"]]
                data_vPvS = [self.results_plot["Anhydrite"]["vPvS"], self.results_plot["Dolomite"]["vPvS"],
                             self.results_plot["Rock Salt"]["vPvS"], self.results_plot["Limestone"]["vPvS"],
                             self.results_plot["Kupferschiefer"]["vPvS"], self.results_plot["Potash"]["vPvS"],
                             self.results_plot["Mudstone"]["vPvS"]]
                data_K = [self.results_plot["Anhydrite"]["K"], self.results_plot["Dolomite"]["K"],
                          self.results_plot["Rock Salt"]["K"], self.results_plot["Limestone"]["K"],
                          self.results_plot["Kupferschiefer"]["K"], self.results_plot["Potash"]["K"],
                          self.results_plot["Mudstone"]["K"]]
                data_G = [self.results_plot["Anhydrite"]["G"], self.results_plot["Dolomite"]["G"],
                          self.results_plot["Rock Salt"]["G"], self.results_plot["Limestone"]["G"],
                          self.results_plot["Kupferschiefer"]["G"], self.results_plot["Potash"]["G"],
                          self.results_plot["Mudstone"]["G"]]
                data_nu = [self.results_plot["Anhydrite"]["Poisson"], self.results_plot["Dolomite"]["Poisson"],
                           self.results_plot["Rock Salt"]["Poisson"], self.results_plot["Limestone"]["Poisson"],
                           self.results_plot["Kupferschiefer"]["Poisson"], self.results_plot["Potash"]["Poisson"],
                           self.results_plot["Mudstone"]["Poisson"]]
                data_phi = [np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                            np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                            np.array(self.results_plot["Rock Salt"]["phi"]) * 100,
                            np.array(self.results_plot["Limestone"]["phi"]) * 100,
                            np.array(self.results_plot["Kupferschiefer"]["phi"]) * 100,
                            np.array(self.results_plot["Potash"]["phi"]) * 100,
                            np.array(self.results_plot["Mudstone"]["phi"]) * 100]
                data_gr = [self.results_plot["Anhydrite"]["GR"], self.results_plot["Dolomite"]["GR"],
                           self.results_plot["Rock Salt"]["GR"], self.results_plot["Limestone"]["GR"],
                           self.results_plot["Kupferschiefer"]["GR"], self.results_plot["Potash"]["GR"],
                           self.results_plot["Mudstone"]["GR"]]
                data_pe = [self.results_plot["Anhydrite"]["PE"], self.results_plot["Dolomite"]["PE"],
                           self.results_plot["Rock Salt"]["PE"], self.results_plot["Limestone"]["PE"],
                           self.results_plot["Kupferschiefer"]["PE"], self.results_plot["Potash"]["PE"],
                           self.results_plot["Mudstone"]["PE"]]
                #
                colors = ["orchid", "lightcyan", "lavender", "skyblue", "gray", "yellowgreen", "olivedrab"]
                litho_list = ["Anhydrite", "Dolomite", "Rock Salt", "Limestone", "Kupferschiefer", "Potash", "Mudstone"]
                #
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_phi, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Porosity $\\phi$ (%)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            #
            self.create_3x3_scatter(parent=self.parent_subsurface, data_x=data_x_rho, data=data, row_id=2, column_id=9,
                                    n_rows=45, n_columns=9, colors=colors, labels=labels,
                                    xlabel="Densitiy $\\varrho$ (g/ccm)", lithos=litho_list)
        #
        elif var_rb.get() == 9: # Porosity
            try:
                self.fig.clf()
                self.ax.cla()
                self.canvas.get_tk_widget().pack_forget()
            except AttributeError:
                pass
            try:
                if self.canvas:
                    self.canvas.destroy()
            except AttributeError:
                pass
            #
            if self.name == "Random":
                data_x = [np.array(self.results_plot["Sandstone"]["phi"]) * 100,
                          np.array(self.results_plot["Shale"]["phi"]) * 100,
                          np.array(self.results_plot[self.list_rocks_short[-1]]["phi"]) * 100]
                #
                data_vP = [self.results_plot["Sandstone"]["vP"], self.results_plot["Shale"]["vP"],
                           self.results_plot[self.list_rocks_short[-1]]["vP"]]
                data_vS = [self.results_plot["Sandstone"]["vS"], self.results_plot["Shale"]["vS"],
                           self.results_plot[self.list_rocks_short[-1]]["vS"]]
                data_vPvS = [self.results_plot["Sandstone"]["vPvS"], self.results_plot["Shale"]["vPvS"],
                             self.results_plot[self.list_rocks_short[-1]]["vPvS"]]
                data_K = [self.results_plot["Sandstone"]["K"], self.results_plot["Shale"]["K"],
                          self.results_plot[self.list_rocks_short[-1]]["K"]]
                data_G = [self.results_plot["Sandstone"]["G"], self.results_plot["Shale"]["G"],
                          self.results_plot[self.list_rocks_short[-1]]["G"]]
                data_nu = [self.results_plot["Sandstone"]["Poisson"], self.results_plot["Shale"]["Poisson"],
                           self.results_plot[self.list_rocks_short[-1]]["Poisson"]]
                data_rho = [np.array(self.results_plot["Sandstone"]["rho"]) / 1000,
                            np.array(self.results_plot["Shale"]["rho"]) / 1000,
                            np.array(self.results_plot[self.list_rocks_short[-1]]["rho"]) / 1000]
                data_gr = [self.results_plot["Sandstone"]["GR"], self.results_plot["Shale"]["GR"],
                           self.results_plot[self.list_rocks_short[-1]]["GR"]]
                data_pe = [self.results_plot["Sandstone"]["PE"], self.results_plot["Shale"]["PE"],
                           self.results_plot[self.list_rocks_short[-1]]["PE"]]
                #
                colors = ["tan", "olivedrab", "darkorange"]
                litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
            elif self.name == "Zechstein":
                data_x = [np.array(self.results_plot["Anhydrite"]["phi"]) * 100,
                          np.array(self.results_plot["Dolomite"]["phi"]) * 100,
                          np.array(self.results_plot["Rock Salt"]["phi"]) * 100,
                          np.array(self.results_plot["Limestone"]["phi"]) * 100,
                          np.array(self.results_plot["Kupferschiefer"]["phi"]) * 100,
                          np.array(self.results_plot["Potash"]["phi"]) * 100,
                          np.array(self.results_plot["Mudstone"]["phi"]) * 100]
                #
                data_vP = [self.results_plot["Anhydrite"]["vP"], self.results_plot["Dolomite"]["vP"],
                           self.results_plot["Rock Salt"]["vP"], self.results_plot["Limestone"]["vP"],
                           self.results_plot["Kupferschiefer"]["vP"], self.results_plot["Potash"]["vP"],
                           self.results_plot["Mudstone"]["vP"]]
                data_vS = [self.results_plot["Anhydrite"]["vS"], self.results_plot["Dolomite"]["vS"],
                           self.results_plot["Rock Salt"]["vS"], self.results_plot["Limestone"]["vS"],
                           self.results_plot["Kupferschiefer"]["vS"], self.results_plot["Potash"]["vS"],
                           self.results_plot["Mudstone"]["vS"]]
                data_vPvS = [self.results_plot["Anhydrite"]["vPvS"], self.results_plot["Dolomite"]["vPvS"],
                             self.results_plot["Rock Salt"]["vPvS"], self.results_plot["Limestone"]["vPvS"],
                             self.results_plot["Kupferschiefer"]["vPvS"], self.results_plot["Potash"]["vPvS"],
                             self.results_plot["Mudstone"]["vPvS"]]
                data_K = [self.results_plot["Anhydrite"]["K"], self.results_plot["Dolomite"]["K"],
                          self.results_plot["Rock Salt"]["K"], self.results_plot["Limestone"]["K"],
                          self.results_plot["Kupferschiefer"]["K"], self.results_plot["Potash"]["K"],
                          self.results_plot["Mudstone"]["K"]]
                data_G = [self.results_plot["Anhydrite"]["G"], self.results_plot["Dolomite"]["G"],
                          self.results_plot["Rock Salt"]["G"], self.results_plot["Limestone"]["G"],
                          self.results_plot["Kupferschiefer"]["G"], self.results_plot["Potash"]["G"],
                          self.results_plot["Mudstone"]["G"]]
                data_nu = [self.results_plot["Anhydrite"]["Poisson"], self.results_plot["Dolomite"]["Poisson"],
                           self.results_plot["Rock Salt"]["Poisson"], self.results_plot["Limestone"]["Poisson"],
                           self.results_plot["Kupferschiefer"]["Poisson"], self.results_plot["Potash"]["Poisson"],
                           self.results_plot["Mudstone"]["Poisson"]]
                data_rho = [np.array(self.results_plot["Anhydrite"]["rho"]) / 1000,
                              np.array(self.results_plot["Dolomite"]["rho"]) / 1000,
                              np.array(self.results_plot["Rock Salt"]["rho"]) / 1000,
                              np.array(self.results_plot["Limestone"]["rho"]) / 1000,
                              np.array(self.results_plot["Kupferschiefer"]["rho"]) / 1000,
                              np.array(self.results_plot["Potash"]["rho"]) / 1000,
                              np.array(self.results_plot["Mudstone"]["rho"]) / 1000]
                data_gr = [self.results_plot["Anhydrite"]["GR"], self.results_plot["Dolomite"]["GR"],
                           self.results_plot["Rock Salt"]["GR"], self.results_plot["Limestone"]["GR"],
                           self.results_plot["Kupferschiefer"]["GR"], self.results_plot["Potash"]["GR"],
                           self.results_plot["Mudstone"]["GR"]]
                data_pe = [self.results_plot["Anhydrite"]["PE"], self.results_plot["Dolomite"]["PE"],
                           self.results_plot["Rock Salt"]["PE"], self.results_plot["Limestone"]["PE"],
                           self.results_plot["Kupferschiefer"]["PE"], self.results_plot["Potash"]["PE"],
                           self.results_plot["Mudstone"]["PE"]]
                #
                colors = ["orchid", "lightcyan", "lavender", "skyblue", "gray", "yellowgreen", "olivedrab"]
                litho_list = ["Anhydrite", "Dolomite", "Rock Salt", "Limestone", "Kupferschiefer", "Potash", "Mudstone"]
                #
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_rho, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Density $\\varrho$ (g/ccm)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            #
            self.create_3x3_scatter(parent=self.parent_subsurface, data_x=data_x, data=data, row_id=2, column_id=9,
                                    n_rows=45, n_columns=9, colors=colors, labels=labels, xlabel="Porosity $\\phi$ (%)",
                                    lithos=litho_list)
        #
    def split_thickness(self, thickness, n_units):
        list_thickness = np.random.multinomial(thickness, np.ones(n_units)/n_units, size=1)[0]
        #
        return list_thickness
    #
    def create_3x3_scatter(self, parent, data_x, data, row_id, column_id, n_rows, n_columns, colors, labels, xlabel, lithos):
        #
        self.canvas = None
        self.fig, self.ax = plt.subplots(ncols=3, nrows=3, figsize=(9, 9), facecolor="#E9ECED")
        #
        for i in range(3):
            for j in range(3):
                for k in range(len(lithos)):
                    self.ax[i][j].scatter(data_x[k], data[i][j][k], color=colors[k], edgecolor="black", alpha=1.0,
                                          label=lithos[k])
                #
                self.ax[i][j].set_xlabel(xlabel, fontsize="small")
                self.ax[i][j].set_ylabel(labels[i][j], labelpad=0.5, fontsize="small")
                self.ax[i][j].grid(True)
                self.ax[i][j].set_axisbelow(True)
        #
        self.ax[0][0].legend(bbox_to_anchor=(0.0, 1.025, 3.5, 0.15), loc=3, ncol=len(lithos), mode="expand",
                             borderaxespad=0, fontsize="x-small")
        self.fig.tight_layout()
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def create_3x3_histo(self, parent, data, row_id, column_id, n_rows, n_columns, colors, labels, lithos):
        #
        self.canvas = None
        self.fig, self.ax = plt.subplots(ncols=3, nrows=3, figsize=(9, 9), facecolor="#E9ECED")
        #
        for i in range(3):
            for j in range(3):
                self.ax[i][j].hist(data[i][j], bins=16, color=colors, histtype="step", fill=True, linewidth=2,
                                   alpha=0.25)
                self.ax[i][j].hist(data[i][j], bins=16, color=colors, histtype="step", linewidth=2, label=lithos)
                #
                self.ax[i][j].set_xlabel(labels[i][j], fontsize="small")
                self.ax[i][j].set_ylabel("Frequency", labelpad=0.5, fontsize="small")
                self.ax[i][j].grid(True)
                self.ax[i][j].set_axisbelow(True)
        self.ax[0][0].legend(bbox_to_anchor=(0.0, 1.025, 3.5, 0.15), loc=3, ncol=len(lithos), mode="expand",
                             borderaxespad=0, fontsize="x-small")
        self.fig.tight_layout()
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def create_well_log_plot(self, parent, data_x, data_y, row_id, column_id, n_rows, n_columns, unit_sections=None):
        #
        self.canvas = None
        max_thickness = max(data_y)
        if max_thickness <= 100:
            step_depth = 10
        elif 100 < max_thickness <= 500:
            step_depth = 50
        elif 500 < max_thickness <= 1500:
            step_depth = 100
        elif max_thickness > 1500:
            step_depth = 200
        self.fig, (self.ax1, self.ax2, self.ax3, self.ax4, self.ax5) = plt.subplots(
            1, 5, sharey="row", gridspec_kw={"wspace": 0.25}, figsize=(12, 24), facecolor="#E9ECED")
        self.fig.subplots_adjust(wspace=0.25)
        # 1
        self.ax1.plot(self.results_sorted["GR"], self.results_sorted["Top"], color="#00549F", linewidth=2)
        self.ax1.set_xlabel("GR [API]")
        self.ax1.set_ylabel("Depth [m]")
        if max(self.results_sorted["GR"]) > 250:
            self.ax1.set_xscale("log")
        else:
            self.ax1.set_xlim(-1, max(self.results_sorted["GR"]))
            self.ax1.set_xticks(np.arange(0, max(self.results_sorted["GR"])+50, 50))
        self.ax1.set_ylim(0, max_thickness)
        self.ax1.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax1.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 2
        vP_edit = np.array(self.results_sorted["vP"])/1000
        vS_edit = np.array(self.results_sorted["vS"])/1000
        self.ax2.plot(vP_edit, self.results_sorted["Top"], color="#00549F", linewidth=2)
        self.ax2.set_xlabel("$v_P$ [km/s]")
        self.ax2.set_xlim(0, max(vP_edit))
        self.ax2.set_xticks(np.arange(0, max(vP_edit)+1.0, 2.0))
        self.ax2.xaxis.label.set_color("#00549F")
        self.ax2.set_ylim(0, max_thickness)
        self.ax2.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax2.grid(color="grey", linestyle="dashed")
        self.ax2_2 = self.ax2.twiny()
        self.ax2_2.plot(vS_edit, self.results_sorted["Top"], color="#CC071E", linewidth=2)
        self.ax2_2.set_xlabel("$v_S$ [km/s]")
        self.ax2_2.set_xlim(0, max(vP_edit))
        self.ax2_2.set_xticks(np.arange(0, max(vP_edit)+1.0, 2.0))
        self.ax2_2.minorticks_on()
        self.ax2_2.xaxis.label.set_color("#CC071E")
        self.ax2_2.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 3
        phi_edit = np.array(self.results_sorted["phi"])*100
        self.ax3.plot(np.array(self.results_sorted["rho"])/1000, self.results_sorted["Top"], color="#57AB27", linewidth=2)
        self.ax3.set_xlabel("$\\varrho$ [g/cm$^3$]")
        self.ax3.set_xlim(1.6, 3.2)
        self.ax3.set_xticks(np.around(np.linspace(1.6, 3.2, 4, endpoint=True), decimals=1))
        self.ax3.xaxis.label.set_color("#57AB27")
        self.ax3.set_ylim(0, max_thickness)
        self.ax3.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax3.grid(color="grey", linestyle="dashed")
        self.ax3_2 = self.ax3.twiny()
        self.ax3_2.plot(phi_edit, self.results_sorted["Top"], color="#00549F", linewidth=2)
        self.ax3_2.set_xlabel("$\\varphi$ [1]")
        self.ax3_2.set_xlim(60, 0)
        self.ax3_2.set_xticks(np.around(np.linspace(60, 0, 6, endpoint=True), decimals=0))
        self.ax3_2.minorticks_on()
        self.ax3_2.xaxis.label.set_color("#00549F")
        self.ax3_2.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 4
        self.ax4.plot(self.results_sorted["PE"], self.results_sorted["Top"], color="#00549F", linewidth=2)
        self.ax4.set_xlabel("PE [barns/electron]")
        self.ax4.set_xlim(min(self.results_sorted["PE"]), max(self.results_sorted["PE"]))
        self.ax4.set_xscale("log")
        self.ax4.set_ylim(0, max_thickness)
        self.ax4.set_yticks(np.arange(0, max_thickness+step_depth, step_depth))
        self.ax4.grid(color="grey", linestyle="dashed", which="both")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 5
        if unit_sections == None:
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
            n_units = len(unit_sections)
            legend_lithology = []
            for key, value in unit_sections.items():
                legend_lithology.append(
                    mpatches.Patch(facecolor=value["Color"], edgecolor="black", hatch="", label=key))
            for key, value in unit_sections.items():
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
        #plt.tight_layout()
        #
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=row_id, column=column_id, rowspan=n_rows, columnspan=n_columns,
                                         sticky="nesw")
    #
    def extract_data(self, data):
        print(data)
    #
if __name__ == "__main__":
    #
    # SYSTEM CHECK
    name_os = os.name
    name_platform = sys.platform
    if name_platform == "darwin":
        color_background = "#FFFFFF"
        button_color_bg = "lightgrey"
    elif name_platform == "linux":
        color_background = "#DADBD9"
    elif name_platform == "win32":
        color_background = "#F0F0EE"
    #
    root = tk.Tk()
    GebPyGUI(root)
    root.mainloop()