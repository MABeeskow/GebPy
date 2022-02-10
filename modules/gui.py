#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gui.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		10.02.2022

#-----------------------------------------------

## MODULES
import os, sys
import tkinter as tk
from modules.gui_elements import SimpleElements as SE
from modules.sulfates import Sulfates
from modules.oxides import Oxides
from modules.sulfides import Sulfides
from modules.carbonates import Carbonates
from modules.halogenes import Halogenes
from modules.silicates import Tectosilicates, Phyllosilicates, Nesosilicates, Sorosilicates, Inosilicates
from modules.phosphates import Phosphates
from modules.minerals import feldspars
from modules.siliciclastics import sandstone, shale
from modules.carbonates import limestone, dolomite, CustomCarbonates
from modules.sequences import SedimentaryBasin
from modules.ore import Ores
from modules.igneous import Plutonic
from modules.evaporites import Evaporites
from modules.sequences import DataProcessing as DP
import numpy as np
import random as rd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import NullFormatter
import matplotlib.patches as mpatches
import time

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
        self.parent.geometry("1700x960")
        self.parent.resizable(False, False)
        self.parent["bg"] = self.color_bg
        #
        self.lbl_w = {}
        self.lbl_w["physics"] = []
        self.lbl_w["chemistry"] = []
        self.entr_w = {}
        self.entr_w["physics"] = []
        self.entr_w["chemistry"] = []
        self.gui_elements = []
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
        self.parent.grid_columnconfigure(0, minsize=170)
        self.parent.grid_columnconfigure(1, minsize=100)
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
        ## Labels
        SE(parent=self.parent, row_id=6, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
           fg=self.color_fg_dark).create_label(text="Minerals", relief=tk.RAISED)
        SE(parent=self.parent, row_id=12, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
           fg=self.color_fg_dark).create_label(text="Rocks", relief=tk.RAISED)
        SE(parent=self.parent, row_id=20, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
           fg=self.color_fg_dark).create_label(text="Subsurface", relief=tk.RAISED)
        SE(parent=self.parent, row_id=0, column_id=3, n_rows=2, n_columns=5, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Statistics", relief=tk.RAISED)
        SE(parent=self.parent, row_id=2, column_id=3, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Parameter", relief=tk.RAISED)
        SE(parent=self.parent, row_id=2, column_id=4, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Minimum", relief=tk.RAISED)
        SE(parent=self.parent, row_id=2, column_id=5, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Maximum", relief=tk.RAISED)
        SE(parent=self.parent, row_id=2, column_id=6, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Mean", relief=tk.RAISED)
        SE(parent=self.parent, row_id=2, column_id=7, n_rows=2, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Standard\n Deviation", relief=tk.RAISED)
        SE(parent=self.parent, row_id=0, column_id=9, n_rows=2, n_columns=9, bg=self.color_bg,
           fg=self.color_fg_dark).create_label(text="Plots", relief=tk.RAISED)
        #
        ## Option Menu
        var_opt_0_0 = tk.StringVar()
        opt_list_0_0 = ["Oxides", "Sulfides", "Carbonates", "Halogenes", "Tectosilicates", "Phyllosilicates",
                        "Sulfates", "Nesosilicates", "Sorosilicates", "Inosilicates", "Phosphates"]
        opt_list_0_0.sort()
        self.opt_mingroup = SE(parent=self.parent, row_id=8, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
            var_opt=var_opt_0_0, var_opt_set="Select Mineral Group", opt_list=opt_list_0_0,
            command=lambda var_opt=var_opt_0_0: self.select_opt(var_opt))
        var_opt_1_0 = tk.StringVar()
        opt_list_1_0 = ["Siliciclastic Rocks", "Carbonate Rocks", "Igneous Rocks", "Metamorphic Rocks",
                        "Evaporite Rocks", "Ore Rocks"]
        self.opt_rocktype = SE(parent=self.parent, row_id=14, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
            var_opt=var_opt_1_0, var_opt_set="Select Rock Type", opt_list=opt_list_1_0,
            command=lambda var_opt=var_opt_1_0: self.select_opt(var_opt))
        var_opt_2_0 = tk.StringVar()
        opt_list_2_0 = ["Zechstein", "Buntsandstein"]
        self.opt_realseq = SE(parent=self.parent, row_id=22, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
            var_opt=var_opt_2_0, var_opt_set="Select Real Sequences", opt_list=opt_list_2_0,
            command=lambda var_opt=var_opt_2_0: self.select_opt(var_opt))
        #
        ## Button
        self.btn_randseq = SE(parent=self.parent, row_id=24, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                              fg=self.color_fg_dark).create_button(text="Create Random Sequence",
                                                                   command=lambda var_btn="random": self.pressed_button(var_btn))
        self.btn_custseq = SE(parent=self.parent, row_id=26, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                              fg=self.color_fg_dark).create_button(text="Create Custom Sequence")
        self.btn_custseq = SE(parent=self.parent, row_id=18, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                              fg=self.color_fg_dark).create_button(text="Create Custom Rock",
                                                                   command=lambda var_btn="custom rock": self.pressed_button(var_btn))
    #
    def pressed_button(self, var_btn):
        if var_btn == "random":
            Subsurface(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                       color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_btn, lbl_w=self.lbl_w,
                       entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_btn == "custom rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock="Custom", lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        else:
            Subsurface(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                       color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_btn, lbl_w=self.lbl_w,
                       entr_w=self.entr_w, gui_elements=self.gui_elements)
    #
    def select_opt(self, var_opt):
        # Minerals
        ## Oxides
        if var_opt in ["Quartz", "Magnetite", "Hematite", "Aluminium Spinels", "Ilmenite", "Cassiterite", "Chromite",
                       "Corundum", "Rutile", "Pyrolusite", "Magnesiochromite", "Zincochromite", "Chromium Spinels",
                       "Cuprospinel", "Jacobsite", "Magnesioferrite", "Trevorite", "Franklinite", "Ulvöspinel",
                       "Iron Spinels", "Uraninite", "Litharge", "Massicot", "Minium", "Plattnerite", "Scrutinyite",
                       "Zincite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        ## Sulfides
        elif var_opt in ["Pyrite", "Chalcopyrite", "Galena", "Acanthite", "Chalcocite", "Bornite", "Sphalerite",
                         "Pyrrhotite", "Millerite", "Pentlandite", "Covellite", "Cinnabar", "Realgar", "Orpiment",
                         "Stibnite", "Marcasite", "Molybdenite", "Fahlore"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        ## Carbonates
        elif var_opt in ["Calcite", "Dolomite", "Magnesite", "Halite", "Fluorite", "Sylvite", "Siderite",
                         "Rhodochrosite", "Aragonite", "Cerussite", "Ankerite", "Azurite", "Malachite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt in ["Barite", "Celestite", "Anglesite", "Anhydrite", "Hanksite", "Gypsum", "Alunite", "Jarosite",
                         "Chalcanthite", "Kieserite", "Scheelite", "Hexahydrite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt in ["Alkalifeldspar", "Plagioclase", "Scapolite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        # PHOSPHATES
        elif var_opt in ["Apatite-F", "Apatite-Cl", "Apatite-OH", "Apatite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt in ["Illite", "Kaolinite", "Montmorillonite", "Chamosite", "Clinochlore", "Pennantite", "Nimite",
                         "Chlorite", "Vermiculite", "Annite", "Phlogopite", "Eastonite", "Siderophyllite", "Biotite",
                         "Muscovite", "Glauconite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        # NESOSILICATES
        elif var_opt in ["Zircon", "Thorite", "Andalusite", "Kyanite", "Sillimanite", "Topaz", "Staurolite",
                         "Fayalite", "Forsterite", "Tephroite", "Calcio-Olivine", "Liebenbergite", "Olivine", "Pyrope",
                         "Almandine", "Grossular", "Andradite", "Uvarovite", "Aluminium Garnet", "Calcium Garnet"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        # SOROSILICATES
        elif var_opt in ["Epidote", "Zoisite", "Gehlenite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        # INOSILICATES
        elif var_opt in ["Enstatite", "Ferrosilite", "Diopside", "Jadeite", "Aegirine", "Spodumene", "Wollastonite",
                         "Tremolite", "Actinolite", "Glaucophane", "Augite", "Riebeckite", "Arfvedsonite",
                         "Calcium Amphibole", "Sodium Amphibole", "Mg-Fe Pyroxene", "Calcium Pyroxene"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w, gui_elements=self.gui_elements)
        # Rocks
        elif var_opt == "Sandstone":
            self.lbl_w, self.entr_w = Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)()
        elif var_opt == "Shale":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        #
        elif var_opt == "Limestone":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Dolomite Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Custom Carbonate Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        #
        elif var_opt == "Felsic Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Intermediate Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Granite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Gabbro":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Diorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Granodiorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Syenite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Tonalite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Monzonite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Quartzolite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Qz-rich Granitoid":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        #
        elif var_opt == "Rock Salt":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt == "Anhydrite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        elif var_opt in ["Kupferschiefer"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w, gui_elements=self.gui_elements)
        #
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
            except:
                pass
            var_opt_0_1 = tk.StringVar()
            opt_list_0_1 = ["Quartz", "Magnetite", "Hematite", "Aluminium Spinels", "Iron Spinels", "Chromium Spinels",
                            "Corundum", "Ilmenite", "Rutile", "Pyrolusite", "Cassiterite", "Chromite",
                            "Magnesiochromite", "Zincochromite", "Cuprospinel", "Jacobsite", "Magnesioferrite",
                            "Trevorite", "Franklinite", "Ulvöspinel", "Uraninite", "Litharge", "Massicot", "Minium",
                            "Plattnerite", "Scrutinyite", "Zincite"]
            opt_list_0_1.sort()
            self.opt_oxide = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_1, var_opt_set="Select Oxide Mineral", opt_list=opt_list_0_1, active_bg=self.color_accent_02,
                command=lambda var_opt=var_opt_0_1: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_2 = tk.StringVar()
            opt_list_0_2 = ["Pyrite", "Chalcopyrite", "Galena", "Acanthite", "Chalcocite", "Bornite", "Sphalerite",
                            "Pyrrhotite", "Millerite", "Pentlandite", "Covellite", "Cinnabar", "Realgar", "Orpiment",
                            "Stibnite", "Marcasite", "Molybdenite", "Fahlore"]
            opt_list_0_2.sort()
            self.opt_sulfide = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                  bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_2, var_opt_set="Select Sulfide Mineral", opt_list=opt_list_0_2,
                command=lambda var_opt=var_opt_0_2: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_3 = tk.StringVar()
            opt_list_0_3 = ["Calcite", "Dolomite", "Magnesite", "Siderite", "Rhodochrosite", "Aragonite", "Cerussite",
                            "Ankerite", "Azurite", "Malachite"]
            opt_list_0_3.sort()
            self.opt_carb = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                               bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_3, var_opt_set="Select Carbonate Mineral", opt_list=opt_list_0_3,
                command=lambda var_opt=var_opt_0_3: self.select_opt(var_opt))
        elif var_opt == "Halogenes":
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
            except:
                pass
            var_opt_0_4 = tk.StringVar()
            opt_list_0_4 = ["Halite", "Fluorite", "Sylvite"]
            self.opt_halogene = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                   bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_4, var_opt_set="Select Halogene Mineral", opt_list=opt_list_0_4,
                command=lambda var_opt=var_opt_0_4: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_5 = tk.StringVar()
            opt_list_0_5 = ["Alkalifeldspar", "Plagioclase", "Scapolite"]
            self.opt_afs = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_5, var_opt_set="Select Tectosilicate Mineral", opt_list=opt_list_0_5,
                command=lambda var_opt=var_opt_0_5: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_6 = tk.StringVar()
            opt_list_0_6 = ["Illite", "Kaolinite", "Montmorillonite", "Chamosite", "Clinochlore", "Pennantite",
                            "Nimite", "Chlorite", "Vermiculite", "Annite", "Phlogopite", "Eastonite", "Siderophyllite",
                            "Biotite", "Muscovite", "Glauconite"]
            opt_list_0_6.sort()
            self.opt_clays = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_6, var_opt_set="Select Phyllosilicate Mineral", opt_list=opt_list_0_6,
                command=lambda var_opt=var_opt_0_6: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_7 = tk.StringVar()
            opt_list_0_7 = ["Barite", "Celestite", "Anglesite", "Anhydrite", "Hanksite", "Gypsum", "Alunite", "Jarosite",
                            "Chalcanthite", "Kieserite", "Scheelite", "Hexahydrite"]
            opt_list_0_7.sort()
            self.opt_sulfate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_7, var_opt_set="Select Sulfate Mineral", opt_list=opt_list_0_7,
                command=lambda var_opt=var_opt_0_7: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_8 = tk.StringVar()
            opt_list_0_8 = ["Zircon", "Thorite", "Andalusite", "Kyanite", "Sillimanite", "Topaz", "Staurolite",
                            "Fayalite", "Forsterite", "Tephroite", "Calcio-Olivine", "Liebenbergite", "Olivine",
                            "Pyrope", "Almandine", "Grossular", "Andradite", "Uvarovite", "Aluminium Garnet",
                            "Calcium Garnet"]
            opt_list_0_8.sort()
            self.opt_nesosilicate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_8, var_opt_set="Select Nesosilicate Mineral", opt_list=opt_list_0_8,
                command=lambda var_opt=var_opt_0_8: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_9 = tk.StringVar()
            opt_list_0_9 = ["Epidote", "Zoisite", "Gehlenite"]
            opt_list_0_9.sort()
            self.opt_sorosilicate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_9, var_opt_set="Select Sorosilicate Mineral", opt_list=opt_list_0_9,
                command=lambda var_opt=var_opt_0_9: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_10 = tk.StringVar()
            opt_list_0_10 = ["Enstatite", "Ferrosilite", "Diopside", "Jadeite", "Aegirine", "Spodumene", "Wollastonite",
                             "Tremolite", "Actinolite", "Glaucophane", "Augite", "Riebeckite", "Arfvedsonite",
                             "Calcium Amphibole", "Sodium Amphibole", "Mg-Fe Pyroxene", "Calcium Pyroxene"]
            opt_list_0_10.sort()
            self.opt_inosilicate = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_10, var_opt_set="Select Inosilicate Mineral", opt_list=opt_list_0_10,
                command=lambda var_opt=var_opt_0_10: self.select_opt(var_opt))
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
            except:
                pass
            var_opt_0_11 = tk.StringVar()
            opt_list_0_11 = ["Apatite", "Apatite-F", "Apatite-Cl", "Apatite-OH"]
            opt_list_0_11.sort()
            self.opt_phosphates = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                       bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_11, var_opt_set="Select Inosilicate Mineral", opt_list=opt_list_0_11,
                command=lambda var_opt=var_opt_0_11: self.select_opt(var_opt))
        elif var_opt == "Siliciclastic Rocks":
            var_opt_1_1 = tk.StringVar()
            opt_list_1_1 = ["Sandstone", "Shale"]
            self.opt_silic = SE(parent=self.parent, row_id=16, column_id=0, n_rows=2, n_columns=2,
                                bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_1, var_opt_set="Select Rock", opt_list=opt_list_1_1,
                command=lambda var_opt=var_opt_1_1: self.select_opt(var_opt))
        elif var_opt == "Carbonate Rocks":
            var_opt_1_2 = tk.StringVar()
            opt_list_1_2 = ["Limestone", "Dolomite Rock"]
            self.opt_carb = SE(parent=self.parent, row_id=16, column_id=0, n_rows=2, n_columns=2,
                               bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_2, var_opt_set="Select Rock", opt_list=opt_list_1_2,
                command=lambda var_opt=var_opt_1_2: self.select_opt(var_opt))
        elif var_opt == "Igneous Rocks":
            var_opt_1_3 = tk.StringVar()
            opt_list_1_3 = ["Felsic Rock", "Intermediate Rock", "Granite", "Gabbro", "Diorite", "Granodiorite",
                            "Monzonite", "Syenite", "Tonalite", "Quartzolite", "Qz-rich Granitoid"]
            self.opt_ign = SE(parent=self.parent, row_id=16, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_3, var_opt_set="Select Rock", opt_list=opt_list_1_3,
                command=lambda var_opt=var_opt_1_3: self.select_opt(var_opt))
        elif var_opt == "Evaporite Rocks":
            var_opt_1_4 = tk.StringVar()
            opt_list_1_4 = ["Rock Salt", "Anhydrite"]
            self.opt_ign = SE(parent=self.parent, row_id=16, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_4, var_opt_set="Select Rock", opt_list=opt_list_1_4,
                command=lambda var_opt=var_opt_1_4: self.select_opt(var_opt))
        elif var_opt == "Ore Rocks":
            var_opt_1_5 = tk.StringVar()
            opt_list_1_5 = ["Kupferschiefer"]
            self.opt_ign = SE(parent=self.parent, row_id=16, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_1_5, var_opt_set="Select Rock", opt_list=opt_list_1_5,
                command=lambda var_opt=var_opt_1_5: self.select_opt(var_opt))
    #
    def change_radiobutton_mode(self, var_rb_mode):
        if var_rb_mode.get() == 0:
            print(var_rb_mode.get())
        elif var_rb_mode.get() == 1:
            print(var_rb_mode.get())
        elif var_rb_mode.get() == 2:
            print(var_rb_mode.get())
#
class Minerals:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, mineral, lbl_w, entr_w, gui_elements):
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
        self.parent_mineral = parent
        self.color_bg = color_bg
        self.color_fg = color_fg
        self.color_acc_01 = color_acc[0]
        self.color_acc_02 = color_acc[1]
        self.mineral = mineral
        self.var_rb = tk.IntVar()
        var_rb_start = 0
        self.var_rb_trace = tk.IntVar()
        var_rb_trace_0 = 2
        self.var_actual_trace = var_rb_trace_0
        self.var_entr = tk.IntVar()
        var_entr_start = 100
        self.lbl_w = lbl_w
        self.entr_w = entr_w
        self.gui_elements = gui_elements
        #
        ## Labels
        lbl_01 = SE(parent=self.parent_mineral, row_id=4, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Molar mass\n (g/mol)", relief=tk.RAISED)
        lbl_02 = SE(parent=self.parent_mineral, row_id=6, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Density\n (kg/cbm)", relief=tk.RAISED)
        lbl_03 = SE(parent=self.parent_mineral, row_id=8, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="P-wave velocity\n (km/s)", relief=tk.RAISED)
        lbl_04 = SE(parent=self.parent_mineral, row_id=10, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="S-wave velocity\n (km/s)", relief=tk.RAISED)
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
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=0, text="Histogram",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_mineral, row_id=31, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=1, text="Scatter plot",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        #
        rb_03 = SE(parent=self.parent_mineral, row_id=32, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_trace, var_rb_set=var_rb_trace_0, value_rb=2,
                                          text="Without Trace Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_trace: self.change_radiobutton(var_rb))
        rb_04 = SE(parent=self.parent_mineral, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_trace, var_rb_set=var_rb_trace_0, value_rb=3,
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
                data = Oxides(impurity="pure").create_magnetite(dict=True)
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
                data = Oxides(impurity="pure").create_hematite(dict=True)
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
            elif self.mineral == "Trevorite":
                data = Oxides(impurity="pure", data_type=True).create_trevorite()
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
            elif self.mineral == "Realgar":
                data = Sulfides(impurity="pure", data_type=True).create_realgar()
            elif self.mineral == "Orpiment":
                data = Sulfides(impurity="pure", data_type=True).create_orpiment()
            elif self.mineral == "Marcasite":
                data = Sulfides(impurity="pure", data_type=True).create_marcasite()
            elif self.mineral == "Fahlore":
                data = Sulfides(impurity="pure", data_type=True).create_fahlore()
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
            elif self.mineral == "Halite":
                data = Halogenes(impurity="pure", dict=True).create_halite()
            elif self.mineral == "Fluorite":
                data = Halogenes(impurity="pure", dict=True).create_fluorite()
            elif self.mineral == "Sylvite":
                data = Halogenes(impurity="pure", dict=True).create_sylvite()
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
            # Nesosilicates
            elif self.mineral == "Zircon":
                data = Nesosilicates(data_type=True).create_zircon()
            elif self.mineral == "Thorite":
                data = Nesosilicates(data_type=True).create_thorite()
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
            if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
                self.w_element = self.element_list["Fe"]
            elif self.mineral in ["Quartz"]:
                self.w_element = self.element_list["Si"]
            elif self.mineral in ["Galena"]:
                self.w_element = self.element_list["Pb"]
            elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
                self.w_element = self.element_list["Cu"]
            elif self.mineral in ["Cassiterite"]:
                self.w_element = self.element_list["Sn"]
            elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
                self.w_element = self.element_list["Ca"]
            elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
                self.w_element = self.element_list["Mg"]
            elif self.mineral in ["Halite"]:
                self.w_element = self.element_list["Na"]
            elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
                self.w_element = self.element_list["Cr"]
            elif self.mineral in ["Pyrolusite"]:
                self.w_element = self.element_list["Mn"]
            elif self.mineral in ["Ilmenite", "Rutile"]:
                self.w_element = self.element_list["Ti"]
            elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
                self.w_element = self.element_list["K"]
            elif self.mineral in ["Illite", "Corundum"]:
                self.w_element = self.element_list["Al"]
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
        self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
                          [self.molar_mass, self.gamma_ray, self.photoelectricity]]
        self.create_3x3_histo(parent=self.parent_mineral, data=self.data_plot, row_id=2, column_id=9, n_rows=45,
                              n_columns=9, color=self.color_mineral, labels=self.labels)
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
    def create_plot(self, parent, data, row_id, column_id, n_rows, n_columns, xlabel, color):
        #
        self.canvas_histo = None
        self.fig_histo = Figure(facecolor="#E9ECED")
        self.ax_histo = self.fig_histo.add_subplot()
        #
        self.ax_histo.axvline(x=np.mean(data), color="#E76F51", linewidth=3, linestyle="dashed")
        self.ax_histo.hist(data, bins=15, color=color, edgecolor="black")
        self.ax_histo.grid(True)
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
        self.fig = Figure(facecolor=self.color_bg)
        self.ax = self.fig.add_subplot()
        data_x = np.array(data_x)
        #
        if self.var_opt_chem.get() in ["No Selection", "Select Element"]:
            self.ax.scatter(data_x, data_y, color=color, edgecolor="black", alpha=0.5)
        else:
            if self.var_opt_chem.get() in self.list_elements:
                plot = self.ax.scatter(data_x, data_y, c=self.data_c, cmap="viridis",
                                       edgecolor="black", alpha=1)
            cbar = self.fig.colorbar(plot, format="%.0f")
            cbar.set_label(self.var_opt_chem.get()+" (%)", rotation=90)
        self.ax.grid(True)
        self.ax.set_xlim(float(0.99*min(data_x)), float(1.01*max(data_x)))
        self.ax.set_xticks(np.around(np.linspace(float(0.99*min(data_x)), float(1.01*max(data_x)), 5), 2))
        self.ax.set_axisbelow(True)
        self.ax.set_xlabel(xlabel, fontsize="small")
        self.ax.set_ylabel(ylabel, labelpad=0.5, fontsize="small")
        self.fig.subplots_adjust(bottom=0.15, left=0.26)
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
        data_all = []
        for i in range(var_entr.get()):
            # Oxides
            if self.mineral == "Quartz":
                data = Oxides(impurity="pure", data_type=True).create_quartz()
            elif self.mineral == "Magnetite":
                data = Oxides(impurity="pure").create_magnetite(dict=True)
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
            elif self.mineral == "Halite":
                data = Halogenes(impurity="pure", dict=True).create_halite()
            elif self.mineral == "Fluorite":
                data = Halogenes(impurity="pure", dict=True).create_fluorite()
            elif self.mineral == "Sylvite":
                data = Halogenes(impurity="pure", dict=True).create_sylvite()
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
            # Nesosilicates
            elif self.mineral == "Zircon":
                data = Nesosilicates(data_type=True).create_zircon()
            elif self.mineral == "Thorite":
                data = Nesosilicates(data_type=True).create_thorite()
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
            if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
                self.w_element = self.element_list["Fe"]
            elif self.mineral in ["Quartz"]:
                self.w_element = self.element_list["Si"]
            elif self.mineral in ["Galena"]:
                self.w_element = self.element_list["Pb"]
            elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
                self.w_element = self.element_list["Cu"]
            elif self.mineral in ["Cassiterite"]:
                self.w_element = self.element_list["Sn"]
            elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
                self.w_element = self.element_list["Ca"]
            elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
                self.w_element = self.element_list["Mg"]
            elif self.mineral in ["Halite"]:
                self.w_element = self.element_list["Na"]
            elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
                self.w_element = self.element_list["Cr"]
            elif self.mineral in ["Ilmenite", "Rutile"]:
                self.w_element = self.element_list["Ti"]
            elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
                self.w_element = self.element_list["K"]
            elif self.mineral in ["Illite", "Corundum"]:
                self.w_element = self.element_list["Al"]
            elif self.mineral in ["Pyrolusite"]:
                self.w_element = self.element_list["Mn"]
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
            self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
                     [self.molar_mass, self.gamma_ray, self.photoelectricity]]
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
            self.create_3x3_scatter(parent=self.parent_mineral, data_x=data_x, data=self.data_plot, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_mineral, labels=self.labels,
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
                        data = Oxides(impurity="pure").create_magnetite(dict=True)
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
                    elif self.mineral == "Halite":
                        data = Halogenes(impurity="pure", dict=True).create_halite()
                    elif self.mineral == "Fluorite":
                        data = Halogenes(impurity="pure", dict=True).create_fluorite()
                    elif self.mineral == "Sylvite":
                        data = Halogenes(impurity="pure", dict=True).create_sylvite()
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
                    # Nesosilicates
                    elif self.mineral == "Zircon":
                        data = Nesosilicates(data_type=True).create_zircon()
                    elif self.mineral == "Thorite":
                        data = Nesosilicates(data_type=True).create_thorite()
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
                    if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
                        self.w_element = self.element_list["Fe"]
                    elif self.mineral in ["Quartz"]:
                        self.w_element = self.element_list["Si"]
                    elif self.mineral in ["Galena"]:
                        self.w_element = self.element_list["Pb"]
                    elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
                        self.w_element = self.element_list["Cu"]
                    elif self.mineral in ["Cassiterite"]:
                        self.w_element = self.element_list["Sn"]
                    elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
                        self.w_element = self.element_list["Ca"]
                    elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
                        self.w_element = self.element_list["Mg"]
                    elif self.mineral in ["Halite"]:
                        self.w_element = self.element_list["Na"]
                    elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
                        self.w_element = self.element_list["Cr"]
                    elif self.mineral in ["Ilmenite", "Rutile"]:
                        self.w_element = self.element_list["Ti"]
                    elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
                        self.w_element = self.element_list["K"]
                    elif self.mineral in ["Illite", "Corundum"]:
                        self.w_element = self.element_list["Al"]
                    elif self.mineral in ["Pyrolusite"]:
                        self.w_element = self.element_list["Mn"]
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
                traces = ["Ti", "Ge", "Sn", "Al", "Fe", "Ga", "As", "B", "H", "Li", "Na", "Ag", "K", "Mg", "Cu", "Be", "Mn"]
                n_trace = rd.randint(1, len(traces))
                selection_trace = rd.sample(traces, n_trace)
            data_all = []
            for i in range(self.var_entr.get()):
                # Oxides
                if self.mineral == "Quartz":
                    data = Oxides(impurity=selection_trace, data_type=True).create_quartz()
                elif self.mineral == "Magnetite":
                    data = Oxides(impurity="pure").create_magnetite(dict=True)
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
                elif self.mineral == "Halite":
                    data = Halogenes(impurity="pure", dict=True).create_halite()
                elif self.mineral == "Fluorite":
                    data = Halogenes(impurity="pure", dict=True).create_fluorite()
                elif self.mineral == "Sylvite":
                    data = Halogenes(impurity="pure", dict=True).create_sylvite()
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
                # Nesosilicates
                elif self.mineral == "Zircon":
                    data = Nesosilicates(data_type=True).create_zircon()
                elif self.mineral == "Thorite":
                    data = Nesosilicates(data_type=True).create_thorite()
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
                if self.mineral in ["Magnetite", "Hematite", "Pyrite", "Siderite"]:
                    self.w_element = self.element_list["Fe"]
                elif self.mineral in ["Quartz"]:
                    self.w_element = self.element_list["Si"]
                elif self.mineral in ["Galena"]:
                    self.w_element = self.element_list["Pb"]
                elif self.mineral in ["Chalcopyrite", "Cuprospinel"]:
                    self.w_element = self.element_list["Cu"]
                elif self.mineral in ["Cassiterite"]:
                    self.w_element = self.element_list["Sn"]
                elif self.mineral in ["Calcite", "Dolomite", "Fluorite", "Plagioclase"]:
                    self.w_element = self.element_list["Ca"]
                elif self.mineral in ["Magnesite", "Aluminium Spinels"]:
                    self.w_element = self.element_list["Mg"]
                elif self.mineral in ["Halite"]:
                    self.w_element = self.element_list["Na"]
                elif self.mineral in ["Chromite", "Magnesiochromite", "Zincochromite", "Chromium Spinels"]:
                    self.w_element = self.element_list["Cr"]
                elif self.mineral in ["Ilmenite", "Rutile"]:
                    self.w_element = self.element_list["Ti"]
                elif self.mineral in ["Sylvite", "Alkalifeldspar"]:
                    self.w_element = self.element_list["K"]
                elif self.mineral in ["Illite", "Corundum"]:
                    self.w_element = self.element_list["Al"]
                elif self.mineral in ["Pyrolusite"]:
                    self.w_element = self.element_list["Mn"]
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
            self.data_plot = [[self.rho_b/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
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
            self.create_3x3_scatter(parent=self.parent_mineral, data_x=data_x, data=self.data_plot, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_mineral, labels=self.labels,
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
            self.create_3x3_scatter(parent=self.parent_mineral, data_x=data_x, data=self.data_plot, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_mineral, labels=self.labels,
                                    xlabel=xlabel)
            #
    #
    def __call__(self):
        return self.lbl_w, self.entr_w
#
class Rocks:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, rock, lbl_w, entr_w, gui_elements):
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
        self.parent_rock = parent
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
        self.rock = rock
        self.lbl_w = lbl_w
        self.entr_w = entr_w
        self.lbl_chem = []
        self.entr_chem = []
        self.gui_elements = gui_elements
        #
        ## Labels
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
        if self.rock == "Sandstone":
            var_phi0_start = 5
            var_phi1_start = 30
        elif self.rock in ["Shale", "Kupferschiefer"]:
            var_phi0_start = 0
            var_phi1_start = 10
        elif self.rock in ["Limestone", "Dolomite Rock"]:
            var_phi0_start = 0
            var_phi1_start = 50
        elif self.rock in ["Rock Salt", "Anhydrite", "Felsic Rock", "Intermediate Rock", "Granite", "Gabbro", "Syenite",
                           "Diorite", "Granodiorite", "Tonalite", "Monzonite", "Quartzolite", "Qz-rich Granitoid"]:
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
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=0, text="Histogram",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_rock, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=1, text="Scatter",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_03 = SE(parent=self.parent_rock, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, var_rb_set=var_rb_geochem_start, value_rb=2,
                                          text="Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        rb_04 = SE(parent=self.parent_rock, row_id=35, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, var_rb_set=var_rb_geochem_start, value_rb=3,
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
            for i in range(var_entr_start):
                if self.rock == "Sandstone":
                    data = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                elif self.rock == "Shale":
                    data = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                elif self.rock == "Limestone":
                    data = limestone(fluid="water", actualThickness=0).create_simple_limestone(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                elif self.rock == "Dolomite Rock":
                    data = dolomite(fluid="water", actualThickness=0).create_simple_dolomite(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
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
                #
                elif self.rock == "Rock Salt":
                    data = Evaporites(fluid="water", actualThickness=0).create_simple_rocksalt(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                elif self.rock == "Anhydrite":
                    data = Evaporites(fluid="water", actualThickness=0).create_simple_anhydrite(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
                elif self.rock == "Kupferschiefer":
                    data = Ores(fluid="water", actualThickness=0, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100), data_type=True).create_kupferschiefer()
                #
                data_all.append(data)
            #
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
            #
            self.list_elements = list(self.chemistry[0].keys())
            self.list_minerals = list(self.mineralogy[0].keys())
            self.elements = {}
            self.minerals = {}
            for element in self.list_elements:
                self.elements[element] = []
                for chemistry_data in self.chemistry:
                    self.elements[element].append(abs(chemistry_data[element]*100))
            for mineral in self.list_minerals:
                self.minerals[mineral] = []
                for mineralogy_data in self.mineralogy:
                    self.minerals[mineral].append(abs(mineralogy_data[mineral]*100))
            #
            self.results = [self.rho, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod, self.poisson, self.phi,
                            self.gamma_ray, self.photoelectricity]
            #
            self.var_opt_chem = tk.StringVar()
            opt_list_chem = ["No Selection"]
            opt_list_chem.extend(self.list_elements)
            opt_list_chem.extend(self.list_minerals)
            self.opt_chem = SE(parent=self.parent_rock, row_id=36, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                               fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element/Mineral",
                                                              opt_list=opt_list_chem, active_bg=self.color_acc_02,
                                                              command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
            self.gui_elements.append(self.opt_chem)
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
            self.data_plot = [[self.rho/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
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
        data_all = []
        for i in range(var_entr.get()):
            if self.rock == "Sandstone":
                data = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
            elif self.rock == "Shale":
                data = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
            elif self.rock == "Limestone":
                data = limestone(fluid="water", actualThickness=0).create_simple_limestone(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
            elif self.rock == "Dolomite Rock":
                data = dolomite(fluid="water", actualThickness=0).create_simple_dolomite(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
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
            #
            elif self.rock == "Rock Salt":
                data = Evaporites(fluid="water", actualThickness=0).create_simple_rocksalt(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
            elif self.rock == "Anhydrite":
                data = Evaporites(fluid="water", actualThickness=0).create_simple_anhydrite(dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100))
            elif self.rock == "Kupferschiefer":
                data = Ores(fluid="water", actualThickness=0, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100), data_type=True).create_kupferschiefer()
            #
            data_all.append(data)
        #
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
        #
        self.list_elements = list(self.chemistry[0].keys())
        self.list_minerals = list(self.mineralogy[0].keys())
        self.elements = {}
        self.minerals = {}
        for element in self.list_elements:
            self.elements[element] = []
            for chemistry_data in self.chemistry:
                self.elements[element].append(abs(chemistry_data[element]*100))
        for mineral in self.list_minerals:
            self.minerals[mineral] = []
            for mineralogy_data in self.mineralogy:
                self.minerals[mineral].append(abs(mineralogy_data[mineral]*100))
        #
        self.results = [self.rho, self.vP, self.vS, self.vPvS, self.bulk_mod, self.shear_mod, self.poisson, self.phi,
                        self.gamma_ray, self.photoelectricity]
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
        self.color_rock = "#7C9097"
        #
        self.data_plot = [[self.rho/1000, self.vP/1000, self.vS/1000], [self.bulk_mod, self.shear_mod, self.poisson],
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
            var_prop_0 = 8
            rb_rho = SE(parent=self.parent_rock, row_id=37, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, var_rb_set=var_prop_0, value_rb=8,
                                                       text="Density", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            rb_phi = SE(parent=self.parent_rock, row_id=38, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, var_rb_set=var_prop_0, value_rb=9,
                                                       text="Porosity", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            self.rb_prop = [rb_rho, rb_phi]
            self.gui_elements.extend([rb_rho, rb_phi])
            #
            data_x_rho = np.array(self.rho)/1000
            xlabel="Densitiy $\\varrho$ (g/ccm)"
            #
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x_rho, data=self.data_plot_scatter_rho, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels,
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
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels,
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
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels,
                                    xlabel=xlabel)
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
        var_prop_0 = 8
        rb_rho = SE(parent=self.parent_rock, row_id=37, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                    fg="black").create_radiobutton(var_rb=self.var_prop, var_rb_set=var_prop_0, value_rb=8,
                                                   text="Density", color_bg=self.color_acc_01,
                                                   command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
        rb_phi = SE(parent=self.parent_rock, row_id=38, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                    fg="black").create_radiobutton(var_rb=self.var_prop, var_rb_set=var_prop_0, value_rb=9,
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
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels,
                                    xlabel=xlabel)
            #
        else:
            self.var_rb.set(1)
            #
            data_c = []
            if var_opt in list(self.elements.keys()):
                for item in self.chemistry:
                    data_c.append(item[var_opt])
            if var_opt in list(self.minerals.keys()):
                for item in self.mineralogy:
                    data_c.append(item[var_opt])
            self.data_c = np.array(data_c)*100
            data_x = np.array(self.rho)/1000
            xlabel = "Density $\\varrho$ g/ccm"
            self.create_3x3_scatter(parent=self.parent_rock, data_x=data_x, data=self.data_plot_scatter_rho, row_id=2,
                                    column_id=9, n_rows=45, n_columns=9, color=self.color_rock, labels=self.labels,
                                    xlabel=xlabel)
    #
    def __call__(self):
        return self.lbl_w, self.entr_w
    #
    def create_custom_rock(self):
        ## Variables
        self.custom_mineralogy = {}
        self.custom_porosities = {}
        #
        ## Buttons
        btn_defmin = SE(parent=self.parent_rock, row_id=36, column_id=0, n_rows=2, n_columns=2, bg=self.color_acc_01,
                              fg="black").create_button(text="Define Mineralogy", command=lambda var_btn="Define Mineralogy": self.press_button(var_btn))
        btn_update = SE(parent=self.parent_rock, row_id=38, column_id=0, n_rows=2, n_columns=2, bg=self.color_acc_01,
                              fg="black").create_button(text="Generate Data")
        self.gui_elements.extend([btn_defmin, btn_update])
    #
    def press_button(self, var_btn):
        print(var_btn)
        if var_btn == "Define Mineralogy":
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
            self.window_custom_mineralogy.geometry("1200x960")
            self.window_custom_mineralogy.resizable(False, False)
            self.window_custom_mineralogy["bg"] = self.color_menu
            #
            self.var_custom_mineralogy = {}
            self.var_custom_mineralogy["checkbox"] = {}
            #
            ## LABELS
            # Oxides
            lbl_oxides = SE(parent=self.window_custom_mineralogy, row_id=0, column_id=0, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Oxides", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=0, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=1, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=2, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=3, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_oxides = ["Quartz", "Magnetite", "Hematite", "Corundum", "Ilmenite", "Rutile", "Pyrolusite", "Cassiterite", "Spinel", "Chromite"]
            for index, oxide in enumerate(list_oxides, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=0, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=oxide, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][oxide] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=1, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][oxide])
            # Tectosilicates
            lbl_tectosilicates = SE(parent=self.window_custom_mineralogy, row_id=0, column_id=4, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Tectosilicates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=4, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=5, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=6, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=7, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_tectosilicates = ["Alkali Feldspar", "Plagioclase", "Scapolite"]
            for index, tectosilicate in enumerate(list_tectosilicates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=4, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=tectosilicate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][tectosilicate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=5, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][tectosilicate])
            # Phyllosilicates
            lbl_phyllosilicates = SE(parent=self.window_custom_mineralogy, row_id=0, column_id=8, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Phyllosilicates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=8, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=9, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=10, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=11, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_phyllosilicates = ["Illite", "Kaolinite", "Montmorillonite", "Chlorite", "Vermiculite", "Biotite", "Muscovite", "Glauconite"]
            for index, phyllosilicate in enumerate(list_phyllosilicates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=8, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=phyllosilicate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][phyllosilicate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=9, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][phyllosilicate])
            # Nesosilicates
            lbl_nesosilicates = SE(parent=self.window_custom_mineralogy, row_id=0, column_id=12, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Nesosilicates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=12, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=13, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=14, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=15, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_nesosilicates = ["Olivine", "Garnet", "Zircon", "Thorite", "Andalusite", "Kyanite", "Sillimanite", "Topaz", "Staurolite"]
            for index, nesosilicate in enumerate(list_nesosilicates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=12, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=nesosilicate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][nesosilicate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=13, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][nesosilicate])
            # Inosilicates
            lbl_inosilicates = SE(parent=self.window_custom_mineralogy, row_id=0, column_id=16, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Inosilicates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=16, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=17, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=18, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=19, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_inosilicates = ["Enstatite", "Ferrosilite", "Diopside", "Jadeite", "Aegirine", "Spodumene",
                                 "Wollastonite", "Pyroxene", "Amphibole", "Cummingtonite", "Tremolite", "Actinolite",
                                 "Hornblende", "Glaucophane"]
            for index, inosilicate in enumerate(list_inosilicates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=16, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=inosilicate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][inosilicate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=17, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][inosilicate])
            # Sorosilicates
            lbl_sorosilicates = SE(parent=self.window_custom_mineralogy, row_id=0, column_id=20, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Sorosilicates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=20, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=21, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=22, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=1, column_id=23, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_sorosilicates = ["Epidote", "Zoisite", "Gehlenite"]
            for index, sorosilicate in enumerate(list_sorosilicates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=20, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=sorosilicate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][sorosilicate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=2+index, column_id=21, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][sorosilicate])
            # Sulfides
            lbl_sulfides = SE(parent=self.window_custom_mineralogy, row_id=16, column_id=4, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Sulfides", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=4, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=5, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=6, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=7, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_sulfides = ["Pyrite", "Chalcopyrite", "Bornite", "Covellite", "Molybdenite", "Sphalerite", "Galena", "Fahlore"]
            for index, sulfide in enumerate(list_sulfides, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=18+index, column_id=4, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=sulfide, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][sulfide] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=18+index, column_id=5, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][sulfide])
            # Carbonates
            lbl_carbonates = SE(parent=self.window_custom_mineralogy, row_id=16, column_id=0, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Carbonates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=0, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=1, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=2, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=3, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_carbonates = ["Calcite", "Dolomite", "Magnesite", "Rhodochrosite", "Siderite", "Aragonite",
                               "Cerussite", "Ankerite", "Azurite", "Malachite"]
            list_carbonates.sort()
            for index, carbonate in enumerate(list_carbonates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=18+index, column_id=0, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=carbonate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][carbonate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=18+index, column_id=1, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][carbonate])
            # Sulfates
            lbl_sulfates = SE(parent=self.window_custom_mineralogy, row_id=16, column_id=8, n_columns=4,
                            bg=self.color_accent_01, fg=self.color_fg_dark).create_label(text="Sulfates", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=8, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Name", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=9, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Part", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=10, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Min", relief=tk.RAISED)
            lbl_name = SE(parent=self.window_custom_mineralogy, row_id=17, column_id=11, bg=self.color_accent_01,
                            fg=self.color_fg_dark).create_label(text="Max", relief=tk.RAISED)
            list_sulfates = ["Barite", "Anhydrite", "Gypsum", "Alunite", "Jarosite", "Anglesite", "Hanksite",
                             "Celestite", "Kieserite", "Chalcanthite", "Scheelite"]
            list_sulfates.sort()
            for index, sulfate in enumerate(list_sulfates, start=0):
                lbl = SE(parent=self.window_custom_mineralogy, row_id=18+index, column_id=8, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_label(text=sulfate, relief=tk.RAISED)
                self.var_custom_mineralogy["checkbox"][sulfate] = tk.IntVar()
                cb = SE(parent=self.window_custom_mineralogy, row_id=18+index, column_id=9, bg=self.color_accent_02,
                         fg=self.color_fg_dark).create_checkbox(text="", var_cb=self.var_custom_mineralogy["checkbox"][sulfate],
                                                                command=lambda var_cb=self.var_custom_mineralogy["checkbox"][sulfate]: self.marked_checkbox(var_cb))
    #
    def marked_checkbox(self, var_cb):
        print(var_cb.get())
#
class Subsurface:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, subsurface, lbl_w, entr_w, gui_elements):
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
            # data = SedimentaryBasin().create_sedimentary_basin()
            # for item in data:
            #     print(item)
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
    def create_random_sequences(self, thickness, style, n_parts=20):
        #
        self.var_rb_stat = tk.IntVar()
        var_rb_stat_start = 0
        self.var_rb_geochem = tk.IntVar()
        var_rb_geochem_start = 3
        self.var_rb_lith = tk.IntVar()
        var_rb_lith_start = 5
        self.current_rb_lith = tk.IntVar()
        #
        results_subsurface = {}
        self.results_sorted = {}
        properties = ["thickness", "rock", "rho", "vP", "vS", "vPvS", "phi", "K", "G", "Poisson", "GR", "PE", "Top", "Bottom"]
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
        self.create_well_log_plot(parent=self.parent_subsurface, data_x=self.results_sorted["GR"], data_y=self.results_sorted["Top"], row_id=2, column_id=9, n_rows=45, n_columns=9)
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
            self.create_well_log_plot(parent=self.parent_subsurface, data_x=self.results_sorted["GR"], data_y=self.results_sorted["Top"], row_id=2, column_id=9, n_rows=45, n_columns=9)
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
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_rho, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Density $\\varrho$ (g/ccm)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            colors = ["tan", "olivedrab", "darkorange"]
            litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
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
            var_prop_0 = 8
            rb_rho = SE(parent=self.parent_subsurface, row_id=40, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, var_rb_set=var_prop_0, value_rb=8,
                                                       text="Density", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            rb_phi = SE(parent=self.parent_subsurface, row_id=41, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                        fg="black").create_radiobutton(var_rb=self.var_prop, var_rb_set=var_prop_0, value_rb=9,
                                                       text="Porosity", color_bg=self.color_acc_01,
                                                       command=lambda var_rb=self.var_prop: self.change_radiobutton(var_rb))
            self.rb_prop = [rb_rho, rb_phi]
            #
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
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_phi, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Porosity $\\phi$ (%)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            colors = ["tan", "olivedrab", "darkorange"]
            litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
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
            if "Sandstone" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Sandstone":
                self.list_elements_0 = self.list_elements_sst
                self.elements_0 = self.elements_sst
            elif "Shale" == self.list_rocks_short[0] or self.list_rocks_short[self.var_rb_lith.get() - 5] == "Shale":
                self.list_elements_0 = self.list_elements_sh
                self.elements_0 = self.elements_sh
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
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_phi, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Porosity $\\phi$ (%)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            colors = ["tan", "olivedrab", "darkorange"]
            litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
            #
            self.create_3x3_scatter(parent=self.parent_subsurface, data_x=data_x_rho, data=data, row_id=2, column_id=9,
                                    n_rows=45, n_columns=9, colors=colors, labels=labels,
                                    xlabel="Densitiy $\\varrho$ (g/ccm)", lithos=litho_list)
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
            data_x = [np.array(self.results_plot["Sandstone"]["phi"])*100,
                      np.array(self.results_plot["Shale"]["phi"])*100,
                      np.array(self.results_plot[self.list_rocks_short[-1]]["phi"])*100]
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
            data_rho = [np.array(self.results_plot["Sandstone"]["rho"])/1000,
                        np.array(self.results_plot["Shale"]["rho"])/1000,
                        np.array(self.results_plot[self.list_rocks_short[-1]]["rho"])/1000]
            data_gr = [self.results_plot["Sandstone"]["GR"], self.results_plot["Shale"]["GR"],
                       self.results_plot[self.list_rocks_short[-1]]["GR"]]
            data_pe = [self.results_plot["Sandstone"]["PE"], self.results_plot["Shale"]["PE"],
                       self.results_plot[self.list_rocks_short[-1]]["PE"]]
            data = [[data_vP, data_vS, data_vPvS], [data_K, data_G, data_nu], [data_rho, data_gr, data_pe]]
            #
            labels = [["Seismic velocity $v_P$ (m/s)", "Seismic velocity $v_S$ (m/s)",
                      "Seismic velocity ratio $v_P/v_S$ (1)"], ["Bulk modulus $K$ (GPa)", "Shear modulus $G$ (GPa)",
                      "Poisson's ratio $\\nu$ (1)"], ["Density $\\varrho$ (g/ccm)", "Gamma Ray $GR$ (API)",
                      "Photoelectric Effect $PE$ (barns\electron)"]]
            colors = ["tan", "olivedrab", "darkorange"]
            litho_list = ["Sandstone", "Shale", self.list_rocks_short[-1]]
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
                for k in range(3):
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
    def create_well_log_plot(self, parent, data_x, data_y, row_id, column_id, n_rows, n_columns):
        #
        self.canvas = None
        max_thickness = max(data_y)
        self.fig, (self.ax1, self.ax2, self.ax3, self.ax4, self.ax5) = plt.subplots(1, 5, sharey="row", gridspec_kw={"wspace": 0.15}, figsize=(12, 16), facecolor="#E9ECED")
        self.fig.subplots_adjust(wspace=0.25)
        # 1
        self.ax1.plot(self.results_sorted["GR"], self.results_sorted["Top"], color="#00549F", linewidth=2)
        self.ax1.set_xlabel("GR [API]")
        self.ax1.set_ylabel("Depth [m]")
        self.ax1.set_xlim(0, 200)
        self.ax1.set_xticks(np.arange(0, 250, 50))
        self.ax1.set_ylim(0, max_thickness)
        self.ax1.set_yticks(np.arange(0, max_thickness+50, 50))
        self.ax1.grid(color="grey", linestyle="dashed")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 2
        vP_edit = np.array(self.results_sorted["vP"])/1000
        vS_edit = np.array(self.results_sorted["vS"])/1000
        self.ax2.plot(vP_edit, self.results_sorted["Top"], color="#00549F", linewidth=2)
        self.ax2.set_xlabel("$v_P$ [km/s]")
        self.ax2.set_xlim(0, 8.5)
        self.ax2.set_xticks(np.arange(0, 8.5, 2.0))
        self.ax2.xaxis.label.set_color("#00549F")
        self.ax2.set_ylim(0, max_thickness)
        self.ax2.set_yticks(np.arange(0, max_thickness+50, 50))
        self.ax2.grid(color="grey", linestyle="dashed")
        self.ax2_2 = self.ax2.twiny()
        self.ax2_2.plot(vS_edit, self.results_sorted["Top"], color="#CC071E", linewidth=2)
        self.ax2_2.set_xlabel("$v_S$ [km/s]")
        self.ax2_2.set_xlim(0, 8.5)
        self.ax2_2.set_xticks(np.arange(0, 8.5, 2.0))
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
        self.ax3.set_yticks(np.arange(0, max_thickness+50, 50))
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
        self.ax4.set_ylim(0, max_thickness)
        self.ax4.set_yticks(np.arange(0, max_thickness+50, 50))
        self.ax4.grid(color="grey", linestyle="dashed", which="both")
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        # 5
        n_units = []
        units_sorted = []
        for rock in self.list_rocks_short:
            units_sorted.append([rock])
            n_units.append(sum(self.results_plot[rock]["thickness"]))
            for index, value in enumerate(self.results_plot[rock]["thickness"], start=0):
                units_sorted[-1].append([self.results_plot[rock]["Top"][index], self.results_plot[rock]["Bottom"][index]])
            if rock == "Sandstone":
                units_sorted[-1].append("tan")
            elif rock == "Shale":
                units_sorted[-1].append("olivedrab")
            elif rock in ["Granite", "Gabbro", "Diorite"]:
                units_sorted[-1].append("darkorange")
        legend_lithology = []
        for i in range(len(units_sorted)):
            legend_lithology.append(mpatches.Patch(facecolor=units_sorted[i][-1], hatch="", label=units_sorted[i][0]))
        for i in range(len(n_units)):
            for j in range(1, len(units_sorted[i])-1):
                self.ax5.hist(x=np.linspace(units_sorted[i][j][0], units_sorted[i][j][1]), bins=len(n_units),
                              color=units_sorted[i][-1], orientation="horizontal")
        self.ax5.set_xlabel("Lithology")
        self.ax5.set_xlim(0, 5)
        self.ax5.set_xticks([])
        self.ax5.set_ylim(0, max_thickness)
        self.ax5.set_yticks(np.arange(0, max_thickness+50, 50))
        self.ax5.margins(0.3, 0.0)
        plt.gca().invert_yaxis()
        plt.rc("axes", axisbelow=True)
        self.ax5.legend(handles=legend_lithology, loc="lower left", bbox_to_anchor=(0, -0.125), shadow=True, ncol=1, prop={'size': 8}, frameon=False)
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