#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gui.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.11.2021

#-----------------------------------------------

## MODULES
import tkinter as tk
from modules.gui_elements import SimpleElements as SE
from modules.oxides import Oxides
from modules.sulfides import Sulfides
from modules.carbonates import Carbonates
from modules.halogenes import Halogenes
from modules.silicates import Tectosilicates
from modules.pyllosilicates import Pyllosilicates
from modules.minerals import feldspars
from modules.siliciclastics import sandstone, shale
from modules.carbonates import limestone, dolomite
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
        self.lbl_w = []
        self.entr_w = []
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
        SE(parent=self.parent, row_id=18, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_01,
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
        opt_list_0_0 = ["Oxides", "Sulfides", "Carbonates", "Halogenes", "Tectosilicates", "Phyllosilicates"]
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
        self.opt_realseq = SE(parent=self.parent, row_id=20, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
            var_opt=var_opt_2_0, var_opt_set="Select Real Sequences", opt_list=opt_list_2_0,
            command=lambda var_opt=var_opt_2_0: self.select_opt(var_opt))
        #
        ## Button
        self.btn_randseq = SE(parent=self.parent, row_id=22, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                              fg=self.color_fg_dark).create_button(text="Create Random Sequence",
                                                                   command=lambda var_btn="random": self.pressed_button(var_btn))
        self.btn_custseq = SE(parent=self.parent, row_id=24, column_id=0, n_rows=2, n_columns=2, bg=self.color_accent_02,
                              fg=self.color_fg_dark).create_button(text="Create Custom Sequence")
    #
    def pressed_button(self, var_btn):
        if var_btn == "random":
            Subsurface(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                       color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_btn, lbl_w=self.lbl_w,
                       entr_w=self.entr_w)
        else:
            Subsurface(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                       color_acc=[self.color_accent_03, self.color_accent_04], subsurface=var_btn, lbl_w=self.lbl_w,
                       entr_w=self.entr_w)
    #
    def select_opt(self, var_opt):
        # Minerals
        if var_opt == "Quartz":
            self.lbl_w, self.entr_w = Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w)()
        elif var_opt in ["Magnetite", "Hematite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w)
        elif var_opt in ["Pyrite", "Chalcopyrite", "Galena"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w)
        elif var_opt in ["Calcite", "Dolomite", "Magnesite", "Halite", "Fluorite", "Sylvite", "Illite"]:
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w)
        elif var_opt == "Alkalifeldspar":
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w)
        elif var_opt == "Plagioclase":
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                     color_acc=[self.color_accent_03, self.color_accent_04], mineral=var_opt, lbl_w=self.lbl_w,
                     entr_w=self.entr_w)
        # Rocks
        elif var_opt == "Sandstone":
            self.lbl_w, self.entr_w = Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)()
        elif var_opt == "Shale":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        #
        elif var_opt == "Limestone":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Dolomite Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        #
        elif var_opt == "Felsic Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Intermediate Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Granite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Gabbro":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Diorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Granodiorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Syenite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Tonalite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Monzonite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Quartzolite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Qz-rich Granitoid":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        #
        elif var_opt == "Rock Salt":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt == "Anhydrite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        elif var_opt in ["Kupferschiefer"]:
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg_light,
                  color_acc=[self.color_accent_03, self.color_accent_04], rock=var_opt, lbl_w=self.lbl_w,
                  entr_w=self.entr_w)
        #
        elif var_opt == "Oxides":
            try:
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
            except:
                pass
            var_opt_0_1 = tk.StringVar()
            opt_list_0_1 = ["Quartz", "Magnetite", "Hematite"]
            self.opt_oxide = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_1, var_opt_set="Select Oxide", opt_list=opt_list_0_1, active_bg=self.color_accent_02,
                command=lambda var_opt=var_opt_0_1: self.select_opt(var_opt))
        elif var_opt == "Sulfides":
            try:
                self.opt_oxide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
            except:
                pass
            var_opt_0_2 = tk.StringVar()
            opt_list_0_2 = ["Pyrite", "Chalcopyrite", "Galena"]
            self.opt_sulfide = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                  bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_2, var_opt_set="Select Sulfide", opt_list=opt_list_0_2,
                command=lambda var_opt=var_opt_0_2: self.select_opt(var_opt))
        elif var_opt == "Carbonates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
            except:
                pass
            var_opt_0_3 = tk.StringVar()
            opt_list_0_3 = ["Calcite", "Dolomite", "Magnesite"]
            self.opt_carb = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                               bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_3, var_opt_set="Select Carbonate", opt_list=opt_list_0_3,
                command=lambda var_opt=var_opt_0_3: self.select_opt(var_opt))
        elif var_opt == "Halogenes":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_clays.grid_remove()
                self.opt_afs.grid_remove()
            except:
                pass
            var_opt_0_4 = tk.StringVar()
            opt_list_0_4 = ["Halite", "Fluorite", "Sylvite"]
            self.opt_halogene = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                                   bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_4, var_opt_set="Select Halogene", opt_list=opt_list_0_4,
                command=lambda var_opt=var_opt_0_4: self.select_opt(var_opt))
        elif var_opt == "Tectosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_clays.grid_remove()
            except:
                pass
            var_opt_0_5 = tk.StringVar()
            opt_list_0_5 = ["Alkalifeldspar", "Plagioclase"]
            self.opt_afs = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_5, var_opt_set="Select Tectosilicate", opt_list=opt_list_0_5,
                command=lambda var_opt=var_opt_0_5: self.select_opt(var_opt))
        elif var_opt == "Phyllosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
                self.opt_afs.grid_remove()
            except:
                pass
            var_opt_0_6 = tk.StringVar()
            opt_list_0_6 = ["Illite"]
            self.opt_clays = SE(parent=self.parent, row_id=10, column_id=0, n_rows=2, n_columns=2,
                              bg=self.color_accent_02, fg=self.color_fg_dark).create_option_menu(
                var_opt=var_opt_0_6, var_opt_set="Select Pyllosilicates", opt_list=opt_list_0_6,
                command=lambda var_opt=var_opt_0_6: self.select_opt(var_opt))
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
    def __init__(self, parent, color_bg, color_fg, color_acc, mineral, lbl_w, entr_w):
        #
        try:
            for lbl in lbl_w:
                lbl.grid_forget()
            for entr in entr_w:
                entr.grid_forget()
            lbl_w.clear()
            entr_w.clear()
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
        self.var_entr = tk.IntVar()
        var_entr_start = 100
        self.lbl_w = lbl_w
        self.entr_w = entr_w
        #
        ## Labels
        lbl_01 = SE(parent=self.parent_mineral, row_id=4, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Molar mass\n (g/mol)", relief=tk.RAISED)
        lbl_02 = SE(parent=self.parent_mineral, row_id=6, column_id=3, n_rows=2, bg=self.color_bg,
           fg="black").create_label(text="Density\n (g/ccm)", relief=tk.RAISED)
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
        lbl_11 = SE(parent=self.parent_mineral, row_id=28, column_id=0, bg=self.color_acc_01,
           fg="black").create_label(text="Number of samples", relief=tk.RAISED)
        #
        lbl_12 = SE(parent=self.parent_mineral, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                 fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        #
        self.lbl_w.extend([lbl_01, lbl_02, lbl_03, lbl_04, lbl_05, lbl_06, lbl_07, lbl_08, lbl_09, lbl_10, lbl_11, lbl_12])
        #
        ## Entry
        entr_01 = SE(parent=self.parent_mineral, row_id=28, column_id=1, bg=self.color_acc_02,
           fg=color_fg).create_entry(var_entr=self.var_entr, var_entr_set=var_entr_start,
                                     command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        self.entr_w.append(entr_01)
        #
        ## Radiobuttons
        rb_01 = SE(parent=self.parent_mineral, row_id=29, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=0, text="Histogram",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_mineral, row_id=30, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=1, text="Scatter plot",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        self.lbl_w.extend([rb_01, rb_02])
        #
        self.var_dict = False
        data_all = []
        for i in range(var_entr_start):
            if self.mineral == "Quartz":
                data = Oxides(impurity="pure").create_quartz()
            elif self.mineral == "Magnetite":
                self.var_dict = True
                data = Oxides(impurity="pure").create_magnetite(dict=self.var_dict)
            elif self.mineral == "Hematite":
                self.var_dict = True
                data = Oxides(impurity="pure").create_hematite(dict=self.var_dict)
            elif self.mineral == "Pyrite":
                self.var_dict = True
                data = Sulfides(impurity="pure", dict=self.var_dict).create_pyrite()
            elif self.mineral == "Chalcopyrite":
                self.var_dict = True
                data = Sulfides(impurity="pure", dict=self.var_dict).create_chalcopyrite()
            elif self.mineral == "Galena":
                self.var_dict = True
                data = Sulfides(impurity="pure", dict=self.var_dict).create_galena()
            elif self.mineral == "Calcite":
                self.var_dict = True
                data = Carbonates(impurity="pure", dict=self.var_dict).create_calcite()
            elif self.mineral == "Dolomite":
                self.var_dict = True
                data = Carbonates(impurity="pure", dict=self.var_dict).create_dolomite()
            elif self.mineral == "Magnesite":
                self.var_dict = True
                data = Carbonates(impurity="pure", dict=self.var_dict).create_magnesite()
            elif self.mineral == "Halite":
                self.var_dict = True
                data = Halogenes(impurity="pure", dict=self.var_dict).create_halite()
            elif self.mineral == "Fluorite":
                self.var_dict = True
                data = Halogenes(impurity="pure", dict=self.var_dict).create_fluorite()
            elif self.mineral == "Sylvite":
                self.var_dict = True
                data = Halogenes(impurity="pure", dict=self.var_dict).create_sylvite()
            elif self.mineral == "Illite":
                self.var_dict = True
                data = Pyllosilicates(impurity="pure", dict=self.var_dict).create_illite()
            elif self.mineral == "Alkalifeldspar":
                data = Tectosilicates(impurity="pure").create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data = Tectosilicates(impurity="pure").create_plagioclase()
            self.color_mineral = "#7C9097"
            #
            data_all.append(data)
        #
        if self.var_dict == False:
            elements = np.array(data_all[0][-1])[:, 0]
            self.element_list = []
            for element in elements:
                self.element_list.append(DP(dataset=data_all).extract_element_amounts(type="mineral", element=element))
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
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Si")
            elif self.mineral == "Alkalifeldspar":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="K")
            elif self.mineral == "Plagioclase":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Ca")
        else:
            self.molar_mass = DP(dataset=data_all).extract_data(keyword="M")
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
                    self.element_list[element].append(chem_dict[element]*100)
            if self.mineral in ["Magnetite", "Hematite", "Pyrite"]:
                self.w_element = self.element_list["Fe"]
            elif self.mineral in ["Galena"]:
                self.w_element = self.element_list["Pb"]
            elif self.mineral in ["Chalcopyrite"]:
                self.w_element = self.element_list["Cu"]
            elif self.mineral in ["Calcite", "Dolomite", "Fluorite"]:
                self.w_element = self.element_list["Ca"]
            elif self.mineral in ["Magnesite"]:
                self.w_element = self.element_list["Mg"]
            elif self.mineral in ["Halite"]:
                self.w_element = self.element_list["Na"]
            elif self.mineral in ["Sylvite"]:
                self.w_element = self.element_list["K"]
            elif self.mineral in ["Illite"]:
                self.w_element = self.element_list["Al"]
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
                self.entr_w.append(entr_min)
                self.entr_w.append(entr_max)
                self.entr_w.append(entr_mean)
                self.entr_w.append(entr_std)
        #
        self.create_plot(parent=self.parent_mineral, data=self.rho_b/1000, row_id=2, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vP/1000, row_id=2, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vS/1000, row_id=2, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.bulk_mod, row_id=17, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.shear_mod, row_id=17, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.poisson, row_id=17, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.molar_mass, row_id=32, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.gamma_ray, row_id=32, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.photoelectricity, row_id=32, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_mineral)
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
        #
        self.ax.scatter(data_x, data_y, color=color, edgecolor="black")
        self.ax.grid(True)
        self.ax.set_xticks(np.around(np.linspace(0, max(data_x), 5), 2))
        self.ax.set_axisbelow(True)
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
            self.fig_histo.clf()
            self.ax_histo.cla()
            self.canvas_histo.get_tk_widget().pack_forget()
        except AttributeError:
            pass

        try:
            if self.canvas_histo:
                self.canvas_histo.destroy()
        except AttributeError:
            pass
        #
        data_all = []
        for i in range(var_entr.get()):
            if self.mineral == "Quartz":
                data = Oxides(impurity="pure").create_quartz()
            elif self.mineral == "Magnetite":
                self.var_dict = True
                data = Oxides(impurity="pure").create_magnetite(dict=self.var_dict)
            elif self.mineral == "Hematite":
                self.var_dict = True
                data = Oxides(impurity="pure").create_hematite(dict=self.var_dict)
            elif self.mineral == "Pyrite":
                self.var_dict = True
                data = Sulfides(impurity="pure", dict=self.var_dict).create_pyrite()
            elif self.mineral == "Chalcopyrite":
                self.var_dict = True
                data = Sulfides(impurity="pure", dict=self.var_dict).create_chalcopyrite()
            elif self.mineral == "Galena":
                self.var_dict = True
                data = Sulfides(impurity="pure", dict=self.var_dict).create_galena()
            elif self.mineral == "Calcite":
                self.var_dict = True
                data = Carbonates(impurity="pure", dict=self.var_dict).create_calcite()
            elif self.mineral == "Dolomite":
                self.var_dict = True
                data = Carbonates(impurity="pure", dict=self.var_dict).create_dolomite()
            elif self.mineral == "Magnesite":
                self.var_dict = True
                data = Carbonates(impurity="pure", dict=self.var_dict).create_magnesite()
            elif self.mineral == "Halite":
                self.var_dict = True
                data = Halogenes(impurity="pure", dict=self.var_dict).create_halite()
            elif self.mineral == "Fluorite":
                self.var_dict = True
                data = Halogenes(impurity="pure", dict=self.var_dict).create_fluorite()
            elif self.mineral == "Sylvite":
                self.var_dict = True
                data = Halogenes(impurity="pure", dict=self.var_dict).create_sylvite()
            elif self.mineral == "Illite":
                self.var_dict = True
                data = Pyllosilicates(impurity="pure", dict=self.var_dict).create_illite()
            elif self.mineral == "Alkalifeldspar":
                data = Tectosilicates(impurity="pure").create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data = Tectosilicates(impurity="pure").create_plagioclase()
            self.color_mineral = "#7C9097"
            #
            data_all.append(data)
        #
        if self.var_dict == False:
            elements = np.array(data_all[0][-1])[:, 0]
            self.element_list = []
            for element in elements:
                self.element_list.append(DP(dataset=data_all).extract_element_amounts(type="mineral", element=element))
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
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Si")
            elif self.mineral == "Alkalifeldspar":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="K")
            elif self.mineral == "Plagioclase":
                self.w_element = DP(dataset=data_all).extract_element_amounts(type="mineral", element="Ca")
        else:
            self.molar_mass = DP(dataset=data_all).extract_data(keyword="M")
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
                    self.element_list[element].append(chem_dict[element]*100)
            if self.mineral in ["Magnetite", "Hematite", "Pyrite"]:
                self.w_element = self.element_list["Fe"]
            elif self.mineral in ["Galena"]:
                self.w_element = self.element_list["Pb"]
            elif self.mineral in ["Chalcopyrite"]:
                self.w_element = self.element_list["Cu"]
            elif self.mineral in ["Calcite", "Dolomite", "Fluorite"]:
                self.w_element = self.element_list["Ca"]
            elif self.mineral in ["Magnesite"]:
                self.w_element = self.element_list["Mg"]
            elif self.mineral in ["Halite"]:
                self.w_element = self.element_list["Na"]
            elif self.mineral in ["Sylvite"]:
                self.w_element = self.element_list["K"]
            elif self.mineral in ["Illite"]:
                self.w_element = self.element_list["Al"]
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
                self.entr_w.append(entr_min)
                self.entr_w.append(entr_max)
                self.entr_w.append(entr_mean)
                self.entr_w.append(entr_std)
        #
        self.create_plot(parent=self.parent_mineral, data=self.rho_b/1000, row_id=2, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vP/1000, row_id=2, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vS/1000, row_id=2, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.bulk_mod, row_id=17, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.shear_mod, row_id=17, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.poisson, row_id=17, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.molar_mass, row_id=32, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.gamma_ray, row_id=32, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.photoelectricity, row_id=32, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_mineral)
        #
        self.var_rb.set(0)
    #
    def change_radiobutton(self, var_rb):
        if var_rb.get() == 0:
            try:
                self.fig_histo.clf()
                self.ax_histo.cla()
                self.canvas_histo.get_tk_widget().pack_forget()
            except AttributeError:
                pass

            try:
                if self.canvas_histo:
                    self.canvas_histo.destroy()
            except AttributeError:
                pass
            #
            self.create_plot(parent=self.parent_mineral, data=self.rho_b/1000, row_id=2, column_id=9, n_rows=15,
                             n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.vP/1000, row_id=2, column_id=12, n_rows=15,
                             n_columns=3, xlabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.vS/1000, row_id=2, column_id=15, n_rows=15,
                             n_columns=3, xlabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.bulk_mod, row_id=17, column_id=9, n_rows=15,
                             n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.shear_mod, row_id=17, column_id=12, n_rows=15,
                             n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.poisson, row_id=17, column_id=15, n_rows=15,
                             n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.molar_mass, row_id=32, column_id=9, n_rows=15,
                             n_columns=3, xlabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.gamma_ray, row_id=32, column_id=12, n_rows=15,
                             n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.photoelectricity, row_id=32, column_id=15, n_rows=15,
                             n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_mineral)
        elif var_rb.get() == 1:
            try:
                self.fig_histo.clf()
                self.ax_histo.cla()
                self.canvas_histo.get_tk_widget().pack_forget()
            except AttributeError:
                pass

            try:
                if self.canvas_histo:
                    self.canvas_histo.destroy()
            except AttributeError:
                pass
            #
            if self.mineral == "Quartz":
                element = "Si"
            elif self.mineral in ["Magnetite", "Hematite", "Pyrite"]:
                element = "Fe"
            elif self.mineral in ["Chalcopyrite"]:
                element = "Cu"
            elif self.mineral in ["Galena"]:
                element = "Pb"
            elif self.mineral in ["Magnesite"]:
                element = "Mg"
            elif self.mineral in ["Halite"]:
                element = "Na"
            elif self.mineral in ["Illite"]:
                element = "Al"
            elif self.mineral in ["Alkalifeldspar", "Sylvite"]:
                element = "K"
            elif self.mineral in ["Plagioclase", "Calcite", "Dolomite", "Fluorite"]:
                element = "Ca"
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.vP/1000, row_id=2,
                                     column_id=9, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.vS/1000, row_id=2,
                                     column_id=12, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.vP/self.vS,
                                     row_id=2, column_id=15, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Velocity ratio $v_P/v_S$ (1)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.bulk_mod, row_id=17,
                                     column_id=9, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.shear_mod,
                                     row_id=17, column_id=12, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.poisson, row_id=17,
                                     column_id=15, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Poisson's ratio $\\mu$ (1)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.molar_mass,
                                     row_id=32, column_id=9, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.gamma_ray,
                                     row_id=32, column_id=12, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Gamma ray GR (API)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.photoelectricity,
                                     row_id=32, column_id=15, n_rows=15, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Photoelectricity PE (barns/electron)", color=self.color_mineral)
    #
    def __call__(self):
        return self.lbl_w, self.entr_w
#
class Rocks:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, rock, lbl_w, entr_w):
        #
        try:
            for lbl in lbl_w:
                lbl.grid_forget()
            for entr in entr_w:
                entr.grid_forget()
            lbl_w.clear()
            entr_w.clear()
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
        self.lbl_w.extend([lbl_01, lbl_02, lbl_03, lbl_04, lbl_05, lbl_06, lbl_07, lbl_08, lbl_09, lbl_10])
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
        lbl_11 = SE(parent=self.parent_rock, row_id=28, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
           fg="black").create_label(text="Number of samples", relief=tk.RAISED)
        SE(parent=self.parent_rock, row_id=28, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
           fg=color_fg).create_entry(var_entr=self.var_entr, var_entr_set=var_entr_start,
                                     command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        lbl_12 = SE(parent=self.parent_rock, row_id=29, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
           fg="black").create_label(text="Minimum Porosity (%)", relief=tk.RAISED)
        entr_01 = SE(parent=self.parent_rock, row_id=29, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
           fg=color_fg).create_entry(var_entr=self.var_phi0, var_entr_set=var_phi0_start,
                                     command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        lbl_13 = SE(parent=self.parent_rock, row_id=30, column_id=0, n_rows=1, n_columns=1, bg=self.color_acc_01,
           fg="black").create_label(text="Maximum Porosity (%)", relief=tk.RAISED)
        entr_02 = SE(parent=self.parent_rock, row_id=30, column_id=1, n_rows=1, n_columns=1, bg=self.color_acc_02,
           fg=color_fg).create_entry(var_entr=self.var_phi1, var_entr_set=var_phi1_start,
                                     command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        rb_01 = SE(parent=self.parent_rock, row_id=31, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=0, text="Histogram",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        rb_02 = SE(parent=self.parent_rock, row_id=32, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=1, text="Scatter",
                                           color_bg=self.color_acc_01,
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        #
        self.lbl_w.extend([lbl_11, lbl_12, lbl_13])
        self.lbl_w.extend([rb_01, rb_02])
        self.entr_w.extend([entr_01, entr_02])
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
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_felsic()
            elif self.rock == "Intermediate Rock":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_intermediate()
            elif self.rock == "Granite":
                data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granite()
            elif self.rock == "Gabbro":
                data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_gabbro()
            elif self.rock == "Diorite":
                data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_diorite()
            elif self.rock == "Granodiorite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granodiorite()
            elif self.rock == "Monzonite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_monzonite()
            elif self.rock == "Quartzolite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzolite()
            elif self.rock == "Qz-rich Granitoid":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzrich_granitoid()
            elif self.rock == "Syenite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_syenite()
            elif self.rock == "Tonalite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_tonalite()
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
        self.opt_chem = SE(parent=self.parent_rock, row_id=35, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
                           fg="black").create_option_menu(var_opt=self.var_opt_chem, var_opt_set="Select Element/Mineral",
                                                          opt_list=opt_list_chem, active_bg=self.color_acc_02,
                                                          command=lambda var_opt=self.var_opt_chem: self.select_opt(var_opt))
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
            self.entr_w.append(entr_min)
            self.entr_w.append(entr_max)
            self.entr_w.append(entr_mean)
            self.entr_w.append(entr_std)
            self.entr_chem.extend([entr_min, entr_max, entr_mean, entr_std])
        #
        rb_03 = SE(parent=self.parent_rock, row_id=33, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, var_rb_set=var_rb_geochem_start, value_rb=2,
                                          text="Elements", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        rb_04 = SE(parent=self.parent_rock, row_id=34, column_id=0, n_rows=1, n_columns=2, bg=self.color_acc_01,
           fg="black").create_radiobutton(var_rb=self.var_rb_geochem, var_rb_set=var_rb_geochem_start, value_rb=3,
                                          text="Minerals", color_bg=self.color_acc_01,
                                          command=lambda var_rb=self.var_rb_geochem: self.change_radiobutton(var_rb))
        #
        self.lbl_w.extend([rb_03, rb_04])
        #
        lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                 fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        self.lbl_w.append(lbl)
        self.lbl_chem.append(lbl)
        for index, element in enumerate(self.list_elements, start=0):
            if element not in ["U"]:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element), relief=tk.RAISED)
            else:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
            self.lbl_w.append(lbl)
            self.lbl_chem.append(lbl)
        #
        self.color_rock = "#7C9097"
        #
        self.create_plot(parent=self.parent_rock, data=self.rho, row_id=2, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.vP, row_id=2, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.vS, row_id=2, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.bulk_mod, row_id=17, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.shear_mod, row_id=17, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.poisson, row_id=17, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.phi*100, row_id=32, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Porosity $\\phi$ (%)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.gamma_ray, row_id=32, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.photoelectricity, row_id=32, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
    #
    def create_plot(self, parent, data, row_id, column_id, n_rows, n_columns, xlabel, color):
        #
        self.canvas_histo = None
        self.fig_histo = Figure(facecolor="#E9ECED")
        self.ax_histo = self.fig_histo.add_subplot()
        #
        self.ax_histo.axvline(x=np.mean(data), color="#E76F51", linewidth=3, linestyle="dashed")
        self.ax_histo.hist(data, bins=15, color=color, edgecolor="black")
        # if self.var_opt_chem.get() in ["No Selection", "Select Element/Mineral"]:
        #     self.ax_histo.hist(data, bins=15, color=color, edgecolor="black")
        # else:
        #     if self.var_opt_chem.get() in self.list_elements:
        #         n, bins, patches = self.ax_histo.hist(data, bins=15, color=color, edgecolor="black")
        #         for i, p in enumerate(patches):
        #             cm = plt.cm.get_cmap("viridis")
        #             plt.setp(p, 'facecolor', cm(i/25))
        #     elif self.var_opt_chem.get() in self.list_minerals:
        #         n, bins, patches = self.ax_histo.hist(data, bins=15, color=color, edgecolor="black")
        #         for i, p in enumerate(patches):
        #             cm = plt.cm.get_cmap("viridis")
        #             plt.setp(p, 'facecolor', cm(i/25))
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
        self.fig = Figure(facecolor="#E9ECED")
        self.ax = self.fig.add_subplot()
        #
        if self.var_opt_chem.get() in ["No Selection", "Select Element/Mineral"]:
            self.ax.scatter(data_x/1000, data_y, color=color, edgecolor="black", alpha=0.5)
        else:
            if self.var_opt_chem.get() in self.list_elements:
                plot = self.ax.scatter(data_x/1000, data_y, c=self.elements[self.var_opt_chem.get()], cmap="viridis",
                                       edgecolor="black", alpha=1)
            elif self.var_opt_chem.get() in self.list_minerals:
                plot = self.ax.scatter(data_x/1000, data_y, c=self.minerals[self.var_opt_chem.get()], cmap="viridis",
                                       edgecolor="black", alpha=1)
            cbar = self.fig.colorbar(plot)
            cbar.set_label(self.var_opt_chem.get()+" (%)", rotation=90)
        self.ax.grid(True)
        self.ax.set_axisbelow(True)
        self.ax.set_xlim(0.9*min(data_x/1000), 1.1*max(data_x/1000))
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
            self.fig_histo.clf()
            self.ax_histo.cla()
            self.canvas_histo.get_tk_widget().pack_forget()
        except AttributeError:
            pass
        #
        try:
            if self.canvas_histo:
                self.canvas_histo.destroy()
        except AttributeError:
            pass
        #
        try:
            for lbl in self.lbl_chem:
                lbl.grid_forget()
            for entr in self.entr_chem:
                entr.grid_forget()
            self.lbl_chem.clear()
            self.entr_chem.clear()
        except:
            pass
        #
        self.var_rb_geochem.set(2)
        lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
               fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
        self.lbl_w.append(lbl)
        self.lbl_chem.append(lbl)
        for index, element in enumerate(self.list_elements, start=0):
            if element not in ["U"]:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element), relief=tk.RAISED)
            else:
                lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                         fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
            self.lbl_w.append(lbl)
            self.lbl_chem.append(lbl)
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
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_felsic()
            elif self.rock == "Intermediate Rock":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_intermediate()
            elif self.rock == "Granite":
                data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granite()
            elif self.rock == "Gabbro":
                data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_gabbro()
            elif self.rock == "Diorite":
                data = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_diorite()
            elif self.rock == "Granodiorite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_granodiorite()
            elif self.rock == "Monzonite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_monzonite()
            elif self.rock == "Quartzolite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzolite()
            elif self.rock == "Qz-rich Granitoid":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_quartzrich_granitoid()
            elif self.rock == "Syenite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_syenite()
            elif self.rock == "Tonalite":
                data = Plutonic(fluid="water", actualThickness=0, dict=True, porosity=rd.uniform(self.var_phi0.get()/100, self.var_phi1.get()/100)).create_simple_tonalite()
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
            self.entr_w.append(entr_min)
            self.entr_w.append(entr_max)
            self.entr_w.append(entr_mean)
            self.entr_w.append(entr_std)
            self.entr_chem.extend([entr_min, entr_max, entr_mean, entr_std])
        #
        self.color_rock = "#7C9097"
        #
        self.create_plot(parent=self.parent_rock, data=self.rho, row_id=2, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.vP, row_id=2, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.vS, row_id=2, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.bulk_mod, row_id=17, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.shear_mod, row_id=17, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.poisson, row_id=17, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.phi*100, row_id=32, column_id=9, n_rows=15,
                         n_columns=3, xlabel="Porosity $\\phi$ (%)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.gamma_ray, row_id=32, column_id=12, n_rows=15,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_rock)
        self.create_plot(parent=self.parent_rock, data=self.photoelectricity, row_id=32, column_id=15, n_rows=15,
                         n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
        #
        self.var_rb.set(0)
    #
    def change_radiobutton(self, var_rb):
        if var_rb.get() == 0:
            try:
                self.fig_histo.clf()
                self.ax_histo.cla()
                self.canvas_histo.get_tk_widget().pack_forget()
            except AttributeError:
                pass

            try:
                if self.canvas_histo:
                    self.canvas_histo.destroy()
            except AttributeError:
                pass
            #
            self.create_plot(parent=self.parent_rock, data=self.rho, row_id=2, column_id=9, n_rows=15,
                             n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.vP, row_id=2, column_id=12, n_rows=15,
                             n_columns=3, xlabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.vS, row_id=2, column_id=15, n_rows=15,
                             n_columns=3, xlabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.bulk_mod, row_id=17, column_id=9, n_rows=15,
                             n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.shear_mod, row_id=17, column_id=12, n_rows=15,
                             n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.poisson, row_id=17, column_id=15, n_rows=15,
                             n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.phi*100, row_id=32, column_id=9, n_rows=15,
                             n_columns=3, xlabel="Porosity $\\phi$ (%)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.gamma_ray, row_id=32, column_id=12, n_rows=15,
                             n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_rock)
            self.create_plot(parent=self.parent_rock, data=self.photoelectricity, row_id=32, column_id=15, n_rows=15,
                             n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
        elif var_rb.get() == 1:
            try:
                self.fig_histo.clf()
                self.ax_histo.cla()
                self.canvas_histo.get_tk_widget().pack_forget()
            except AttributeError:
                pass

            try:
                if self.canvas_histo:
                    self.canvas_histo.destroy()
            except AttributeError:
                pass
            #
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.vP, row_id=2,
                                     column_id=9, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.vS, row_id=2,
                                     column_id=12, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.vP/self.vS, row_id=2,
                                     column_id=15, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Velocity ratio $v_P/v_S$ (1)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.bulk_mod, row_id=17,
                                     column_id=9, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.shear_mod, row_id=17,
                                     column_id=12, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Shear modulus $G$ (GPa)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.poisson, row_id=17,
                                     column_id=15, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.phi*100, row_id=32,
                                     column_id=9, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Porosity $\\phi$ (%)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.gamma_ray, row_id=32,
                                     column_id=12, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Gamma ray GR (API)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.photoelectricity,
                                     row_id=32, column_id=15, n_rows=15, n_columns=3,
                                     xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
        elif var_rb.get() == 2:
            #
            try:
                for lbl in self.lbl_chem:
                    lbl.grid_forget()
                for entr in self.entr_chem:
                    entr.grid_forget()
                self.lbl_chem.clear()
                self.entr_chem.clear()
            except:
                pass
            #
            lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                   fg="black").create_label(text="Chemical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w.append(lbl)
            self.lbl_chem.append(lbl)
            for index, element in enumerate(self.list_elements, start=0):
                if element not in ["U"]:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(element)+" (ppm)", relief=tk.RAISED)
                self.lbl_w.append(lbl)
                self.lbl_chem.append(lbl)
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
                self.entr_w.append(entr_min)
                self.entr_w.append(entr_max)
                self.entr_w.append(entr_mean)
                self.entr_w.append(entr_std)
                self.entr_chem.extend([entr_min, entr_max, entr_mean, entr_std])
        #
        elif var_rb.get() == 3:
            #
            try:
                for lbl in self.lbl_chem:
                    lbl.grid_forget()
                for entr in self.entr_chem:
                    entr.grid_forget()
                self.lbl_chem.clear()
                self.entr_chem.clear()
            except:
                pass
            #
            lbl = SE(parent=self.parent_rock, row_id=24, column_id=3, n_columns=5, bg=self.color_bg,
                   fg="black").create_label(text="Mineralogical composition (weight amounts %)", relief=tk.RAISED)
            self.lbl_w.append(lbl)
            self.lbl_chem.append(lbl)
            for index, mineral in enumerate(self.list_minerals, start=0):
                if mineral not in ["Urn"]:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(mineral), relief=tk.RAISED)
                else:
                    lbl = SE(parent=self.parent_rock, row_id=25+index, column_id=3, bg=self.color_bg,
                             fg="black").create_label(text=str(mineral)+" (ppm)", relief=tk.RAISED)
                self.lbl_w.append(lbl)
                self.lbl_chem.append(lbl)
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
                self.entr_w.append(entr_min)
                self.entr_w.append(entr_max)
                self.entr_w.append(entr_mean)
                self.entr_w.append(entr_std)
                self.entr_chem.extend([entr_min, entr_max, entr_mean, entr_std])
    #
    def select_opt(self, var_opt):
        try:
            self.fig_histo.clf()
            self.ax_histo.cla()
            self.canvas_histo.get_tk_widget().pack_forget()
        except AttributeError:
            pass

        try:
            if self.canvas_histo:
                self.canvas_histo.destroy()
        except AttributeError:
            pass
        #
        self.var_rb.set(1)
        #
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.vP, row_id=2,
                                 column_id=9, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.vS, row_id=2,
                                 column_id=12, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.vP/self.vS, row_id=2,
                                 column_id=15, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Velocity ratio $v_P/v_S$ (1)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.bulk_mod, row_id=17,
                                 column_id=9, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.shear_mod, row_id=17,
                                 column_id=12, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Shear modulus $G$ (GPa)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.poisson, row_id=17,
                                 column_id=15, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.phi*100, row_id=32,
                                 column_id=9, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Porosity $\\phi$ (%)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.gamma_ray, row_id=32,
                                 column_id=12, n_rows=15, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Gamma ray GR (API)", color=self.color_rock)
        self.create_scatter_plot(parent=self.parent_rock, data_x=self.rho, data_y=self.photoelectricity,
                                 row_id=32, column_id=15, n_rows=15, n_columns=3,
                                 xlabel="Densitiy $\\varrho$ (g/ccm)",
                                 ylabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
    #
    def __call__(self):
        return self.lbl_w, self.entr_w
#
class Subsurface:
    #
    def __init__(self, parent, color_bg, color_fg, color_acc, subsurface, lbl_w, entr_w):
        #
        try:
            for lbl in lbl_w:
                lbl.grid_forget()
            for entr in entr_w:
                entr.grid_forget()
            lbl_w.clear()
            entr_w.clear()
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
        self.lbl_chem = []
        self.entr_chem = []
        #
        if self.subsurface == "random":
            self.create_random_sequences(thickness=1000, style="siliciclastic")
            # data = SedimentaryBasin().create_sedimentary_basin()
            # for item in data:
            #     print(item)
        #
    def create_random_sequences(self, thickness, style):
        #
        results_subsurface = {}
        #
        if style == "siliciclastic":
            rocks = ["sandstone", "shale"]
            basement = ["granite", "gabbro", "diorite"]
            list_rocks = []
            #
            n_units = rd.randint(10, 20)
            for i in range(n_units):
                if len(list_rocks) < n_units - 1:
                    magicnumber = rd.randint(0, len(rocks)-1)
                    list_rocks.append(rocks[magicnumber])
                else:
                    magicnumber = rd.randint(0, len(basement)-1)
                    list_rocks.append(basement[magicnumber])
            list_thickness = self.split_thickness(thickness=thickness, n_units=n_units)
            #
            for index, rock in enumerate(list_rocks, start=0):
                results_subsurface[index] = {}
                if rock == "sandstone":
                    if index == 0:
                        thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=10)
                        for i, d in enumerate(thickness_unit, start=0):
                            if i == 0:
                                results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                results_subsurface[index][i]["thickness"] = d
                            else:
                                results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                   amounts=results_subsurface[index][i-1]["mineralogy"])
                                results_subsurface[index][i]["thickness"] = d
                    else:
                        thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=10)
                        for i, d in enumerate(thickness_unit, start=0):
                            if i == 0:
                                if list_rocks[index-1] == "shale":
                                    magicnumber = rd.randint(0, 2)
                                    if magicnumber == 0:
                                        results_subsurface[index][i] = sandstone(fluid="gas", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                    elif magicnumber == 1:
                                        results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                    elif magicnumber == 2:
                                        results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                elif list_rocks[index-1] == "sandstone" and results_subsurface[index-1][i]["fluid"] == "gas":
                                    results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                elif list_rocks[index-1] == "sandstone" and results_subsurface[index-1][i]["fluid"] == "oil":
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                else:
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3))
                                #
                                results_subsurface[index][i]["thickness"] = d
                            else:
                                if list_rocks[index-1] == "shale":
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
                                elif list_rocks[index-1] == "sandstone" and results_subsurface[index-1][i]["fluid"] == "gas":
                                    results_subsurface[index][i] = sandstone(fluid="oil", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                     amounts=results_subsurface[index][i-1]["mineralogy"])
                                elif list_rocks[index-1] == "sandstone" and results_subsurface[index-1][i]["fluid"] == "oil":
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                       amounts=results_subsurface[index][i-1]["mineralogy"])
                                else:
                                    results_subsurface[index][i] = sandstone(fluid="water", actualThickness=0).create_simple_sandstone(dict_output=True, porosity=rd.uniform(0.05, 0.3),
                                                                                                                                       amounts=results_subsurface[index][i-1]["mineralogy"])
                                #
                                results_subsurface[index][i]["thickness"] = d
                #
                elif rock == "shale":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=10)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(0, 0.1))
                        else:
                            results_subsurface[index][i] = shale(fluid="water").create_simple_shale(dict_output=True, porosity=rd.uniform(0, 0.1), amounts=results_subsurface[index][i-1]["mineralogy"])
                #
                elif rock == "granite":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=10)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_granite()
                        else:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_granite(amounts=results_subsurface[index][i-1]["mineralogy"])
                elif rock == "gabbro":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=10)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_gabbro()
                        else:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_gabbro(amounts=results_subsurface[index][i-1]["mineralogy"])
                elif rock == "diorite":
                    thickness_unit = self.split_thickness(thickness=list_thickness[index], n_units=10)
                    for i, d in enumerate(thickness_unit, start=0):
                        if i == 0:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_diorite()
                        else:
                            results_subsurface[index][i] = Plutonic(fluid="water", actualThickness=0, dict_output=True, porosity=rd.uniform(0, 0.05)).create_simple_diorite(amounts=results_subsurface[index][i-1]["mineralogy"])
                #
                results_subsurface[index][i]["thickness"] = list_thickness[index]
        #
        print("Summary:")
        for item in results_subsurface.values():
            for item_part in item.values():
                print(item_part["rock"], item_part["rho"])

        #
    def split_thickness(self, thickness, n_units):
        list_thickness = np.random.multinomial(thickness, np.ones(n_units)/n_units, size=1)[0]
        #
        return list_thickness
    #
if __name__ == "__main__":
    root = tk.Tk()
    GebPyGUI(root)
    root.mainloop()