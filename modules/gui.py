#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		gui.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		30.10.2021

#-----------------------------------------------

## MODULES
import tkinter as tk
from modules.gui_elements import SimpleElements as SE
from modules.oxides import Oxides
from modules.silicates import Tectosilicates
from modules.minerals import feldspars
from modules.siliciclastics import sandstone, shale
from modules.carbonates import limestone, dolomite
from modules.igneous import Plutonic
from modules.evaporites import Evaporites
from modules.sequences import DataProcessing as DP
import numpy as np
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
        self.color_bg = "#E9ECED"
        self.color_fg = "black"
        #
        self.parent = parent
        self.parent.title("GebPy")
        self.parent.geometry("1120x920")
        self.parent["bg"] = self.color_bg
        #
        self.var_rb_mode = tk.IntVar()
        var_rb_mode_start = 0
        #
        for x in range(12):
            tk.Grid.columnconfigure(self.parent, x, weight=1)
        for y in range(26):
            tk.Grid.rowconfigure(self.parent, y, weight=1)
        #
        # Rows
        for i in range(0, 20):
            self.parent.grid_rowconfigure(i, minsize=30)
        for i in range(20, 25):
            self.parent.grid_rowconfigure(i, minsize=60)
        self.parent.grid_rowconfigure(25, minsize=20)
        # Columns
        self.parent.grid_columnconfigure(0, minsize=180)
        self.parent.grid_columnconfigure(1, minsize=15)
        for i in range(2, 11):
            self.parent.grid_columnconfigure(i, minsize=100)
        self.parent.grid_columnconfigure(11, minsize=15)
        #
        menu_frame = tk.Frame(self.parent, bg="#264653", borderwidth=0, highlightthickness=0)
        menu_frame.grid(row=0, column=0, rowspan=26, sticky="nesw")
        menu_frame2 = tk.Frame(self.parent, bg="#7C9097", borderwidth=0, highlightthickness=0)
        menu_frame2.grid(row=0, column=1, rowspan=26, sticky="nesw")
        #
        ## Logo
        gebpy_logo = tk.PhotoImage(file="../documents/readme_images/GebPy_Logo.png")
        gebpy_logo = gebpy_logo.subsample(6, 6)
        img = tk.Label(self.parent, image=gebpy_logo, bg="#264653")
        img.image = gebpy_logo
        img.grid(row=0, column=0, rowspan=2, columnspan=1, sticky="nesw")
        #
        ## Radiobuttons
        SE(parent=self.parent, row_id=2, column_id=0, bg=self.color_bg, fg="white").create_radiobutton(
            var_rb=self.var_rb_mode, var_rb_set=var_rb_mode_start, value_rb=0, text="Minerals", color_bg="#264653",
            command=lambda var_rb_mode=self.var_rb_mode: self.change_radiobutton_mode(var_rb_mode))
        SE(parent=self.parent, row_id=5, column_id=0, bg=self.color_bg, fg="white").create_radiobutton(
            var_rb=self.var_rb_mode, var_rb_set=var_rb_mode_start, value_rb=1, text="Rocks", color_bg="#264653",
            command=lambda var_rb_mode=self.var_rb_mode: self.change_radiobutton_mode(var_rb_mode))
        SE(parent=self.parent, row_id=8, column_id=0, bg=self.color_bg, fg="white").create_radiobutton(
            var_rb=self.var_rb_mode, var_rb_set=var_rb_mode_start, value_rb=2, text="Sequences", color_bg="#264653",
            command=lambda var_rb_mode=self.var_rb_mode: self.change_radiobutton_mode(var_rb_mode))
        #
        ## Option Menu
        var_opt_0_0 = tk.StringVar()
        opt_list_0_0 = ["Oxides", "Sulfides", "Carbonates", "Halogenes", "Tectosilicates"]
        self.opt_mingroup = SE(parent=self.parent, row_id=3, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
            var_opt=var_opt_0_0, var_opt_set="Select Mineral Group", opt_list=opt_list_0_0,
            command=lambda var_opt=var_opt_0_0: self.select_opt(var_opt))
        var_opt_1_0 = tk.StringVar()
        opt_list_1_0 = ["Siliciclastic Rocks", "Carbonate Rocks", "Igneous Rocks", "Metamorphic Rocks",
                        "Evaporite Rocks"]
        self.opt_rocktype = SE(parent=self.parent, row_id=6, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
            var_opt=var_opt_1_0, var_opt_set="Select Rock Type", opt_list=opt_list_1_0,
            command=lambda var_opt=var_opt_1_0: self.select_opt(var_opt))
        var_opt_2_0 = tk.StringVar()
        opt_list_2_0 = ["Zechstein", "Buntsandstein"]
        self.opt_realseq = SE(parent=self.parent, row_id=9, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
            var_opt=var_opt_2_0, var_opt_set="Select Real Sequences", opt_list=opt_list_2_0,
            command=lambda var_opt=var_opt_2_0: self.select_opt(var_opt))
        #
        ## Button
        self.btn_randseq = SE(parent=self.parent, row_id=10, column_id=0, bg="#264653",
                              fg=self.color_fg).create_button(text="Create Random Sequence")
        self.btn_custseq = SE(parent=self.parent, row_id=11, column_id=0, bg="#264653",
                              fg=self.color_fg).create_button(text="Create Custom Sequence")
    #
    def select_opt(self, var_opt):
        # Minerals
        if var_opt == "Quartz":
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, mineral=var_opt)
        elif var_opt == "Alkalifeldspar":
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, mineral=var_opt)
        elif var_opt == "Plagioclase":
            Minerals(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, mineral=var_opt)
        # Rocks
        elif var_opt == "Sandstone":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Shale":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        #
        elif var_opt == "Limestone":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Dolomite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        #
        elif var_opt == "Felsic Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Intermediate Rock":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Granite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Gabbro":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Diorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Granodiorite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Syenite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Tonalite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Monzonite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Quartzolite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Qz-rich Granitoid":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        #
        elif var_opt == "Halite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        elif var_opt == "Anhydrite":
            Rocks(parent=self.parent, color_bg=self.color_bg, color_fg=self.color_fg, rock=var_opt)
        #
        elif var_opt == "Oxides":
            try:
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
            except:
                pass
            var_opt_0_1 = tk.StringVar()
            opt_list_0_1 = ["Quartz", "Magnetite", "Hematite"]
            self.opt_oxide = SE(parent=self.parent, row_id=4, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_0_1, var_opt_set="Select Oxide", opt_list=opt_list_0_1,
                command=lambda var_opt=var_opt_0_1: self.select_opt(var_opt))
        elif var_opt == "Sulfides":
            try:
                self.opt_oxide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
            except:
                pass
            var_opt_0_2 = tk.StringVar()
            opt_list_0_2 = ["Pyrite", "Chalcopyrite", "Galena"]
            self.opt_sulfide = SE(parent=self.parent, row_id=4, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_0_2, var_opt_set="Select Sulfide", opt_list=opt_list_0_2,
                command=lambda var_opt=var_opt_0_2: self.select_opt(var_opt))
        elif var_opt == "Carbonates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_halogene.grid_remove()
            except:
                pass
            var_opt_0_3 = tk.StringVar()
            opt_list_0_3 = ["Calcite", "Dolomite", "Magnesite"]
            self.opt_carb = SE(parent=self.parent, row_id=4, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_0_3, var_opt_set="Select Carbonate", opt_list=opt_list_0_3,
                command=lambda var_opt=var_opt_0_3: self.select_opt(var_opt))
        elif var_opt == "Halogenes":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
            except:
                pass
            var_opt_0_4 = tk.StringVar()
            opt_list_0_4 = ["Halite", "Fluorite", "Sylvite"]
            self.opt_halogene = SE(parent=self.parent, row_id=4, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_0_4, var_opt_set="Select Halogene", opt_list=opt_list_0_4,
                command=lambda var_opt=var_opt_0_4: self.select_opt(var_opt))
        elif var_opt == "Tectosilicates":
            try:
                self.opt_oxide.grid_remove()
                self.opt_sulfide.grid_remove()
                self.opt_carb.grid_remove()
                self.opt_halogene.grid_remove()
            except:
                pass
            var_opt_0_5 = tk.StringVar()
            opt_list_0_5 = ["Alkalifeldspar", "Plagioclase"]
            self.opt_afs = SE(parent=self.parent, row_id=4, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_0_5, var_opt_set="Select Tectosilicate", opt_list=opt_list_0_5,
                command=lambda var_opt=var_opt_0_5: self.select_opt(var_opt))
        elif var_opt == "Siliciclastic Rocks":
            var_opt_1_1 = tk.StringVar()
            opt_list_1_1 = ["Sandstone", "Shale"]
            self.opt_silic = SE(parent=self.parent, row_id=7, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_1_1, var_opt_set="Select Rock", opt_list=opt_list_1_1,
                command=lambda var_opt=var_opt_1_1: self.select_opt(var_opt))
        elif var_opt == "Carbonate Rocks":
            var_opt_1_2 = tk.StringVar()
            opt_list_1_2 = ["Limestone", "Dolomite"]
            self.opt_carb = SE(parent=self.parent, row_id=7, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_1_2, var_opt_set="Select Rock", opt_list=opt_list_1_2,
                command=lambda var_opt=var_opt_1_2: self.select_opt(var_opt))
        elif var_opt == "Igneous Rocks":
            var_opt_1_3 = tk.StringVar()
            opt_list_1_3 = ["Felsic Rock", "Intermediate Rock", "Granite", "Gabbro", "Diorite", "Granodiorite",
                            "Monzonite", "Syenite", "Tonalite", "Quartzolite", "Qz-rich Granitoid"]
            self.opt_ign = SE(parent=self.parent, row_id=7, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_1_3, var_opt_set="Select Rock", opt_list=opt_list_1_3,
                command=lambda var_opt=var_opt_1_3: self.select_opt(var_opt))
        elif var_opt == "Evaporite Rocks":
            var_opt_1_4 = tk.StringVar()
            opt_list_1_4 = ["Halite", "Anhydrite"]
            self.opt_ign = SE(parent=self.parent, row_id=7, column_id=0, bg="#264653", fg=self.color_fg).create_option_menu(
                var_opt=var_opt_1_4, var_opt_set="Select Rock", opt_list=opt_list_1_4,
                command=lambda var_opt=var_opt_1_4: self.select_opt(var_opt))
    #
    def change_radiobutton_mode(self, var_rb_mode):
        if var_rb_mode.get() == 0:
            # try:
            #     self.opt_rocktype.grid_remove()
            #     self.opt_realseq.grid_remove()
            #     self.btn_randseq.grid_remove()
            #     self.btn_custseq.grid_remove()
            #     self.opt_oxide.grid_remove()
            #     self.opt_sulfide.grid_remove()
            #     self.opt_carb.grid_remove()
            #     self.opt_halogene.grid_remove()
            # except:
            #     pass
            print(var_rb_mode.get())
        elif var_rb_mode.get() == 1:
            print(var_rb_mode.get())
        elif var_rb_mode.get() == 2:
            print(var_rb_mode.get())
#
class Minerals:
    #
    def __init__(self, parent, color_bg, color_fg, mineral):
        self.parent_mineral = parent
        self.color_bg =color_bg
        self.color_fg = color_fg
        self.mineral = mineral
        self.var_rb = tk.IntVar()
        var_rb_start = 0
        self.var_entr = tk.IntVar()
        var_entr_start = 100
        #
        SE(parent=self.parent_mineral, row_id=15, column_id=0, bg=self.color_bg,
           fg="white").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=0, text="Histogram",
                                           color_bg="#264653",
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        SE(parent=self.parent_mineral, row_id=16, column_id=0, bg=self.color_bg,
           fg="white").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=1, text="Scatter",
                                           color_bg="#264653",
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        #
        SE(parent=self.parent_mineral, row_id=13, column_id=0, bg="#264653",
           fg="white").create_label(text="Number of samples", relief=tk.RAISED)
        SE(parent=self.parent_mineral, row_id=14, column_id=0, bg=self.color_bg,
           fg=color_fg).create_entry(var_entr=self.var_entr, var_entr_set=var_entr_start, command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        #
        data_all = []
        for i in range(var_entr_start):
            if self.mineral == "Quartz":
                data = Oxides(impurity="pure").create_quartz()
            elif self.mineral == "Alkalifeldspar":
                data = Tectosilicates(impurity="pure").create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data = Tectosilicates(impurity="pure").create_plagioclase()
            self.color_mineral = "#7C9097"
            #
            data_all.append(data)
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
        #
        self.create_plot(parent=self.parent_mineral, data=self.rho_b/1000, row_id=0, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vP/1000, row_id=0, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vS/1000, row_id=0, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.bulk_mod, row_id=10, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.shear_mod, row_id=10, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.poisson, row_id=10, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.molar_mass, row_id=20, column_id=2, n_rows=5,
                         n_columns=3, xlabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.gamma_ray, row_id=20, column_id=5, n_rows=5,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.photoelectricity, row_id=20, column_id=8, n_rows=5,
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
            elif self.mineral == "Alkalifeldspar":
                data = Tectosilicates(impurity="pure").create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data = Tectosilicates(impurity="pure").create_plagioclase()
            self.color_mineral = "#7C9097"
            #
            data_all.append(data)
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
        #
        self.create_plot(parent=self.parent_mineral, data=self.rho_b/1000, row_id=0, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vP/1000, row_id=0, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.vS/1000, row_id=0, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.bulk_mod, row_id=10, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.shear_mod, row_id=10, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.poisson, row_id=10, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.molar_mass, row_id=20, column_id=2, n_rows=5,
                         n_columns=3, xlabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.gamma_ray, row_id=20, column_id=5, n_rows=5,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_mineral)
        self.create_plot(parent=self.parent_mineral, data=self.photoelectricity, row_id=20, column_id=8, n_rows=5,
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
            self.create_plot(parent=self.parent_mineral, data=self.rho_b/1000, row_id=0, column_id=2, n_rows=10,
                             n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.vP/1000, row_id=0, column_id=5, n_rows=10,
                             n_columns=3, xlabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.vS/1000, row_id=0, column_id=8, n_rows=10,
                             n_columns=3, xlabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.bulk_mod, row_id=10, column_id=2, n_rows=10,
                             n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.shear_mod, row_id=10, column_id=5, n_rows=10,
                             n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.poisson, row_id=10, column_id=8, n_rows=10,
                             n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.molar_mass, row_id=20, column_id=2, n_rows=5,
                             n_columns=3, xlabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.gamma_ray, row_id=20, column_id=5, n_rows=5,
                             n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_mineral)
            self.create_plot(parent=self.parent_mineral, data=self.photoelectricity, row_id=20, column_id=8, n_rows=5,
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
            elif self.mineral == "Alkalifeldspar":
                element = "K"
            elif self.mineral == "Plagioclase":
                element = "Ca"
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.vP/1000, row_id=0,
                                     column_id=2, n_rows=10, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Seismic velocity $v_P$ (km/s)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.vS/1000, row_id=0,
                                     column_id=5, n_rows=10, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Seismic velocity $v_S$ (km/s)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.vP/self.vS,
                                     row_id=0, column_id=8, n_rows=10, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Velocity ratio $v_P/v_S$ (1)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.bulk_mod, row_id=10,
                                     column_id=2, n_rows=10, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Bulk modulus $K$ (GPa)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.shear_mod,
                                     row_id=10, column_id=5, n_rows=10, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Shear modulus $G$ (GPa)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.poisson, row_id=10,
                                     column_id=8, n_rows=10, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Poisson's ratio $\\mu$ (1)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.molar_mass,
                                     row_id=20, column_id=2, n_rows=5, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Molar mass $M$ (g/mol)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.gamma_ray,
                                     row_id=20, column_id=5, n_rows=5, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Gamma ray GR (API)", color=self.color_mineral)
            self.create_scatter_plot(parent=self.parent_mineral, data_x=self.w_element, data_y=self.photoelectricity,
                                     row_id=20, column_id=8, n_rows=5, n_columns=3,
                                     xlabel=str(element)+" amount $w_{"+str(element)+"}$ (1)",
                                     ylabel="Photoelectricity PE (barns/electron)", color=self.color_mineral)
    #
#
class Rocks:
    #
    def __init__(self, parent, color_bg, color_fg, rock):
        self.parent_sandstone = parent
        self.color_bg = color_bg
        self.var_entr = tk.IntVar()
        var_entr_start = 100
        self.var_rb = tk.IntVar()
        var_rb_start = 0
        self.rock = rock
        #
        SE(parent=self.parent_sandstone, row_id=15, column_id=0, bg="#264653",
           fg="white").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=0, text="Histogram",
                                           color_bg="#264653",
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        SE(parent=self.parent_sandstone, row_id=16, column_id=0, bg="#264653",
           fg="white").create_radiobutton(var_rb=self.var_rb, var_rb_set=var_rb_start, value_rb=1, text="Scatter",
                                           color_bg="#264653",
                                           command=lambda var_rb=self.var_rb: self.change_radiobutton(var_rb))
        #
        SE(parent=self.parent_sandstone, row_id=13, column_id=0, bg="#264653",
           fg="white").create_label(text="Number of samples")
        SE(parent=self.parent_sandstone, row_id=14, column_id=0, bg="#264653",
           fg=color_fg).create_entry(var_entr=self.var_entr, var_entr_set=var_entr_start, command=lambda event, var_entr=self.var_entr: self.enter_samples(var_entr, event))
        #
        data_all = []
        for i in range(var_entr_start):
            if self.rock == "Sandstone":
                data = sandstone(fluid="water", actualThickness=0).create_simple_sandstone()
                self.color_rock = "tan"
            elif self.rock == "Shale":
                data = shale().create_simple_shale()
                self.color_rock = "olivedrab"
            elif self.rock == "Limestone":
                data = limestone(fluid="water", actualThickness=0).create_simple_limestone()
                self.color_rock = "lightblue"
            elif self.rock == "Dolomite":
                data = dolomite(fluid="water", actualThickness=0).create_simple_dolomite()
                self.color_rock = "violet"
            #
            elif self.rock == "Felsic Rock":
                data = Plutonic(fluid="water", actualThickness=0).create_felsic()
                self.color_rock = "coral"
            elif self.rock == "Intermediate Rock":
                data = Plutonic(fluid="water", actualThickness=0).create_intermediate()
                self.color_rock = "coral"
            elif self.rock == "Granite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_granite()
                self.color_rock = "coral"
            elif self.rock == "Gabbro":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_gabbro()
                self.color_rock = "coral"
            elif self.rock == "Diorite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_diorite()
                self.color_rock = "coral"
            elif self.rock == "Granodiorite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_granodiorite()
                self.color_rock = "coral"
            elif self.rock == "Monzonite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_monzonite()
                self.color_rock = "coral"
            elif self.rock == "Quartzolite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_quartzolite()
                self.color_rock = "coral"
            elif self.rock == "Qz-rich Granitoid":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_quartzrich_granitoid()
                self.color_rock = "coral"
            elif self.rock == "Syenite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_syenite()
                self.color_rock = "coral"
            elif self.rock == "Tonalite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_tonalite()
                self.color_rock = "coral"
            #
            elif self.rock == "Halite":
                data = Evaporites(fluid="water", actualThickness=0).create_simple_rocksalt()
                self.color_rock = "thistle"
            elif self.rock == "Anhydrite":
                data = Evaporites(fluid="water", actualThickness=0).create_simple_anhydrite()
                self.color_rock = "pink"
            #
            data_all.append(data)
        self.rho_b = DP(dataset=data_all).extract_densities(type="random", keyword="bulk")
        self.vP = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vP")
        self.vS = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vS")
        self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="bulk")
        self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="shear")
        self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="poisson")
        self.phi = DP(dataset=data_all).extract_porosity(type="random")
        self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="random")
        self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="random")
        #
        self.create_plot(parent=self.parent_sandstone, data=self.rho_b, row_id=0, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.vP, row_id=0, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.vS, row_id=0, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.bulk_mod, row_id=10, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.shear_mod, row_id=10, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.poisson, row_id=10, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.phi*100, row_id=20, column_id=2, n_rows=5,
                         n_columns=3, xlabel="Porosity $\\phi$ (%)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.gamma_ray, row_id=20, column_id=5, n_rows=5,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.photoelectricity, row_id=20, column_id=8, n_rows=5,
                         n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
    #
    def create_plot(self, parent, data, row_id, column_id, n_rows, n_columns, xlabel, color):
        #
        self.canvas_histo = None
        self.fig_histo = Figure(facecolor="#E9ECED")
        self.ax_histo = self.fig_histo.add_subplot()
        #
        self.ax_histo.axvline(x=np.mean(data), color="black", linestyle="dashed")
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
        self.fig = Figure(facecolor="#E9ECED")
        self.ax = self.fig.add_subplot()
        #
        self.ax.scatter(data_x, data_y, color=color, edgecolor="black")
        self.ax.grid(True)
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
        data_all = []
        for i in range(var_entr.get()):
            if self.rock == "Sandstone":
                data = sandstone(fluid="water", actualThickness=0).create_simple_sandstone()
                self.color_rock = "tan"
            elif self.rock == "Shale":
                data = shale().create_simple_shale()
                self.color_rock = "olivedrab"
            elif self.rock == "Limestone":
                data = limestone(fluid="water", actualThickness=0).create_simple_limestone()
                self.color_rock = "lightblue"
            elif self.rock == "Dolomite":
                data = dolomite(fluid="water", actualThickness=0).create_simple_dolomite()
                self.color_rock = "violet"
            #
            elif self.rock == "Felsic Rock":
                data = Plutonic(fluid="water", actualThickness=0).create_felsic()
                self.color_rock = "coral"
            elif self.rock == "Intermediate Rock":
                data = Plutonic(fluid="water", actualThickness=0).create_intermediate()
                self.color_rock = "coral"
            elif self.rock == "Granite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_granite()
                self.color_rock = "coral"
            elif self.rock == "Gabbro":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_gabbro()
                self.color_rock = "coral"
            elif self.rock == "Diorite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_diorite()
                self.color_rock = "coral"
            elif self.rock == "Granodiorite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_granodiorite()
                self.color_rock = "coral"
            elif self.rock == "Monzonite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_monzonite()
                self.color_rock = "coral"
            elif self.rock == "Quartzolite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_quartzolite()
                self.color_rock = "coral"
            elif self.rock == "Qz-rich Granitoid":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_quartzrich_granitoid()
                self.color_rock = "coral"
            elif self.rock == "Syenite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_syenite()
                self.color_rock = "coral"
            elif self.rock == "Tonalite":
                data = Plutonic(fluid="water", actualThickness=0).create_simple_tonalite()
                self.color_rock = "coral"
            #
            elif self.rock == "Halite":
                data = Evaporites(fluid="water", actualThickness=0).create_simple_rocksalt()
                self.color_rock = "thistle"
            elif self.rock == "Anhydrite":
                data = Evaporites(fluid="water", actualThickness=0).create_simple_anhydrite()
                self.color_rock = "pink"
            #
            data_all.append(data)
        self.rho_b = DP(dataset=data_all).extract_densities(type="random", keyword="bulk")
        self.vP = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vP")
        self.vS = DP(dataset=data_all).extract_seismic_velocities(type="random", keyword="vS")
        self.bulk_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="bulk")
        self.shear_mod = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="shear")
        self.poisson = DP(dataset=data_all).extract_elastic_moduli(type="random", keyword="poisson")
        self.phi = DP(dataset=data_all).extract_porosity(type="random")
        self.gamma_ray = DP(dataset=data_all).extract_gamma_ray(type="random")
        self.photoelectricity = DP(dataset=data_all).extract_photoelectricity(type="random")
        #
        self.create_plot(parent=self.parent_sandstone, data=self.rho_b, row_id=0, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.vP, row_id=0, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.vS, row_id=0, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.bulk_mod, row_id=10, column_id=2, n_rows=10,
                         n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.shear_mod, row_id=10, column_id=5, n_rows=10,
                         n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.poisson, row_id=10, column_id=8, n_rows=10,
                         n_columns=3, xlabel="Poisson's ratio $\\mu$", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.phi*100, row_id=20, column_id=2, n_rows=5,
                         n_columns=3, xlabel="Porosity $\\phi$ (%)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.gamma_ray, row_id=20, column_id=5, n_rows=5,
                         n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_rock)
        self.create_plot(parent=self.parent_sandstone, data=self.photoelectricity, row_id=20, column_id=8, n_rows=5,
                         n_columns=3, xlabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
        #
        self.var_rb.set(0)
    #
    def change_radiobutton(self, var_rb):
        print(var_rb.get())
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
            self.create_plot(parent=self.parent_sandstone, data=self.rho_b, row_id=0, column_id=2, n_rows=10,
                             n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.vP, row_id=0, column_id=5, n_rows=10,
                             n_columns=3, xlabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.vS, row_id=0, column_id=8, n_rows=10,
                             n_columns=3, xlabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.bulk_mod, row_id=10, column_id=2, n_rows=10,
                             n_columns=3, xlabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.shear_mod, row_id=10, column_id=5, n_rows=10,
                             n_columns=3, xlabel="Shear modulus $G$ (GPa)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.poisson, row_id=10, column_id=8, n_rows=10,
                             n_columns=3, xlabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.phi*100, row_id=20, column_id=2, n_rows=5,
                             n_columns=3, xlabel="Porosity $\\phi$ (%)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.gamma_ray, row_id=20, column_id=5, n_rows=5,
                             n_columns=3, xlabel="Gamma ray GR (API)", color=self.color_rock)
            self.create_plot(parent=self.parent_sandstone, data=self.photoelectricity, row_id=20, column_id=8, n_rows=5,
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
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.vP, row_id=0,
                                     column_id=2, n_rows=10, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Seismic velocity $v_P$ (m/s)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.vS, row_id=0,
                                     column_id=5, n_rows=10, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Seismic velocity $v_S$ (m/s)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.vP/self.vS, row_id=0,
                                     column_id=8, n_rows=10, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Velocity ratio $v_P/v_S$ (1)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.bulk_mod, row_id=10,
                                     column_id=2, n_rows=10, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Bulk modulus $K$ (GPa)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.shear_mod, row_id=10,
                                     column_id=5, n_rows=10, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Shear modulus $G$ (GPa)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.poisson, row_id=10,
                                     column_id=8, n_rows=10, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Poisson's ratio $\\mu$ (1)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.phi*100, row_id=20,
                                     column_id=2, n_rows=5, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Porosity $\\phi$ (%)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.gamma_ray, row_id=20,
                                     column_id=5, n_rows=5, n_columns=3, xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Gamma ray GR (API)", color=self.color_rock)
            self.create_scatter_plot(parent=self.parent_sandstone, data_x=self.rho_b, data_y=self.photoelectricity,
                                     row_id=20, column_id=8, n_rows=5, n_columns=3,
                                     xlabel="Densitiy $\\varrho$ (g/ccm)",
                                     ylabel="Photoelectricity PE (barns/electron)", color=self.color_rock)
    #

if __name__ == "__main__":
    root = tk.Tk()
    GebPyGUI(root)
    root.mainloop()