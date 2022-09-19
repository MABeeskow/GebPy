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
        gui_elements = ["Frame", "Label", "Button", "Radiobutton", "Checkbox", "Entry", "Option Menu"]
        gui_priority = ["Static", "Temporary"]
        gui_subwindows = ["Trace Elements"]
        for subwindow in gui_subwindows:
            self.gui_elements_sub[subwindow] = {}
            for priority in gui_priority:
                self.gui_elements_sub[subwindow][priority] = {}
                for gui_element in gui_elements:
                    self.gui_elements_sub[subwindow][priority][gui_element] = []
        for priority in gui_priority:
            self.gui_elements[priority] = {}
            for gui_element in gui_elements:
                self.gui_elements[priority][gui_element] = []
        #
        ### Colors
        self.colors_gebpy = {"Background": "#FFFCF2", "Navigation": "#252422", "Accent": "#EB5E28", "Option": "#CCC5B9",
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
            "Cyclosilicates", "Inosilicates", "Phyllosilicates"]
        mineral_groups.sort()
        for mineral_group in mineral_groups:
            if mineral_group == "Oxides":
                sub_oxides = tk.Menu(sub_mineral_groups, tearoff=0)
                self.oxide_minerals = ["Quartz"]
                self.oxide_minerals.sort()
                for mineral in self.oxide_minerals:
                    sub_oxides.add_command(
                        label=mineral, command=lambda name=mineral: self.select_mineral(name))
                #
                sub_mineral_groups.add_cascade(
                    label="Oxides",
                    menu=sub_oxides)
            else:
                sub_mineral_groups.add_command(
                    label=mineral_group)
        #
        sub_special_groups = tk.Menu(mineralogy_menu, tearoff=0)
        special_groups = [
            "Spinel Group", "Hematite Group"]
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
        lbl_analysis = SimpleElements(
            parent=self.parent, row_id=16, column_id=0, n_rows=4, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Analysis Mode", font_option="sans 10 bold", relief=tk.FLAT)
        lbl_traces = SimpleElements(
            parent=self.parent, row_id=20, column_id=0, n_rows=6, n_columns=16, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_label(
            text="Trace Elements", font_option="sans 10 bold", relief=tk.FLAT)
        #
        self.gui_elements["Static"]["Label"].extend(
            [lbl_mineralogy, lbl_name, lbl_samples, lbl_analysis, lbl_traces])
        #
        ## Entries
        entr_samples = SimpleElements(
            parent=self.parent, row_id=14, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Background"],
            fg=self.colors_gebpy["Navigation"]).create_entry(var_entr=self.gui_variables["Entry"]["Number Samples"])
        #
        self.gui_elements["Static"]["Entry"].extend([entr_samples])
        #
        ## Radiobuttons
        rb_geophysics = SimpleElements(
            parent=self.parent, row_id=16, column_id=16, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Physics", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=0,
            color_bg=self.colors_gebpy["Background"], command=self.change_rb_analysis)
        rb_geochemistry = SimpleElements(
            parent=self.parent, row_id=18, column_id=16, n_rows=2, n_columns=14, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Mineral Chemistry", var_rb=self.gui_variables["Radiobutton"]["Analysis Mode"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=self.change_rb_analysis)
        rb_trace_without = SimpleElements(
            parent=self.parent, row_id=20, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="Without Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=0,
            color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
        rb_trace_with = SimpleElements(
            parent=self.parent, row_id=22, column_id=16, n_rows=2, n_columns=15, bg=self.colors_gebpy["Navigation"],
            fg=self.colors_gebpy["Background"]).create_radiobutton(
            text="With Trace Elements", var_rb=self.gui_variables["Radiobutton"]["Trace Elements"], value_rb=1,
            color_bg=self.colors_gebpy["Navigation"], command=lambda var_name=name: self.change_rb_traces(var_name))
        #
        self.gui_elements["Static"]["Radiobutton"].extend(
            [rb_geophysics, rb_geochemistry, rb_trace_without, rb_trace_with])
        #
        ## Buttons
        btn_generate_data = SimpleElements(
            parent=self.parent, row_id=28, column_id=16, n_rows=4, n_columns=15,
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
                    for gui_item in gui_items:
                        gui_item.grid_remove()
                    gui_items.clear()
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
                "M\n (mol)", "V\n (A3/mol)", "rho\n (kg/m3)", "vP\n (m/s)", "vS\n (m/s)", "vP/vS\n (1)", "K\n (GPa)",
                "G\n (GPa)", "E\n (GPa)", "nu\n (1)", "GR\n (API)", "PE\n (barns/e^-)"]
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
                entr_min = SimpleElements(
                    parent=self.parent, row_id=2*(2*index + 4), column_id=start_column + 9, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Minimum"][categories_short[index]])
                entr_max = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 18, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Maximum"][categories_short[index]])
                entr_mean = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 27, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mean"][categories_short[index]])
                entr_error = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 36, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Error"][categories_short[index]])
                #
                self.gui_elements["Temporary"]["Entry"].extend([entr_min, entr_max, entr_mean, entr_error])
                #
        elif self.gui_variables["Radiobutton"]["Analysis Mode"].get() == 1:   # Mineral Chemistry
            ## Cleaning
            for key, gui_items in self.gui_elements["Temporary"].items():
                if len(gui_items) > 0:
                    for gui_item in gui_items:
                        gui_item.grid_remove()
                    gui_items.clear()
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
                entr_min = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 9, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Minimum"][element])
                entr_max = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 18, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Maximum"][element])
                entr_mean = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 27, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Mean"][element])
                entr_error = SimpleElements(
                    parent=self.parent, row_id=2 * (2 * index + 4), column_id=start_column + 36, n_rows=4, n_columns=9,
                    bg=self.colors_gebpy["White"], fg=self.colors_gebpy["Navigation"]).create_entry(
                    var_entr=self.gui_variables["Entry"]["Error"][element])
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
                    parent=self.parent, row_id=24, column_id=16, n_rows=2, n_columns=15,
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
            data_mineral = Oxides(
                mineral=var_name, data_type=True, traces_list=self.traces_list).generate_dataset(
                number=self.gui_variables["Entry"]["Number Samples"].get())
        for key, dataset in data_mineral.items():
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
    GebPyGUI(parent=root)
    root.mainloop()