#!/usr/bin/env python
# -*-coding: utf-8 -*-
# ----------------------
# gui_elements.py
# Maximilian Beeskow
# 01.11.2021
# ----------------------
#
## MODULES
import tkinter as tk
#
## CLASSES
class SimpleElements:
    #
    def __init__(self, parent, row_id, column_id, fg, bg, n_rows=1, n_columns=1):
        self.parent = parent
        self.row_id = row_id
        self.column_id = column_id
        self.n_rows = n_rows
        self.n_columns = n_columns
        self.fg = fg
        self.bg = bg
    #
    def create_label(self, text, fontsize=None, relief=tk.GROOVE):
        if fontsize != None:
            lbl = tk.Label(self.parent, text=text, relief=relief, bg=self.bg, fg=self.fg, font=(fontsize))
        else:
            lbl = tk.Label(self.parent, text=text, relief=relief, bg=self.bg, fg=self.fg)
        lbl.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return lbl
    #
    def create_option_menu(self, var_opt, var_opt_set, opt_list, active_bg="#F9DED7", command=None):
        var_opt.set(var_opt_set)
        if command == None:
            opt_menu = tk.OptionMenu(self.parent, var_opt, *opt_list)
        else:
            opt_menu = tk.OptionMenu(self.parent, var_opt, *opt_list, command=command)
        opt_menu.config(bg=self.bg)
        opt_menu.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                      sticky="nesw")
        opt_menu.config(bg=self.bg, activebackground=active_bg, highlightthickness=0)
        opt_menu["menu"].config(bg=self.bg, activebackground=active_bg)
        #
        return opt_menu
    #
    def create_entry(self, var_entr, var_entr_set, command=None):
        var_entr.set(var_entr_set)
        entry = tk.Entry(self.parent, textvariable=var_entr, background=self.bg, highlightthickness=0)
        entry.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                   sticky="nesw")
        if command != None:
            entry.bind("<Return>", command)
        #
        return entry
    #
    def create_radiobutton(self, var_rb, var_rb_set, value_rb, color_bg, relief=tk.FLAT, command=None, text=""):
        var_rb.set(var_rb_set)
        if command == None:
            rb = tk.Radiobutton(self.parent, text=text, variable=var_rb, value=value_rb, bg=color_bg,
                                activebackground=color_bg, highlightthickness=0, relief=relief)
        else:
            rb = tk.Radiobutton(self.parent, text=text, variable=var_rb, value=value_rb, bg=color_bg, fg=self.fg,
                                activebackground=color_bg, highlightthickness=0, relief=relief, command=command)
        rb.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return  rb
    #
    def create_button(self, text, command=None):
        if command == None:
            btn = tk.Button(self.parent, text=text, bg=self.bg, fg=self.fg, activebackground="#F9DED7", highlightbackground=self.bg)
        else:
            btn = tk.Button(self.parent, text=text, bg=self.bg, fg=self.fg, activebackground="#F9DED7", highlightbackground=self.bg,
                            highlightthickness=0, command=command)
        btn.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return btn
    #