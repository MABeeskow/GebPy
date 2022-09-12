#!/usr/bin/env python
# -*-coding: utf-8 -*-
# ----------------------
# gui_elements.py
# Maximilian Beeskow
# 31.08.2022
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
        opt_menu.config(bg=self.bg, fg=self.fg, activebackground=active_bg, highlightthickness=0)
        opt_menu["menu"].config(bg=self.bg, fg=self.fg, activebackground=active_bg)
        #
        return opt_menu
    #
    def create_entry(self, var_entr, var_entr_set, width=None, command=None):
        var_entr.set(var_entr_set)
        if width == None:
            entry = tk.Entry(self.parent, textvariable=var_entr, background=self.bg, highlightthickness=0)
        else:
            entry = tk.Entry(self.parent, textvariable=var_entr, background=self.bg, highlightthickness=0, width=width)
        entry.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                   sticky="nesw")
        if command != None:
            entry.bind("<Return>", command)
        #
        return entry
    #
    def create_radiobutton(self, var_rb, value_rb, color_bg, relief=tk.FLAT, command=None, text=""):
        if command == None:
            rb = tk.Radiobutton(self.parent, text=text, variable=var_rb, value=value_rb, bg=self.bg, fg=self.fg,
                                activebackground=color_bg, highlightthickness=0, relief=relief, selectcolor=color_bg)
        else:
            rb = tk.Radiobutton(self.parent, text=text, variable=var_rb, value=value_rb, bg=self.bg, fg=self.fg,
                                activebackground=color_bg, highlightthickness=0, relief=relief, selectcolor=color_bg,
                                command=command)
        rb.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return  rb
    #
    def create_button(self, text, command=None):
        if command == None:
            btn = tk.Button(
                self.parent, text=text, bg=self.bg, fg=self.fg, activebackground="#F9DED7", highlightbackground=self.bg,
                highlightthickness=0)
        else:
            btn = tk.Button(
                self.parent, text=text, bg=self.bg, fg=self.fg, activebackground="#F9DED7", highlightbackground=self.bg,
                highlightthickness=0, command=command)
        #
        btn.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return btn
    #
    def create_checkbox(self, text, var_cb, command=None):
        if command == None:
            cb = tk.Checkbutton(self.parent, text=text, variable=var_cb, bg=self.bg, fg=self.fg,
                                activebackground="#F9DED7", highlightbackground=self.bg, highlightthickness=0)
        else:
            cb = tk.Checkbutton(self.parent, text=text, variable=var_cb, bg=self.bg, fg=self.fg,
                                activebackground="#F9DED7", highlightbackground=self.bg, highlightthickness=0,
                                command=command)
        cb.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return cb