#!/usr/bin/env python
# -*-coding: utf-8 -*-
# ----------------------
# gui_elements.py
# Maximilian Beeskow
# 22.03.2023
# ----------------------
#
## MODULES
import tkinter as tk
from tkinter import ttk
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
    def create_frame(self, relief=tk.FLAT):
        frm = tk.Frame(self.parent, bg=self.bg, borderwidth=0, highlightthickness=0, relief=relief)
        frm.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns, sticky="nesw")
        #
        return frm
    #
    def create_label(self, text, font_option=None, relief=tk.GROOVE, sticky_option="nesw", anchor_option=tk.CENTER):
        if font_option != None:
            lbl = tk.Label(self.parent, text=text, relief=relief, bg=self.bg, fg=self.fg, font=(font_option),
                           anchor=anchor_option)
        else:
            lbl = tk.Label(self.parent, text=text, relief=relief, bg=self.bg, fg=self.fg, font=(font_option),
                           anchor=anchor_option)
        #
        lbl.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                 sticky=sticky_option)
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
    def create_entry(self, var_entr, var_entr_set=None, width=None, command=None):
        if var_entr_set != None:
            var_entr.set(var_entr_set)
        if width == None:
            entry = tk.Entry(self.parent, textvariable=var_entr, background=self.bg, fg=self.fg, highlightthickness=0)
        else:
            entry = tk.Entry(self.parent, textvariable=var_entr, background=self.bg, fg=self.fg, highlightthickness=0,
                             width=width)
        entry.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                   sticky="nesw")
        if command != None:
            entry.bind("<Return>", command)
        #
        return entry
    #
    def create_radiobutton(self, var_rb, value_rb, color_bg, relief=tk.FLAT, command=None, text="", anchor_option=tk.W):
        if command == None:
            rb = tk.Radiobutton(self.parent, text=text, variable=var_rb, value=value_rb, bg=self.bg, fg=self.fg,
                                activebackground=self.bg, highlightthickness=0, relief=relief, selectcolor=self.bg,
                                anchor=anchor_option)
        else:
            rb = tk.Radiobutton(self.parent, text=text, variable=var_rb, value=value_rb, bg=self.bg, fg=self.fg,
                                activebackground=self.bg, highlightthickness=0, relief=relief, selectcolor=self.bg,
                                anchor=anchor_option, command=command)
        #
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
    #
    def create_treeview(self, n_categories=2, text_1="Ratio", text_2="\u03BC", width_1 = 90, width_2 = 90, text_n=[],
                        width_n=[], individual=False):
        ttk.Style().configure("Treeview", background=self.bg, foreground=self.fg, fieldbackground=self.bg)
        style = ttk.Style()
        style.configure("Treeview.Heading", background=self.bg, pressed_color=self.bg,
                        highlight_color=self.bg, foreground=self.fg)
        #
        if n_categories == 2 and individual == False:
            columns = ("#1", "#2")
            treeview = ttk.Treeview(self.parent, columns=columns, show="headings")
            treeview.heading("#1", text=text_1)
            treeview.column("#1", minwidth=0, width=width_1, stretch=tk.NO, anchor=tk.CENTER)
            treeview.heading("#2", text=text_2)
            treeview.column("#2", minwidth=0, width=width_2, stretch=tk.YES, anchor=tk.CENTER)
            treeview.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                          sticky="nesw")
        #
        if n_categories > 1 and individual == True:
            columns = []
            for n in range(n_categories):
                var_i = "#" + str(n + 1)
                columns.append(var_i)
            #
            treeview = ttk.Treeview(self.parent, columns=columns, show="headings")
            #
            for index, element in enumerate(columns):
                treeview.heading(element, text=text_n[index])
                treeview.column(element, minwidth=0, width=width_n[index], stretch=tk.NO, anchor=tk.CENTER)
            #
            treeview.grid(row=self.row_id, column=self.column_id, rowspan=self.n_rows, columnspan=self.n_columns,
                          sticky="nesw")
        #
        return treeview