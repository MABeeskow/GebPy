#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------------------------------------------------------------------------------

# Name:		main
# Author:	Maximilian A. Beeskow
# Version:	pre-release
# Date:		10.12.2024

# -----------------------------------------------------------------------------------------------------------------------

## MODULES
# external
import os, sys
import tkinter as tk
# internal
from gebpy_app import GebPyGUI

def pysills():
	root = tk.Tk()
	root.title("GebPy - Synthetic mineral and rock data")
	path = os.path.dirname(os.path.realpath(sys.argv[0]))
	screen_width = root.winfo_screenwidth()
	screen_height = root.winfo_screenheight()

	GebPyGUI(parent=root, var_screen_width=screen_width, var_screen_height=screen_height, var_path=path)

	root.mainloop()