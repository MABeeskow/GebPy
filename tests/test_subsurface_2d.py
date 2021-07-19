#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_subsurface_2d.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.07.2021

# -----------------------------------------------

## MODULES
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from modules import subsurface_2d

## CLASS
class TestingSubsurface:
    #
    def __init__(self, coordinates, degree, kind=None):
        self.coordinates = coordinates
        self.degree = degree

        x_left = self.coordinates[0][0]
        x_right = self.coordinates[-1][0]
        x = np.linspace(x_left, x_right, 50+1)
        condition = False
        while condition == False:
            x_outcrops = np.sort(rd.sample(range(int(x_left*(1-0.05)), int(x_right*(1-0.05))), 10))
            x_diff = np.diff(x_outcrops)
            check = any(item in np.arange(0, 5) for item in x_diff)
            if check == False:
                condition = True
            else:
                continue

        fig, ax = plt.subplots(dpi=100)

        if self.degree == "linear":
            a, b = subsurface_2d.Surface(coordinates=self.coordinates).create_linear_surface()
            y = a*x + b
            y_outcrops = a*x_outcrops + b

            ax.plot(x, y, label="Surface", color="C0", zorder=1)
            ax.scatter(x_outcrops, y_outcrops, label="Outcrops", zorder=2, color="C1")
        elif self.degree == "quadratic":
            a, b, c = subsurface_2d.Surface(coordinates=self.coordinates).create_quadratic_surface()
            y = a*x**2 + b*x + c
            y_outcrops = a*x_outcrops**2 + b*x_outcrops + c

            ax.plot(x, y, label="Surface", color="C0", zorder=1)
            ax.scatter(x_outcrops, y_outcrops, label="Outcrops", zorder=2, color="C1")
        elif self.degree == "interpolated":
            f = subsurface_2d.Surface(coordinates=self.coordinates).create_interpolated_surface()
            y = f(x)
            if kind == None:
                x_units, y_units, x_outcrops = subsurface_2d.Units(x_limits=[x[0], x[-1]], y_limits=[y[0], y[-1]],
                                                               number=10).create_units(f_surface=f)
            else:
                x_units, y_units, x_outcrops = subsurface_2d.Units(x_limits=[x[0], x[-1]], y_limits=[y[0], y[-1]],
                                                               number=10).create_units(f_surface=f, kind=kind)
            y_outcrops = f(x_outcrops)

            ax.plot(x, y, label="Surface", color="black", zorder=2)
            ax.scatter(x_outcrops, y_outcrops, label="Outcrops", zorder=3, color="C1")

            for i in range(len(y_units)):
                ax.plot(x_units[i], y_units[i], linestyle="dashed", zorder=1,
                        label=format("Top Unit "+str(i+1)))

        ax.grid()
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_xticks(np.arange(x_left, x_right+10, 20))
        ax.legend(loc="upper left", prop={'size': 8})
        ax.set_axisbelow(True)
        plt.show()

# PROGRAM
TestingSubsurface(coordinates=[[-100, 0], [100, 10]], degree="linear")
TestingSubsurface(coordinates=[[-100, 0], [0, 20], [100, 0]], degree="quadratic")
TestingSubsurface(coordinates=[[-100, 0], [-10, 20], [20, 25], [100, 5]], degree="interpolated")
TestingSubsurface(coordinates=[[-100, 0], [100, 10]], degree="interpolated")
TestingSubsurface(coordinates=[[-100, 0], [-10, 20], [20, 25], [100, 5]], degree="interpolated", kind="horizontal")