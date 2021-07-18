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
    def __init__(self, coordinates, degree):
        self.coordinates = coordinates
        self.degree = degree

        x_left = self.coordinates[0][0]
        x_right = self.coordinates[-1][0]
        x = np.linspace(x_left, x_right, 25+1)

        condition = False
        while condition == False:
            x_outcrops = np.sort(rd.sample(range(int(x_left*(1+0.05)), int(x_right*(1-0.05))), 5))
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
            ax.scatter(x, y, label="Surface", zorder=2)
            ax.plot(x, y, zorder=1)
            ax.scatter(x_outcrops, y_outcrops, label="Outcrops", zorder=2)
        elif self.degree == "quadratic":
            a, b, c = subsurface_2d.Surface(coordinates=self.coordinates).create_quadratic_surface()
            y = a*x**2 + b*x + c
            y_outcrops = a*x_outcrops**2 + b*x_outcrops + c
            ax.scatter(x, y, label="Surface", zorder=2)
            ax.plot(x, y, zorder=1)
            ax.scatter(x_outcrops, y_outcrops, label="Outcrops", zorder=2)
        elif self.degree == "interpolated":
            f = subsurface_2d.Surface(coordinates=self.coordinates).create_interpolated_surface()
            subsurface_2d.Units(x_values=x_outcrops).create_units(f_surface=f)
            y = f(x)
            y_outcrops = f(x_outcrops)

            ax.plot(x, y, label="Surface", zorder=1)
            ax.scatter(x_outcrops, y_outcrops, label="Outcrops", zorder=2, color="orangered")

        ax.grid()
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_xticks(np.arange(x_left, x_right+10, 10))
        ax.legend(loc="upper left")
        ax.set_axisbelow(True)
        plt.show()

# PROGRAM
TestingSubsurface(coordinates=[[-50, 0], [50, 10]], degree="linear")
TestingSubsurface(coordinates=[[-50, 0], [0, 20], [50, 0]], degree="quadratic")
TestingSubsurface(coordinates=[[-50, 0], [-10, 20], [20, 25], [50, 5]], degree="interpolated")
TestingSubsurface(coordinates=[[-50, 0], [50, 10]], degree="interpolated")