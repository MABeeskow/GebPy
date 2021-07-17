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

        if self.degree == "linear":
            a, b = subsurface_2d.Surface(coordinates=self.coordinates).create_linear_surface()
            y = a*x + b
        elif self.degree == "quadratic":
            a, b, c = subsurface_2d.Surface(coordinates=self.coordinates).create_quadratic_surface()
            y = a*x**2 + b*x + c
        elif self.degree == "interpolated":
            f = subsurface_2d.Surface(coordinates=self.coordinates).create_interpolated_surface()
            y = f(x)

        fig, ax = plt.subplots(dpi=100)
        ax.scatter(x, y)
        ax.plot(x, y)
        ax.grid()
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_xticks(np.arange(x_left, x_right+10, 10))
        plt.show()

# PROGRAM
TestingSubsurface(coordinates=[[-50, 0], [50, 10]], degree="linear")
TestingSubsurface(coordinates=[[-50, 0], [0, 20], [50, 0]], degree="quadratic")
TestingSubsurface(coordinates=[[-50, 0], [-10, 20], [20, 25], [50, 5]], degree="interpolated")
