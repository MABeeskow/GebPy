#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		common.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		16.12.2025

#-----------------------------------------------

"""
Module: common.py
This module contains common functions calculating geophysical parameters.
"""

# PACKAGES
import numpy as np

# MODULES

# CLASSES
class Geophysics:
    def __init__(self):
        pass

    def calculate_seismic_velocities(self, val_K, val_G, val_rho):
        K = val_K*10**9
        G = val_G*10**9
        rho = val_rho
        vS = (G/rho)**0.5
        vP = ((K + 4/3*G)/rho)**0.5
        vPvS = vP/vS

        return vP, vS, vPvS

    def calculate_elastic_parameter_data(self, val_K, val_G):
        K = val_K*10**9
        G = val_G*10**9
        E = ((9*K*G)/(3*K + G))*1e-9
        poisson = (3*K - 2*G)/(2*(3*K + G))
        lame = (K - (2*G)/3)*1e-9

        return E, poisson, lame

    def calculate_bulk_density_data(self, v_phi, val_rho, val_rho_f, val_n):
        porosity = v_phi
        rho_s = val_rho
        rho = (1 - porosity)*val_rho + porosity*val_rho_f
        rho_f = np.ones(val_n)*val_rho_f

        return rho, rho_s, rho_f