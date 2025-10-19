#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		common.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		19.10.2025

#-----------------------------------------------

"""
Module: common.py
This module contains several routines that are commonly used by the different mineral-related modules.
"""

# PACKAGES

# MODULES
from modules.minerals import CrystalPhysics

# CODE
class GeophysicalProperties:
    def __init__(self):
        pass

    def calculate_elastic_properties(self, bulk_mod, shear_mod):
        E = (9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod)           # Young's modulus
        nu = (3*bulk_mod - 2*shear_mod)/(2*(3*bulk_mod + shear_mod))  # Poisson's ratio
        return E, nu

    def calculate_seismic_velocities(self, bulk_mod, shear_mod, rho):
        vPvS = ((bulk_mod + 4/3*shear_mod)/shear_mod)**0.5  # vP/vS
        vP = ((bulk_mod + 4/3*shear_mod)/rho)**0.5          # P-wave velocity vP
        vS = (shear_mod/rho)**0.5                           # S-wave velocity vS
        return vPvS, vP, vS

    def calculate_radiation_properties(self, constr_radiation, rho_electron):
        gamma_ray = constr_radiation.calculate_gr() # Gamma ray
        pe = constr_radiation.calculate_pe()        # Photoelectricity Pe
        U = pe*rho_electron*10**(-3)                # Photoelectricity U
        return gamma_ray, pe, U

class CrystallographicProperties:
    def __init__(self):
        pass

    def calculate_molar_volume(self, constr_volume, constr_molar_volume, cell_z):
        dataV = constr_volume
        V = dataV.calculate_volume()                                                # Cell volume
        V_m = constr_molar_volume.calculate_molar_volume(volume_cell=V, z=cell_z)   # Molar volume
        return V, V_m

    def calculate_mineral_density(self, constr_density):
        rho = constr_density.calculate_bulk_density()   # Density
        return rho

    def calculate_electron_density(self, constr_electron_density):
        rho_e = constr_electron_density.calculate_electron_density()    # Electron density
        return rho_e