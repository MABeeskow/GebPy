#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		common.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		27.11.2025

#-----------------------------------------------

"""
Module: common.py
This module contains several routines that are commonly used by the different mineral-related modules.
"""

# PACKAGES
import numpy as np
import scipy

# MODULES

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

class CrystalPhysics:
    def __init__(self, properties):
        self.properties = properties
        self.avogadro = scipy.constants.physical_constants["Avogadro constant"]

    def calculate_bulk_density(self):
        # properties = [ molar mass, formula unit, unit cell volume ]
        M = self.properties[0]  # in g/mol
        Z = self.properties[1]
        V = self.properties[2]  # in cm^3

        # density rho in kg/m^3
        rho = (Z*M)/(V*self.avogadro[0])*1000

        return rho

    def calculate_electron_density(self):
        # properties = [ elements, amounts, bulk density ]
        Z = [self.properties[1][i]*self.properties[0][i][1]
             for i in range(len(self.properties[0]))][0]/np.sum(self.properties[1])
        A = [self.properties[1][i]*self.properties[0][i][2]
             for i in range(len(self.properties[0]))][0]/np.sum(self.properties[1])
        rho_b = self.properties[2]

        rho_e = 2*Z/A * rho_b

        return rho_e

    def calculate_volume(self):
        # properties = [ list of lattice lengths, list of lattice angles, crystal system ]
        lenghts = self.properties[0]    # in angstrom
        angles = self.properties[1]     # in degree
        crystalsystem = self.properties[2]

        if crystalsystem == "cubic":
            a = lenghts[0]*10**(-8)
            V = a**3
            return V
        elif crystalsystem == "tetragonal":
            a = lenghts[0]*10**(-8)
            c = lenghts[1]*10**(-8)
            V = a**2 * c
            return V
        elif crystalsystem == "hexagonal":
            a = lenghts[0]*10**(-8)
            c = lenghts[1]*10**(-8)
            angle = 60
            V = (a**2 * c)*np.sin(angle*np.pi/180)
            return V
        elif crystalsystem == "trigonal":
            a = lenghts[0]*10**(-8)
            c = lenghts[1]*10**(-8)
            angle = 60
            V = (a**2 * c)*np.sin(angle*np.pi/180)
            return V
        elif crystalsystem == "orthorhombic":
            a = lenghts[0]*10**(-8)
            b = lenghts[1]*10**(-8)
            c = lenghts[2]*10**(-8)
            V = a * b * c
            return V
        elif crystalsystem == "monoclinic":
            a = lenghts[0]*10**(-8)
            b = lenghts[1]*10**(-8)
            c = lenghts[2]*10**(-8)
            beta = angles[0]
            V = (a * b * c)*np.sin(beta*np.pi/180)
            return V
        elif crystalsystem == "triclinic":
            a = lenghts[0]*10**(-8)
            b = lenghts[1]*10**(-8)
            c = lenghts[2]*10**(-8)
            alpha = angles[0]
            beta = angles[1]
            gamma = angles[2]
            V = (a * b * c)*(1 - np.cos(alpha*np.pi/180)**2 - np.cos(beta*np.pi/180)**2 - np.cos(gamma*np.pi/180)**2 +
                             2*(abs(np.cos(alpha*np.pi/180)*np.cos(beta*np.pi/180)*np.cos(gamma*np.pi/180))))**(0.5)
            return V