#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		petrophysics.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		24.10.2023

#-----------------------------------------------

## MODULES
import numpy as np
import random as rd

## METHODS
class SolidProperties:
    #
    def __init__(self, value_amount, value_property):
        self.value_amount = value_amount
        self.value_property = value_property
    #
    def calculate_solid_density(self):
        pass

class SeismicVelocities:
    #
    def __init__(self, rho_solid, rho_fluid):
        self.rho_solid = rho_solid
        self.rho_fluid = rho_fluid
    #
    def calculate_seismic_velocities(self, rho_limits, vP_limits, vS_limits, delta, porosity):
        ## Density
        rho_min = rho_limits[0]
        rho_max = rho_limits[1]
        ## Seismic Velocities
        center = 1.00
        upper_limit = center + delta
        lower_limit = center - delta
        #
        # P-Wave
        vP_min_theo = vP_limits[0]
        vP_max_theo = vP_limits[1]
        vP_min_sim_upper = int(rd.uniform(center, upper_limit)*vP_min_theo)
        vP_mib_sim_lower = int(rd.uniform(lower_limit, center)*vP_min_theo)
        vP_max_sim_upper = int(rd.uniform(center, upper_limit)*vP_max_theo)
        vP_max_sim_lower = int(rd.uniform(lower_limit, center)*vP_max_theo)
        vP_min = rd.randint(vP_mib_sim_lower, vP_min_sim_upper)
        vP_max = rd.randint(vP_max_sim_lower, vP_max_sim_upper)
        #
        # S-Wave
        vS_min_theo = vS_limits[0]
        vS_max_theo = vS_limits[1]
        vS_min_sim_upper = int(rd.uniform(center, upper_limit)*vS_min_theo)
        vS_mib_sim_lower = int(rd.uniform(lower_limit, center)*vS_min_theo)
        vS_max_sim_upper = int(rd.uniform(center, upper_limit)*vS_max_theo)
        vS_max_sim_lower = int(rd.uniform(lower_limit, center)*vS_max_theo)
        vS_min = rd.randint(vS_mib_sim_lower, vS_min_sim_upper)
        vS_max = rd.randint(vS_max_sim_lower, vS_max_sim_upper)
        #
        constant_a_P = (vP_max - vP_min)/(rho_max - rho_min)
        constant_b_P = vP_max - constant_a_P*rho_max
        constant_a_S = (vS_max - vS_min)/(rho_max - rho_min)
        constant_b_S = vS_max - constant_a_S*rho_max
        #
        condition_v = False
        while condition_v == False:
            ## Porosity
            var_porosity = round(rd.uniform(porosity[0], porosity[1]), 4)
            #
            ## Density
            rho = round((1 - var_porosity)*self.rho_solid + var_porosity*self.rho_fluid, 3)
            #
            ## Seismic Velocities
            vP = constant_a_P * rho + constant_b_P
            #
            if vP_min <= vP <= vP_max:
                vS = constant_a_S * rho + constant_b_S
                if vS_min <= vS <= vS_max:
                    condition_v = True
            #     else:
            #         print("rho:", rho, "vS:", vS)
            # else:
            #    print("rho:", rho, "vP:", vP, vP_max, vP_min)
        #
        vPvS = round(vP/vS, 4)
        #
        return vP, vS, vPvS, rho, var_porosity
    #
    def calculate_elastic_properties(self, rho, vP, vS):
        ## Elastic Parameters
        bulk_modulus = round(rho*(vP**2 - 4/3*vS**2)*10**(-9), 3)
        shear_modulus = round((rho*vS**2)*10**(-9), 3)
        youngs_modulus = round((9*bulk_modulus*shear_modulus)/(3*bulk_modulus + shear_modulus), 3)
        poisson_ratio = round((3*bulk_modulus - 2*shear_modulus)/(6*bulk_modulus + 2*shear_modulus), 4)
        #
        return bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio