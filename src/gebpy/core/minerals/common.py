#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		common.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		29.11.2025

#-----------------------------------------------

"""
Module: common.py
This module contains several routines that are commonly used by the different mineral-related modules.
"""

# PACKAGES
import numpy as np
import scipy, re

# MODULES
from modules.geochemistry import MineralChemistry
from modules.geophysics import WellLog as wg

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
        self.avogadro = scipy.constants.Avogadro

    def calculate_bulk_density(self):
        # properties = [ molar mass, formula unit, unit cell volume ]
        M = self.properties[0]  # in g/mol
        Z = self.properties[1]
        V = self.properties[2]  # in cm^3

        # density rho in kg/m^3
        rho = (Z*M)/(V*self.avogadro)*1000

        return rho

    def calculate_electron_density(self):
        # properties = [ elements, amounts, bulk density ]
        Z = np.sum([self.properties[1][i]*self.properties[0][i][1] for i in range(len(self.properties[0]))])/np.sum(
            self.properties[1])
        A = np.sum([self.properties[1][i]*self.properties[0][i][2] for i in range(len(self.properties[0]))])/np.sum(
            self.properties[1])

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
        elif crystalsystem in ["hexagonal", "trigonal"]:
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

class MineralGeneration:
    """
    This class controls the mineral data generation for minerals with fixed composition.
    """
    def __init__(self, name, yaml_data, elements, cache, geophysical_properties, rounding):
        self.name = name
        self.yaml_data = yaml_data
        self.elements = elements
        self.cache = cache
        self.geophysical_properties = geophysical_properties
        self.rounding = rounding

    def _get_value(self, data: dict, path: list[str], default=None):
        """Safely extract a float or string value from nested YAML data."""
        try:
            for key in path:
                data = data[key]
            if isinstance(data, dict) and "value" in data:
                return float(data["value"])
            else:
                return data  # kann z.B. ein String oder eine Zahl sein
        except (KeyError, TypeError):
            return default

    def _determine_volume_constructor(self, vals):
        val_a = vals["a"]
        val_system = vals["crystal_system"]
        if val_system in ["isometric", "cubic"]:
            constr_vol = CrystalPhysics([[val_a], [], val_system])
        elif val_system in ["tetragonal", "hexagonal", "trigonal"]:
            val_c = vals["c"]
            constr_vol = CrystalPhysics([[val_a, val_c], [], val_system])
        elif val_system in ["orthorhombic"]:
            val_b = vals["b"]
            val_c = vals["c"]
            constr_vol = CrystalPhysics([[val_a, val_b, val_c], [], val_system])
        elif val_system in ["monoclinic"]:
            val_b = vals["b"]
            val_c = vals["c"]
            val_beta = vals["beta"]
            constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_beta], val_system])
        elif val_system in ["triclinic"]:
            val_b = vals["b"]
            val_c = vals["c"]
            val_alpha = vals["alpha"]
            val_beta = vals["beta"]
            val_gamma = vals["gamma"]
            constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_alpha, val_beta, val_gamma], val_system])
        return constr_vol

    def _parse_formula(self, formula: str):
        pattern = r"([A-Z][a-z]?)(\d*)"
        matches = re.findall(pattern, formula)

        composition = {}
        for elem, amount in matches:
            amount = int(amount) if amount else 1
            composition[elem] = composition.get(elem, 0) + amount

        return composition

    def _get_cation_element(self, oxide: str) -> str:
        first = oxide[0]
        if len(oxide) > 1 and oxide[1].islower():
            return oxide[:2]

        return first

    def _determine_oxide_conversion_factors(self):
        list_oxides = [
            "H2O", "CO", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "Cl2O", "K2O", "CaO", "MnO", "Mn2O3", "MnO2", "MnO3",
            "Mn2O7", "FeO", "Fe2O3", "FeO3", "NiO", "Ni2O3"]
        mass_oxygen = self.elements["O"][2]
        _conversion_factors = {}
        for oxide in list_oxides:
            _conversion_factors[oxide] = self._parse_formula(formula=oxide)
            cation = self._get_cation_element(oxide=oxide)
            mass_cation = self.elements[cation][2]
            _conversion_factors[oxide]["factor"] = (_conversion_factors[oxide][cation]*mass_cation +
                                                    _conversion_factors[oxide]["O"]*mass_oxygen)/(
                    _conversion_factors[oxide][cation]*mass_cation)

        return _conversion_factors

    def create_mineral_data_fixed_composition(self):
        """
        Synthetic mineral data generation for an user-selected mineral.
        All mechanical properties (K, G, E) are stored in Pascals internally.
        For output, they are converted to GPa.
        """
        name_lower = self.name.lower()
        # Chemistry
        val_state = "fixed"
        traces_data = []
        # Molar mass, elemental amounts
        majors_data = []
        molar_mass_pure = 0
        for element, amount in self.yaml_data["chemistry"].items():
            n_order = int(self.elements[element][1])
            val_amount = float(amount)
            molar_mass = float(self.elements[element][2])
            majors_data.append([element, n_order, val_amount, molar_mass])
            molar_mass_pure += val_amount*molar_mass
        majors_data.sort(key=lambda x: x[1])

        if name_lower not in self.cache:
            vals = {}
            for key in ["K", "G", "a", "b", "c", "alpha", "beta", "gamma", "Z"]:
                if key in self.yaml_data["cell_data"] or key in self.yaml_data["physical_properties"]:
                    vals[key] = self._get_value(self.yaml_data, ["physical_properties", key]) \
                                if key in ["K", "G"] else \
                                self._get_value(self.yaml_data, ["cell_data", key])
            for key in ["key", "crystal_system"]:
                vals[key] = self._get_value(self.yaml_data, ["metadata", key])

            constr_minchem = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data)

            self.cache[name_lower] = {
                "majors_data": majors_data,
                "molar_mass_pure": molar_mass_pure,
                "constants": vals,
                "MineralChemistry": constr_minchem
            }
        else:
            vals = self.cache[name_lower]["constants"]
            constr_minchem = self.cache[name_lower]["MineralChemistry"]
            constr_electr_density = self.cache[name_lower]["const_electron_density"]
            constr_vol = self.cache[name_lower]["constr_volume"]
            constr_density = self.cache[name_lower]["constr_density"]
            constr_radiation = self.cache[name_lower]["constr_radiation"]

        # Reading and assigning the mineral-specific information from the YAML file
        val_key = vals["key"]
        val_K = vals["K"]
        val_G = vals["G"]
        val_Z = vals["Z"]

        molar_mass, amounts = constr_minchem.calculate_molar_mass()
        element = [self.elements[name] for name, *_ in amounts]
        # (Molar) Volume
        if "constr_volume" not in self.cache[name_lower]:
            constr_vol = self._determine_volume_constructor(vals=vals)
            self.cache[name_lower]["constr_volume"] = constr_vol

        V, V_m = CrystallographicProperties().calculate_molar_volume(
            constr_volume=constr_vol, constr_molar_volume=constr_minchem, cell_z=val_Z)
        # Density
        if "constr_density" not in self.cache[name_lower]:
            constr_density = CrystalPhysics([molar_mass, val_Z, V])
            self.cache[name_lower]["constr_density"] = constr_density

        rho = CrystallographicProperties().calculate_mineral_density(constr_density=constr_density)

        if "const_electron_density" not in self.cache[name_lower]:
            constr_electr_density = wg(amounts=amounts, elements=element, rho_b=rho)
            self.cache[name_lower]["const_electron_density"] = constr_electr_density

        rho_e = CrystallographicProperties().calculate_electron_density(constr_electron_density=constr_electr_density)
        # Elastic properties
        E, nu = self.geophysical_properties.calculate_elastic_properties(bulk_mod=val_K, shear_mod=val_G)
        # Seismic properties
        vPvS, vP, vS = self.geophysical_properties.calculate_seismic_velocities(
            bulk_mod=val_K, shear_mod=val_G, rho=rho)
        # Radiation properties
        if "constr_radiation" not in self.cache[name_lower]:
            constr_radiation = wg(amounts=amounts, elements=element)
            self.cache[name_lower]["constr_radiation"] = constr_radiation

        gamma_ray, pe, U = self.geophysical_properties.calculate_radiation_properties(
            constr_radiation=constr_radiation, rho_electron=rho_e)
        # Electrical resistivity
        p = None
        # Results
        results = {
            "mineral": val_key, "state": val_state, "M": round(molar_mass, self.rounding),
            "chemistry": {name: round(val[1], 6) for name, *val in amounts}, "rho": round(rho, self.rounding),
            "rho_e": round(rho_e, self.rounding), "V": round(V_m, self.rounding), "vP": round(vP, self.rounding),
            "vS": round(vS, self.rounding), "vP/vS": round(vPvS, self.rounding),
            "K": round(val_K*10**(-9), self.rounding), "G": round(val_G*10**(-9), self.rounding),
            "E": round(E*10**(-9), self.rounding), "nu": round(nu, 6), "GR": round(gamma_ray, self.rounding),
            "PE": round(pe, self.rounding), "U": round(U, self.rounding), "p": p}
        return results