#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		common.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: common.py
This module contains several routines that are commonly used by the different mineral-related modules.
"""

# PACKAGES
import numpy as np
import scipy, re

# MODULES
from ..chemistry.geochemistry import MineralChemistry
from ..physics.geophysics import WellLog as wg

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

class OxideComposition:
    """
    This class calculates the oxide composition of a certain mineral.
    """
    def __init__(self):
        pass

    def _get_cation_element(self, oxide: str) -> str:
        first = oxide[0]
        if len(oxide) > 1 and oxide[1].islower():
            return oxide[:2]

        return first

    def _parse_formula(self, formula: str):
        pattern = r"([A-Z][a-z]?)(\d*)"
        matches = re.findall(pattern, formula)

        composition = {}
        for elem, amount in matches:
            amount = int(amount) if amount else 1
            composition[elem] = composition.get(elem, 0) + amount

        return composition

    def _determine_oxide_conversion_factors(self, elements):
        list_oxides = [
            "H2O", "CO", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "Cl2O", "K2O", "CaO", "MnO", "Mn2O3", "MnO2", "MnO3",
            "Mn2O7", "FeO", "Fe2O3", "FeO3", "NiO", "Ni2O3", "TiO2", "Ti2O3", "VO", "V2O3", "VO2", "V2O10", "CrO",
            "Cr2O3", "CrO3", "CoO", "Co2O3", "Cu2O", "CuO", "ZnO", "GeO2", "As2O3", "As2O10", "ZrO2", "Nb2O3", "Nb2O10",
            "MoO", "Mo2O3", "MoO2", "Mo2O10", "MoO3", "CdO", "SnO", "SnO2", "Sb2O3", "Sb2O10", "TeO2", "TeO3", "Ta2O10",
            "WO", "W2O3", "WO2", "W2O10", "WO3", "Au2O", "Au2O3", "Tl2O", "Tl2O3", "PbO", "PbO2", "Bi2O3", "Bi2O10",
            "U2O3", "UO2", "U2O10", "UO3"]
        mass_oxygen = elements["O"][2]
        _conversion_factors = {}
        for oxide in list_oxides:
            _conversion_factors[oxide] = self._parse_formula(formula=oxide)
            cation = self._get_cation_element(oxide=oxide)
            if cation in elements:
                mass_cation = elements[cation][2]
                _conversion_factors[oxide]["factor"] = (_conversion_factors[oxide][cation]*mass_cation +
                                                        _conversion_factors[oxide]["O"]*mass_oxygen)/(
                        _conversion_factors[oxide][cation]*mass_cation)
            else:
                pass
                #print(self.name, ": cation", cation, "not found in chemical container.")

        return _conversion_factors

    def _element_amounts_as_dict(self, amounts):
        helper_dict = {}
        for item in amounts:
            helper_dict[item[0]] = item[2]

        return helper_dict

class MineralGeneration:
    """
    This class controls the mineral data generation for minerals with fixed composition.
    """
    def __init__(
            self, name, yaml_data, elements, cache, geophysical_properties, rounding, rng, variability=False,
            uncertainty=1.0):
        self.name = name
        self.yaml_data = yaml_data
        self.elements = elements
        self.cache = cache
        self.geophysical_properties = geophysical_properties
        self.rounding = rounding
        self.rng = rng
        self.variability = variability
        self.uncertainty = uncertainty
        self.conversion_factors = self._determine_oxide_conversion_factors()

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

    def _determine_sulfide_conversion_factors(self):
        list_sulfides = ["FeS", "FeS2"]
        mass_sulfur = self.elements["S"][2]
        _conversion_factors = {}
        for sulfide in list_sulfides:
            _conversion_factors[sulfide] = self._parse_formula(formula=sulfide)
            cation = self._get_cation_element(oxide=sulfide)
            if cation in self.elements:
                mass_cation = self.elements[cation][2]
                _conversion_factors[sulfide]["factor"] = (_conversion_factors[sulfide][cation]*mass_cation +
                                                        _conversion_factors[sulfide]["S"]*mass_sulfur)/(
                        _conversion_factors[sulfide][cation]*mass_cation)
            else:
                pass
                #print(self.name, ": cation", cation, "not found in chemical container.")

        return _conversion_factors

    def _determine_oxide_conversion_factors(self):
        list_oxides = [
            "H2O", "CO", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "Cl2O", "K2O", "CaO", "MnO", "Mn2O3", "MnO2", "MnO3",
            "Mn2O7", "FeO", "Fe2O3", "FeO3", "NiO", "Ni2O3", "TiO2", "Ti2O3", "VO", "V2O3", "VO2", "V2O10", "CrO",
            "Cr2O3", "CrO3", "CoO", "Co2O3", "Cu2O", "CuO", "ZnO", "GeO2", "As2O3", "As2O10", "ZrO2", "Nb2O3", "Nb2O10",
            "MoO", "Mo2O3", "MoO2", "Mo2O10", "MoO3", "CdO", "SnO", "SnO2", "Sb2O3", "Sb2O10", "TeO2", "TeO3", "Ta2O10",
            "WO", "W2O3", "WO2", "W2O10", "WO3", "Au2O", "Au2O3", "Tl2O", "Tl2O3", "PbO", "PbO2", "Bi2O3", "Bi2O10",
            "U2O3", "UO2", "U2O10", "UO3", "Nb2O5", "Ta2O5"]
        mass_oxygen = self.elements["O"][2]
        _conversion_factors = {}
        for oxide in list_oxides:
            _conversion_factors[oxide] = self._parse_formula(formula=oxide)
            cation = self._get_cation_element(oxide=oxide)
            if cation in self.elements:
                mass_cation = self.elements[cation][2]
                _conversion_factors[oxide]["factor"] = (_conversion_factors[oxide][cation]*mass_cation +
                                                        _conversion_factors[oxide]["O"]*mass_oxygen)/(
                        _conversion_factors[oxide][cation]*mass_cation)
            else:
                pass
                #print(self.name, ": cation", cation, "not found in chemical container.")

        return _conversion_factors

    def _element_amounts_as_dict(self, amounts):
        helper_dict = {}
        for item in amounts:
            helper_dict[item[0]] = item[2]

        return helper_dict

    def _determine_majors_data(self):
        majors_data = []
        molar_mass_pure = 0
        for element, amount in self.yaml_data["chemistry"].items():
            n_order = int(self.elements[element][1])
            val_amount = float(amount)
            molar_mass = float(self.elements[element][2])
            majors_data.append([element, n_order, val_amount, molar_mass])
            molar_mass_pure += val_amount*molar_mass
        majors_data.sort(key=lambda x: x[1])
        return majors_data, molar_mass_pure

    def _extract_values_from_yaml(self):
        vals = {}
        # Physical parameters
        for key in ["K", "G", "a_K", "b_K", "a_G", "b_G"]:
            if key in self.yaml_data.get("physical_properties", {}):
                vals[key] = float(self.yaml_data["physical_properties"][key]["value"])
        # Cell parameters
        for key in ["a", "b", "c", "alpha", "beta", "gamma", "Z"]:
            if key in self.yaml_data.get("cell_data", {}):
                vals[key] = float(self.yaml_data["cell_data"][key]["value"])
        # Meta data
        vals["key"] = self.yaml_data["metadata"]["key"]
        vals["crystal_system"] = self.yaml_data["metadata"]["crystal_system"]

        return vals

    def vary_elastic_moduli(self, K, G, uncertainty):
        delta_K = K*uncertainty/100
        delta_G = G*uncertainty/100
        lower_K = K - delta_K
        upper_K = K + delta_K
        lower_G = G - delta_G
        upper_G = G + delta_G
        K = self.rng.uniform(lower_K, upper_K)
        G = self.rng.uniform(lower_G, upper_G)

        return K, G

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
        majors_data, molar_mass_pure = self._determine_majors_data()
        if "oxides" in self.yaml_data:
            oxides_data = {}
            for oxide, amount in self.yaml_data["oxides"].items():
                cation = self._get_cation_element(oxide=oxide)
                oxides_data[oxide] = [cation, None, amount]
        elif "sulfides" in self.yaml_data:
            self.conversion_factors = self._determine_sulfide_conversion_factors()
            oxides_data = {}
            for oxide, amount in self.yaml_data["sulfides"].items():
                cation = self._get_cation_element(oxide=oxide)
                oxides_data[oxide] = [cation, None, amount]

        if name_lower not in self.cache:
            vals = self._extract_values_from_yaml()
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
        val_Z = vals["Z"]
        if "K" in vals:
            val_K = vals["K"]
            val_G = vals["G"]
        else:
            val_a_K = float(vals["a_K"])
            val_b_K = float(vals["b_K"])
            val_a_G = float(vals["a_G"])
            val_b_G = float(vals["b_G"])

        molar_mass, amounts = constr_minchem.calculate_molar_mass()
        amounts_dict = self._element_amounts_as_dict(amounts=amounts)
        element = [self.elements[name] for name, *_ in amounts]
        for oxide in oxides_data.keys():
            cation = self._get_cation_element(oxide=oxide)
            weight = oxides_data[oxide][2]
            if weight != 1:
                value = weight
            else:
                try:
                    value = amounts_dict[cation]*self.conversion_factors[oxide]["factor"]
                except:
                    value = 0
            oxides_data[oxide][1] = value
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
        if "K" not in vals:
            val_K = (val_a_K*rho + val_b_K)*10**9
            val_G = (val_a_G*rho + val_b_G)*10**9

        if self.variability == True:
            val_K, val_G = self.vary_elastic_moduli(K=val_K, G=val_G, uncertainty=self.uncertainty)

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
        if "oxides" in self.yaml_data:
            results["compounds"] = {name: round(val[1], 6) for name, val in oxides_data.items()}
        elif "sulfides" in self.yaml_data:
            results["compounds"] = {name: round(val[1], 6) for name, val in oxides_data.items()}

        return results

    def create_mineral_data_endmember_series(
            self, endmember_series, var_class, current_seed):
        """
        Synthetic mineral data generation for an user-selected mineral.
        All mechanical properties (K, G, E) are stored in Pascals internally.
        For output, they are converted to GPa.
        """
        val_state = "variable"

        if not hasattr(self, "cache"):
            self.cache = {}

        name_lower = endmember_series[self.name]["name_lower"]
        val_key = endmember_series[self.name]["key"]
        endmember = endmember_series[self.name]["endmembers"]
        oxides_data = {}
        for oxide in endmember_series[self.name]["oxides"]:
            cation = self._get_cation_element(oxide=oxide)
            oxides_data[oxide] = [cation, None]

        if "endmembers" not in self.cache:
            self.cache["endmembers"] = {}

        if name_lower not in self.cache:
            endmember_data = {}
            list_elements = []
            for mineral in endmember:
                if mineral not in self.cache["endmembers"]:
                    mineral_data = var_class(name=mineral, random_seed=current_seed).generate_dataset(
                        number=1)
                    self.cache["endmembers"][mineral] = mineral_data
                endmember_data[mineral] = self.cache["endmembers"][mineral]
                mineral_data = endmember_data[mineral]
                for element in mineral_data["chemistry"]:
                    if element not in list_elements:
                        list_elements.append(element)
            constr_OxComp = OxideComposition()

            self.cache[name_lower] = {
                "endmember_data": endmember_data, "list_elements": list_elements, "OxComp": constr_OxComp}
        else:
            endmember_data = self.cache[name_lower]["endmember_data"]
            list_elements = self.cache[name_lower]["list_elements"]
            constr_OxComp = self.cache[name_lower]["OxComp"]
        weights = self.rng.dirichlet(np.ones(len(endmember)))
        fraction_endmember = dict(zip(endmember, weights))

        properties = ["M", "rho", "rho_e", "V", "K", "G"]
        helper_results = {
            prop: sum(fraction_endmember[m]*endmember_data[m][prop][0] for m in endmember)
            for prop in properties
        }
        # Amounts
        amounts = []
        for element in list_elements:
            amount = sum(fraction_endmember[mineral]*endmember_data[mineral]["chemistry"].get(element, [0])[0]
                         for mineral in endmember)
            amounts.append([element, self.elements[element][1], amount])
        element = [self.elements[name] for name, *_ in amounts]
        # Oxide amounts
        amounts_dict = constr_OxComp._element_amounts_as_dict(amounts=amounts)
        try:
            for oxide in oxides_data.keys():
                cation = constr_OxComp._get_cation_element(oxide=oxide)
                value = amounts_dict[cation]*self.conversion_factors[oxide]["factor"]
                oxides_data[oxide][1] = value
        except:
            print("No oxide data available!")
        # Elastic properties
        val_K = helper_results["K"]*10**9
        val_G = helper_results["G"]*10**9

        if self.variability == True:
            val_K, val_G = self.vary_elastic_moduli(K=val_K, G=val_G, uncertainty=self.uncertainty)

        rho = helper_results["rho"]
        rho_e = helper_results["rho_e"]
        E, nu = self.geophysical_properties.calculate_elastic_properties(bulk_mod=val_K, shear_mod=val_G)
        # Seismic properties
        vPvS, vP, vS = self.geophysical_properties.calculate_seismic_velocities(
            bulk_mod=val_K, shear_mod=val_G, rho=rho)
        # Radiation properties
        constr_radiation = wg(amounts=amounts, elements=element)
        gamma_ray, pe, U = self.geophysical_properties.calculate_radiation_properties(
            constr_radiation=constr_radiation, rho_electron=rho_e)
        # Electrical resistivity
        p = None
        # Results
        results = {
            "mineral": val_key, "state": val_state, "M": round(helper_results["M"], self.rounding),
            "chemistry": {name: round(val[1], 6) for name, *val in amounts}, "rho": round(rho, self.rounding),
            "rho_e": round(rho_e, self.rounding), "V": round(helper_results["V"], self.rounding),
            "vP": round(vP, self.rounding), "vS": round(vS, self.rounding),
            "vP/vS": round(vPvS, self.rounding), "K": round(val_K*10**(-9), self.rounding),
            "G": round(val_G*10**(-9), self.rounding), "E": round(E*10**(-9), self.rounding), "nu": round(nu, 6),
            "GR": round(gamma_ray, self.rounding), "PE": round(pe, self.rounding), "U": round(U, self.rounding), "p": p}
        try:
            results["compounds"] = {name: round(val[1], 6) for name, val in oxides_data.items()}
        except:
            print("No oxide/sulfide data available!")
        return results