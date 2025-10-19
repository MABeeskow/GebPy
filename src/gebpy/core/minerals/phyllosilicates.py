#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		phyllosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		20.10.2025

#-----------------------------------------------

"""
Module: phyllosilicates.py
This module controls the generation of the synthetic data of the phyllosilicate minerals.
"""

# PACKAGES
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# MODULES
from modules.chemistry import PeriodicSystem
from modules.geochemistry import MineralChemistry
from modules.geophysics import WellLog as wg
from modules.minerals import CrystalPhysics
from src.gebpy.core.minerals.common import GeophysicalProperties, CrystallographicProperties

# CODE
class Phyllosilicates:
    def __init__(self, name, random_seed) -> None:
        self.name = name
        self.random_seed = random_seed
        self.rng = np.random.default_rng(random_seed)
        self.current_seed = int(np.round(self.rng.uniform(0, 1000), 0))
        self.data_path = Path(__file__).resolve().parents[2] / "data"

        # Chemistry
        self.elements = {
            "H": PeriodicSystem(name="H").get_data(),
            "O": PeriodicSystem(name="O").get_data(),
            "Mg": PeriodicSystem(name="Mg").get_data(),
            "Al": PeriodicSystem(name="Al").get_data(),
            "Si": PeriodicSystem(name="Si").get_data(),
            "K": PeriodicSystem(name="K").get_data(),
            "Fe": PeriodicSystem(name="Fe").get_data(),
        }

        # Geophysics
        self.geophysical_properties = GeophysicalProperties()

        # Mineral-specific data
        if self.name == "Annite":
            self.yaml_data = self._load_yaml("annite")
        elif self.name == "Phlogopite":
            self.yaml_data = self._load_yaml("phlogopite")
        elif self.name == "Siderophyllite":
            self.yaml_data = self._load_yaml("siderophyllite")
        elif self.name == "Eastonite":
            self.yaml_data = self._load_yaml("eastonite")

    def _load_yaml(self, mineral_name: str) -> dict:
        yaml_file = self.data_path / f"{mineral_name.lower()}.yaml"
        if not yaml_file.exists():
            raise FileNotFoundError(f"No YAML file found for {mineral_name}.")
        with open(yaml_file, "r") as f:
            return yaml.safe_load(f)

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

    def generate_dataset(self, number: int = 1) -> None:
        generators = {
        "Annite": self.create_mineral_data_fixed_composition,
        "Biotite": self.create_biotite,
        "Eastonite": self.create_mineral_data_fixed_composition,
        "Phlogopite": self.create_mineral_data_fixed_composition,
        "Siderophyllite": self.create_mineral_data_fixed_composition,
        }

        if self.name not in generators:
            raise ValueError(f"Mineral '{self.name}' not recognized.")

        dataset = {}
        for index in range(number):
            self.current_seed = np.uint32(self.random_seed + index)
            data_mineral = generators[self.name]()

            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U", "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif key in ["mineral", "state", "trace elements"] and key not in dataset:
                    dataset[key] = value
                elif key in ["chemistry", "compounds"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)

        return dataset

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
            val_amount = int(amount)
            molar_mass = float(self.elements[element][2])
            majors_data.append([element, n_order, val_amount, molar_mass])
            molar_mass_pure += val_amount*molar_mass
        majors_data.sort(key=lambda x: x[1])

        if not hasattr(self, "cache"):
            self.cache = {}

        if name_lower not in self.cache:
            vals = {}
            for key in ["K", "G", "a", "b", "c", "beta", "Z"]:
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
        val_system = vals["crystal_system"]
        val_K = vals["K"]
        val_G = vals["G"]
        val_a = vals["a"]
        val_b = vals["b"]
        val_c = vals["c"]
        val_beta = vals["beta"]
        val_Z = vals["Z"]

        molar_mass, amounts = constr_minchem.calculate_molar_mass()
        element = [self.elements[name] for name, *_ in amounts]
        # (Molar) Volume
        if "constr_volume" not in self.cache[name_lower]:
            constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_beta], val_system])
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
            "mineral": val_key, "state": val_state, "M": molar_mass,
            "chemistry": {name: val[1] for name, *val in amounts}, "rho": rho, "rho_e": rho_e, "V": V_m, "vP": vP,
            "vS": vS, "vP/vS": vPvS, "K": val_K*10**(-9), "G": val_G*10**(-9), "E": E*10**(-9), "nu": nu,
            "GR": gamma_ray, "PE": pe, "U": U, "p": p}
        return results

    def create_biotite(self) -> None:
        pass

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = Phyllosilicates(name="Annite", random_seed=42).generate_dataset(number=10)