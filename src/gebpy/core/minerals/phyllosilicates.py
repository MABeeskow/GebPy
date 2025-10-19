#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		phyllosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		19.10.2025

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
            "Al": PeriodicSystem(name="Al").get_data(),
            "Si": PeriodicSystem(name="Si").get_data(),
            "K": PeriodicSystem(name="K").get_data(),
            "Fe": PeriodicSystem(name="Fe").get_data(),
        }

        # Mineral-specific data
        if self.name == "Annite":
            self.yaml_data = self._load_yaml("annite")

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
        "Annite": self.create_annite,
        "Biotite": self.create_biotite,
        "Eastonite": self.create_eastonite,
        "Phlogopite": self.create_phlogopite,
        "Siderophyllite": self.create_siderophyllite
        }

        if self.name not in generators:
            raise ValueError(f"Mineral '{self.name}' not recognized.")

        dataset = {}
        for index in range(number):
            self.current_seed = hash((self.random_seed, index))%(2**32)
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

    def create_annite(self) -> None:
        """
        Synthetic mineral data generation for annite.
        All mechanical properties (K, G, E) are stored in Pascals internally.
        For output, they are converted to GPa.
        """
        # Random number generator
        rng = np.random.default_rng(self.current_seed)

        # Reading and assigning the mineral-specific information from the YAML file
        val_key = self._get_value(self.yaml_data, ["metadata", "key"])
        val_system = self._get_value(self.yaml_data, ["metadata", "crystal_system"])
        val_K = self._get_value(self.yaml_data, ["physical_properties", "K"])
        val_G = self._get_value(self.yaml_data, ["physical_properties", "G"])
        val_a = self._get_value(self.yaml_data, ["cell_data", "a"])
        val_b = self._get_value(self.yaml_data, ["cell_data", "b"])
        val_c = self._get_value(self.yaml_data, ["cell_data", "c"])
        val_beta = self._get_value(self.yaml_data, ["cell_data", "beta"])
        val_Z = self._get_value(self.yaml_data, ["cell_data", "Z"])

        # Chemistry
        val_state = "fixed"
        traces_data = []
        # Major elements
        hydrogen = self.elements["H"]
        oxygen = self.elements["O"]
        aluminium = self.elements["Al"]
        silicon = self.elements["Si"]
        potassium = self.elements["K"]
        iron = self.elements["Fe"]
        majors_name = ["H", "O", "Al", "Si", "K", "Fe"]

        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 10+2, oxygen[2]],
                                ["Al", aluminium[1], 1, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["K", potassium[1], 1, potassium[2]],
                                ["Fe", iron[1], 3, iron[2]]], dtype=object)

        # Molar mass, elemental amounts
        molar_mass_pure = potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(oxygen[2]+hydrogen[2])

        if not hasattr(self, "cache"):
            self.cache = {}
        if "annite" not in self.cache:
            self.cache["annite"] = {
                "majors_data": majors_data,
                "molar_mass_pure": molar_mass_pure
            }
        else:
            majors_data = self.cache["annite"]["majors_data"]

        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [self.elements[name] for name, *_ in amounts]
        # Density, volume
        dataV = CrystalPhysics([[val_a, val_b, val_c], [val_beta], val_system])
        V = dataV.calculate_volume()
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=val_Z)
        dataRho = CrystalPhysics([molar_mass, val_Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Young's modulus
        E = (9*val_K*val_G)/(3*val_K + val_G)
        # Poisson's ratio
        nu = (3*val_K - 2*val_G)/(2*(3*val_K + val_G))
        # vP/vS
        vPvS = ((val_K + 4/3*val_G)/val_G)**0.5
        # P-wave velocity
        vP = ((val_K + 4/3*val_G)/rho)**0.5
        # S-wave velocity
        vS = (val_G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None

        results = {}
        results["mineral"] = val_key
        results["state"] = val_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = rho
        results["rho_e"] = rho_e
        results["V"] = V_m
        results["vP"] = vP
        results["vS"] = vS
        results["vP/vS"] = vPvS
        results["G"] = val_G*10**(-9)
        results["K"] = val_K*10**(-9)
        results["E"] = E*10**(-9)
        results["nu"] = nu
        results["GR"] = gamma_ray
        results["PE"] = pe
        results["U"] = U
        if p != None:
            results["p"] = p
        else:
            results["p"] = p

        return results

    def create_biotite(self) -> None:
        pass

    def create_eastonite(self) -> None:
        pass

    def create_phlogopite(self) -> None:
        pass

    def create_siderophyllite(self) -> None:
        pass

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = Phyllosilicates(name="Annite", random_seed=42).generate_dataset(number=10)