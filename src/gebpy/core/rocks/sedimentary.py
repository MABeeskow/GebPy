#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sedimentary.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		13.12.2025

#-----------------------------------------------

"""
Module: sedimentary.py
This module controls the generation of the synthetic data of the oxides minerals.
"""

# PACKAGES
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# MODULES
from src.gebpy.core.minerals.synthesis import MineralDataGeneration
from src.gebpy.core.rocks.common import RockGeneration

# Code
BASE_PATH = Path(__file__).resolve().parents[2]
DATA_PATH = BASE_PATH / "data_rocks"

class SedimentaryRocks:
    _yaml_cache = {}
    _mineralogy_cache = {}
    _rocks = {"Sandstone", "Limestone", "Shale"}

    def __init__(self, name, random_seed, rounding: int = 3) -> None:
        self.name = name
        self.random_seed = random_seed
        self.rng = np.random.default_rng(random_seed)
        self.current_seed = int(np.round(self.rng.uniform(0, 1000), 0))
        self.rounding = rounding
        self.data_path = DATA_PATH
        self.conversion_factors = RockGeneration()._determine_oxide_conversion_factors()

    def _load_yaml(self, rock_name: str) -> dict:
        # 1) Cache-Hit
        if rock_name in SedimentaryRocks._yaml_cache:
            return SedimentaryRocks._yaml_cache[rock_name]

        # 2) Laden von Disk
        yaml_file = self.data_path/f"{rock_name}.yaml"
        if not yaml_file.exists():
            raise FileNotFoundError(f"No YAML file found for {rock_name}.")

        with open(yaml_file, "r") as f:
            data = yaml.safe_load(f)

        if "mineralogy" in data and rock_name not in SedimentaryRocks._mineralogy_cache:
            self._compile_mineralogy(rock_name, data["mineralogy"])

        # 3) Cache schreiben
        SedimentaryRocks._yaml_cache[rock_name] = data

        return data

    def _compile_mineralogy(self, rock_name: str, mineralogy_dict: dict):
        """
        Extracts and compiles all chemistry formulas from the YAML file.
        Stores the compiled ASTs in the global formula cache.
        """
        if rock_name not in SedimentaryRocks._mineralogy_cache:
            SedimentaryRocks._mineralogy_cache[rock_name] = {}

        for element, entry in mineralogy_dict.items():
            mineral = element
            interval = list(entry.values())
            lower_limit = interval[0]
            upper_limit = interval[1]
            compiled = [lower_limit, upper_limit]
            SedimentaryRocks._mineralogy_cache[rock_name][mineral] = compiled

    def _sample_mineralogy(self, rock_name: str):
        mineral_ranges = SedimentaryRocks._mineralogy_cache[rock_name]
        sampled = {}

        for mineral, (lower, upper) in mineral_ranges.items():
            sampled[mineral] = np.random.uniform(lower, upper)

        # Normalisieren auf Summe = 1
        total = sum(sampled.values())
        sampled = {m: v/total for m, v in sampled.items()}

        return sampled

    def _sample_bounded_simplex(self, min_vals, max_vals):
        min_vals = np.array(min_vals)
        max_vals = np.array(max_vals)
        n = len(min_vals)

        # Schritt 1: freie Menge
        R = 1 - min_vals.sum()
        if R < 0:
            raise ValueError("Sum of minimum fractions > 1, impossible mixture")

        # Dirichlet-Sample
        y = self.rng.dirichlet(alpha=np.ones(n))

        # Skaliere y auf die maximal zulässige Stretch-Länge
        stretch = max_vals - min_vals
        y = y * R

        # enforce max boundary by scaling down
        overshoot = y > stretch
        if overshoot.any():
            # Skaliere nur die Komponenten runter, die zu groß wären
            scale = np.min(stretch[overshoot] / y[overshoot])
            y *= scale

        # final mix
        x = y + min_vals
        x /= x.sum()  # kleine numerische Korrektur

        return x

    def _collect_mineral_data(self, list_minerals, number):
        _mineral_data = []
        for index, mineral in enumerate(list_minerals):
            data_init = MineralDataGeneration(mineral, number)
            data_mineral = data_init.generate_data()
            is_fixed = data_mineral.shape[0] == 1
            if is_fixed and number > 1:
                data_mineral = pd.concat([data_mineral]*number, ignore_index=True)
            elif not is_fixed and data_mineral.shape[0] != number:
                raise ValueError(
                    f"Mineral '{mineral}' returned {data_mineral.shape[0]} rows, "
                    f"but expected {number}.")
            _mineral_data.append(data_mineral)

        return _mineral_data

    def _extract_mineral_property_data(self, list_minerals, data_mineral, property):
        _helper = []
        for index, mineral in enumerate(list_minerals):
            dataset = data_mineral[index][property]
            _helper.append(list(dataset))

        return np.array(_helper)

    def _update_bulk_density_data(self, _bulk_data, density_fluid, number):
        porosity = _bulk_data["porosity"]
        _bulk_data["rho_s"] = _bulk_data["rho"]
        _bulk_data["rho"] = (1 - porosity)*_bulk_data["rho"] + porosity*density_fluid
        _bulk_data["rho_f"] = np.ones(number)*density_fluid

        return _bulk_data

    def _update_seismic_velocities(self, _bulk_data):
        val_G = _bulk_data["G"]*10**9
        val_rho = _bulk_data["rho"]
        val_K = _bulk_data["K"]*10**9
        _bulk_data["vS"] = (val_G/val_rho)**0.5
        _bulk_data["vP"] = ((val_K + 4/3*val_G)/val_rho)**0.5
        _bulk_data["vP/vS"] = _bulk_data["vP"]/_bulk_data["vS"]

        return _bulk_data

    def _update_elastic_parameter_data(self, _bulk_data):
        val_G = _bulk_data["G"]*10**9
        val_K = _bulk_data["K"]*10**9
        _bulk_data["E"] = ((9*val_K*val_G)/(3*val_K + val_G))*1e-9
        _bulk_data["poisson"] = (3*val_K - 2*val_G)/(2*(3*val_K + val_G))
        _bulk_data["lame"] = (val_K - (2*val_G)/3)*1e-9

        return _bulk_data

    def _extract_element_data(self, data_minerals, list_elements):
        for index, dataset in enumerate(data_minerals):
            list_keys_i = list(dataset.keys())
            filtered = [x[len("chemistry."):] for x in list_keys_i if x.startswith("chemistry.")]
            for element in filtered:
                if element not in list_elements:
                    list_elements.append(element)

        return list_elements

    def _extract_oxide_data(self, data_minerals, list_oxides):
        for index, dataset in enumerate(data_minerals):
            list_keys_i = list(dataset.keys())
            filtered = [x[len("compounds."):] for x in list_keys_i if x.startswith("compounds.")]
            for element in filtered:
                if element not in list_oxides:
                    list_oxides.append(element)
            if dataset["mineral"][0] in ["Py"]:
                list_oxides.append("Fe2O3")
                list_oxides.append("SO3")
                list_oxides.remove("FeS2")

        return list_oxides

    def _update_chemistry_data(self, _bulk_data, data_minerals, data_composition, element, number):
        _helper = [[] for _ in range(len(data_minerals))]
        key_element = "chemistry." + element
        for index, dataset in enumerate(data_minerals):
            if key_element in list(dataset.keys()):
                data = list(dataset[key_element])
                _helper[index] = data
            else:
                data = [0]*number
                _helper[index] = data
        _helper = np.array(_helper)
        _helper_bulk = np.sum(data_composition*_helper.T, axis=1)
        key_element = "w." + element
        _bulk_data[key_element] = np.array(_helper_bulk)

        return _bulk_data

    def _update_oxide_data(self, bulk_data, list_oxides):
        for oxide in list_oxides:
            cation = RockGeneration()._get_cation_element(oxide=oxide)
            anion = RockGeneration()._get_anion_element(compound=oxide)
            if anion == "O":
                key_cation = "w." + cation
                values = self.conversion_factors[oxide]["factor"]*bulk_data[key_cation]
                key_oxide = "w." + oxide
                bulk_data[key_oxide] = values
            else:
                print("There is a non-oxide compound part of the list.")

        return bulk_data

    def _assign_mineral_amounts(self, bulk_data, data_amounts):
        for mineral, values in data_amounts.items():
            key_mineral = "phi." + mineral
            bulk_data[key_mineral] = values

        return bulk_data

    def generate_dataset(self, number: int = 1, fluid: str = "water", density_fluid=None) -> None:
        if density_fluid is None:
            if fluid == "water":
                density_fluid = 1000
            elif fluid == "oil":
                density_fluid = 800
            elif fluid == "natural gas":
                density_fluid = 750

        siliciclastics = {"Sandstone"}
        data_yaml = self._load_yaml(rock_name=self.name)
        min_porosity = data_yaml["physical_properties"]["porosity"]["min"]
        max_porosity = data_yaml["physical_properties"]["porosity"]["max"]
        porosity = self.rng.uniform(min_porosity, max_porosity, number)
        mineral_limits = SedimentaryRocks._mineralogy_cache[self.name]
        _limits = {"lower": [], "upper": []}

        for mineral, values in mineral_limits.items():
            _limits["lower"].append(values[0])
            _limits["upper"].append(values[1])

        list_minerals = list(SedimentaryRocks._mineralogy_cache[self.name].keys())
        _properties = ["rho", "vP", "vS", "K", "G", "GR", "PE"]
        _bulk_data = {}
        # Collect mineralogical composition data
        _helper_amounts = {}
        _helper_composition = []
        _helper_mineral_amounts = {}
        for index_main in range(number):
            mineral_amounts = self._sample_bounded_simplex(min_vals=_limits["lower"], max_vals=_limits["upper"])
            _helper_composition.append(list(mineral_amounts))
            _helper_dict = {}
            for j, mineral in enumerate(list_minerals):
                _helper_dict[mineral] = mineral_amounts[j]
                if mineral not in _helper_mineral_amounts:
                    _helper_mineral_amounts[mineral] = []
                _helper_mineral_amounts[mineral].append(mineral_amounts[j])
            _helper_amounts[index_main] = _helper_dict
        _helper_composition = np.array(_helper_composition)
        # Collect mineral data
        _helper_elements = []
        _helper_oxides = []
        _mineral_data = self._collect_mineral_data(list_minerals=list_minerals, number=number)
        _helper_elements = self._extract_element_data(data_minerals=_mineral_data, list_elements=_helper_elements)
        _helper_oxides = self._extract_oxide_data(data_minerals=_mineral_data, list_oxides=_helper_oxides)
        _helper_bulk_data = {}
        for property in ["rho", "K", "G", "GR", "PE"]:
            _helper_property = self._extract_mineral_property_data(
                list_minerals=list_minerals, data_mineral=_mineral_data, property=property)
            _helper_bulk_data[property] = np.sum(_helper_composition*_helper_property.T, axis=1)
        # Assign porosity data
        _helper_bulk_data["porosity"] = porosity
        # Update bulk density data
        _helper_bulk_data = self._update_bulk_density_data(
            _bulk_data=_helper_bulk_data, density_fluid=density_fluid, number=number)
        # Update bulk seismic velocity data
        _helper_bulk_data = self._update_seismic_velocities(_bulk_data=_helper_bulk_data)
        # Update elastic parameter data
        _helper_bulk_data = self._update_elastic_parameter_data(_bulk_data=_helper_bulk_data)
        # Update chemistry data
        for element in _helper_elements:
            _helper_bulk_data = self._update_chemistry_data(
                _bulk_data=_helper_bulk_data, data_minerals=_mineral_data, data_composition=_helper_composition,
                element=element, number=number)
        # Update oxide data
        _helper_bulk_data = self._update_oxide_data(bulk_data=_helper_bulk_data, list_oxides=_helper_oxides)
        # Update rock composition data
        _helper_bulk_data = self._assign_mineral_amounts(
            bulk_data=_helper_bulk_data, data_amounts=_helper_mineral_amounts)
        # Conversion to a pandas dataframe object
        _bulk_data = pd.DataFrame(_helper_bulk_data)

        return _bulk_data

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = SedimentaryRocks(name="Sandstone", random_seed=42).generate_dataset(number=10)