#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		isotropic_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: isotropic_rocks.py
This module controls the generation of synthetic data for isotropic rocks.
"""

# PACKAGES
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# MODULES
from ..minerals.synthesis import MineralDataGeneration
from ..rocks.common import RockGeneration
from ..physics.common import Geophysics

# Code
BASE_PATH = Path(__file__).resolve().parents[2]
DATA_PATH = BASE_PATH / "data_rocks"

class IsotropicRocks:
    _yaml_cache = {}
    _mineralogy_cache = {}
    _mineral_groups_cache = {}
    _rocks = {"Sandstone", "Limestone", "Dolostone", "Marl"}

    def __init__(self, name, random_seed, alpha_K=1.0, alpha_G=1.5, variability=False, uncertainty=1.0) -> None:
        self.name = name
        self.random_seed = random_seed
        self.alpha_K = alpha_K
        self.alpha_G = alpha_G
        self.variability = variability
        self.uncertainty = uncertainty
        self.rng = np.random.default_rng(random_seed)
        self.current_seed = int(np.round(self.rng.uniform(0, 1000), 0))
        self.data_path = DATA_PATH
        self.rock_gen = RockGeneration()
        self.geophysics = Geophysics()
        self.conversion_factors = self.rock_gen._determine_oxide_conversion_factors()
        self.cache = {}

    def _load_yaml(self, rock_name: str) -> dict:
        # 1) Cache-Hit
        if rock_name in IsotropicRocks._yaml_cache:
            return IsotropicRocks._yaml_cache[rock_name]

        # 2) Laden von Disk
        yaml_file = self.data_path/f"{rock_name}.yaml"
        if not yaml_file.exists():
            raise FileNotFoundError(f"No YAML file found for {rock_name}.")

        with open(yaml_file, "r") as f:
            data = yaml.safe_load(f)

        if "mineralogy" in data and rock_name not in IsotropicRocks._mineralogy_cache:
            self._compile_mineralogy(rock_name, data["mineralogy"])

        if "mineral_groups" in data:
            self._compile_mineral_groups(rock_name, data["mineral_groups"])

        # 3) Cache schreiben
        IsotropicRocks._yaml_cache[rock_name] = data

        return data

    def _compile_mineralogy(self, rock_name: str, mineralogy_dict: dict):
        """
        Extracts and compiles all chemistry formulas from the YAML file.
        Stores the compiled ASTs in the global formula cache.
        """
        if rock_name not in IsotropicRocks._mineralogy_cache:
            IsotropicRocks._mineralogy_cache[rock_name] = {}

        for element, entry in mineralogy_dict.items():
            mineral = element
            interval = list(entry.values())
            lower_limit = interval[0]
            upper_limit = interval[1]
            compiled = [lower_limit, upper_limit]
            IsotropicRocks._mineralogy_cache[rock_name][mineral] = compiled

    def _compile_mineral_groups(self, rock_name: str, group_dict: dict):
        if rock_name not in IsotropicRocks._mineral_groups_cache:
            IsotropicRocks._mineral_groups_cache[rock_name] = {}

        for group, entry in group_dict.items():
            minerals = entry["minerals"]
            min_val = entry["min"]
            max_val = entry["max"]

            IsotropicRocks._mineral_groups_cache[rock_name][group] = {
                "minerals": minerals,
                "min": min_val,
                "max": max_val
            }

    def _sample_mineralogy(self, rock_name: str, number: int):
        has_groups = rock_name in IsotropicRocks._mineral_groups_cache

        if not has_groups:
            # ---- Modus A: flach ----
            mineral_limits = IsotropicRocks._mineralogy_cache[rock_name]
            mins = [v[0] for v in mineral_limits.values()]
            maxs = [v[1] for v in mineral_limits.values()]
            minerals = list(mineral_limits.keys())

            comp = self._sample_bounded_simplex_batch(mins, maxs, number)
            return minerals, comp

        # ---- Modus B: gruppiert ----
        groups = IsotropicRocks._mineral_groups_cache[rock_name]
        mineral_limits = IsotropicRocks._mineralogy_cache[rock_name]

        # 1️⃣ Gruppen + freie Minerale
        group_names = list(groups.keys())
        group_mins = [groups[g]["min"] for g in group_names]
        group_maxs = [groups[g]["max"] for g in group_names]

        free_minerals = [
            m for m in mineral_limits
            if not any(m in groups[g]["minerals"] for g in groups)
        ]

        free_mins = [mineral_limits[m][0] for m in free_minerals]
        free_maxs = [mineral_limits[m][1] for m in free_minerals]

        labels = group_names + free_minerals
        mins = group_mins + free_mins
        maxs = group_maxs + free_maxs

        # 2️⃣ Top-Level-Sampling
        top_comp = self._sample_bounded_simplex_batch(mins, maxs, number)

        # 3️⃣ Gruppen intern auflösen
        mineral_list = []
        mineral_comp = []

        for i, label in enumerate(labels):
            frac = top_comp[:, i]

            if label in groups:
                minerals = groups[label]["minerals"]
                n = len(minerals)
                split = self.rng.dirichlet(np.ones(n), size=number)

                for j, m in enumerate(minerals):
                    mineral_list.append(m)
                    mineral_comp.append(frac * split[:, j])
            else:
                mineral_list.append(label)
                mineral_comp.append(frac)

        comp = np.vstack(mineral_comp).T
        return mineral_list, comp

    def _sample_bounded_simplex_batch(self, min_vals, max_vals, number):
        min_vals = np.asarray(min_vals, dtype=float)
        max_vals = np.asarray(max_vals, dtype=float)
        n = len(min_vals)

        # --- Consistency ---
        if min_vals.sum() > 1:
            raise ValueError("Sum of minimum fractions > 1")

        if max_vals.sum() < 1:
            raise ValueError("Sum of maximum fractions < 1")

        span = max_vals - min_vals
        samples = np.zeros((number, n))

        for i in range(number):
            remaining = 1.0 - min_vals.sum()
            order = np.arange(n)

            self.rng.shuffle(order)
            x = np.zeros(n)

            for j in order[:-1]:
                upper = min(span[j], remaining)
                val = self.rng.uniform(0, upper)
                x[j] = val
                remaining -= val

            x[order[-1]] = remaining
            samples[i] = min_vals + x

        return samples

    def _collect_mineral_data(self, list_minerals, number):
        _mineral_data = []
        for index, mineral in enumerate(list_minerals):
            data_init = MineralDataGeneration(
                name=mineral, n_datapoints=number, variability=self.variability, uncertainty=self.uncertainty)
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
        arrays = [data_mineral[i][property].to_numpy() for i in range(len(list_minerals))]
        return np.vstack(arrays)

    def _extract_element_data(self, data_minerals, list_elements):
        seen = set(list_elements)
        ordered = list(list_elements)

        for dataset in data_minerals:
            for key in dataset.keys():
                if key.startswith("chemistry."):
                    element = key[len("chemistry."):]
                    if element not in seen:
                        seen.add(element)
                        ordered.append(element)

        return ordered

    def _extract_oxide_data(self, data_minerals, list_oxides):
        seen = set(list_oxides)
        ordered = list(list_oxides)

        for dataset in data_minerals:
            for key in dataset.keys():
                if key.startswith("compounds."):
                    oxide = key[len("compounds."):]
                    if oxide not in seen:
                        seen.add(oxide)
                        ordered.append(oxide)

            # Spezialfall Pyrit
            if dataset["mineral"][0] == "Py":
                for oxide in ("Fe2O3", "SO3"):
                    if oxide not in seen:
                        seen.add(oxide)
                        ordered.append(oxide)
                if "FeS2" in seen:
                    seen.remove("FeS2")
                    ordered = [o for o in ordered if o != "FeS2"]

        return ordered

    def _update_chemistry_data(self, _bulk_data, data_minerals, data_composition, element, number):
        n_minerals = len(data_minerals)
        helper = np.zeros((n_minerals, number))
        key_element = "chemistry." + element

        for i, dataset in enumerate(data_minerals):
            if key_element in dataset:
                helper[i, :] = dataset[key_element].to_numpy()

        bulk_values = np.sum(data_composition * helper.T, axis=1)
        _bulk_data["w." + element] = bulk_values

        return _bulk_data

    def _update_oxide_data(self, bulk_data, list_oxides):
        for oxide in list_oxides:
            cation, anion = self.rock_gen._get_elements_of_compound(compound=oxide)
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

    def collect_geophysical_properties(self, _helper_bulk_data, rho_f, n):
        # Update bulk density data
        (_helper_bulk_data["rho"], _helper_bulk_data["rho_s"],
         _helper_bulk_data["rho_f"]) = self.geophysics.calculate_bulk_density_data(
            v_phi=_helper_bulk_data["porosity"], val_rho=_helper_bulk_data["rho"], val_rho_f=rho_f,
            val_n=n)
        # Update bulk seismic velocity data
        _helper_bulk_data["K"] = _helper_bulk_data["K"]*(1 - _helper_bulk_data["porosity"])**self.alpha_K
        _helper_bulk_data["G"] = _helper_bulk_data["G"]*(1 - _helper_bulk_data["porosity"])**self.alpha_G
        (_helper_bulk_data["vP"], _helper_bulk_data["vS"],
         _helper_bulk_data["vP/vS"]) = self.geophysics.calculate_seismic_velocities(
            val_K=_helper_bulk_data["K"], val_G=_helper_bulk_data["G"], val_rho=_helper_bulk_data["rho"])
        # Update elastic parameter data
        (_helper_bulk_data["E"], _helper_bulk_data["poisson"],
         _helper_bulk_data["lame"]) = self.geophysics.calculate_elastic_parameter_data(
            val_K=_helper_bulk_data["K"], val_G=_helper_bulk_data["G"])

        return _helper_bulk_data

    def collect_initial_compositional_data(self, list_minerals, n):
        _helper_elements = []
        _helper_oxides = []
        _mineral_data = self._collect_mineral_data(list_minerals=list_minerals, number=n)
        _helper_elements = self._extract_element_data(data_minerals=_mineral_data, list_elements=_helper_elements)
        _helper_oxides = self._extract_oxide_data(data_minerals=_mineral_data, list_oxides=_helper_oxides)

        return _mineral_data, _helper_elements, _helper_oxides

    def collect_initial_bulk_data(self, list_minerals, _mineral_data, _helper_composition):
        _helper_bulk_data = {}
        for property in ["rho", "K", "G", "GR", "PE"]:
            _helper_property = self._extract_mineral_property_data(
                list_minerals=list_minerals, data_mineral=_mineral_data, property=property)
            _helper_bulk_data[property] = np.sum(_helper_composition*_helper_property.T, axis=1)

        return _helper_bulk_data

    def update_compositional_bulk_data(
            self, _helper_elements, _helper_bulk_data, _mineral_data, _helper_composition, _helper_oxides,
            _helper_mineral_amounts, n):
        for element in _helper_elements:
            _helper_bulk_data = self._update_chemistry_data(
                _bulk_data=_helper_bulk_data, data_minerals=_mineral_data, data_composition=_helper_composition,
                element=element, number=n)
        # Update oxide data
        _helper_bulk_data = self._update_oxide_data(bulk_data=_helper_bulk_data, list_oxides=_helper_oxides)
        # Update rock composition data
        _helper_bulk_data = self._assign_mineral_amounts(
            bulk_data=_helper_bulk_data, data_amounts=_helper_mineral_amounts)

        return _helper_bulk_data

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
        mineral_limits = IsotropicRocks._mineralogy_cache[self.name]
        _limits = {"lower": [], "upper": []}

        for mineral, values in mineral_limits.items():
            _limits["lower"].append(values[0])
            _limits["upper"].append(values[1])

        list_minerals = list(IsotropicRocks._mineralogy_cache[self.name].keys())
        _properties = ["rho", "vP", "vS", "K", "G", "GR", "PE"]
        _bulk_data = {}
        # Collect mineralogical composition data
        n_minerals = len(list_minerals)
        _helper_composition = np.zeros((number, n_minerals))
        _helper_mineral_amounts = {mineral: np.zeros(number) for mineral in list_minerals}

        _helper_composition = self._sample_bounded_simplex_batch(
            min_vals=_limits["lower"], max_vals=_limits["upper"], number=number)
        _helper_mineral_amounts = {mineral: _helper_composition[:, j] for j, mineral in enumerate(list_minerals)}

        # Collect mineral data
        _mineral_data, _helper_elements, _helper_oxides = self.collect_initial_compositional_data(
            list_minerals=list_minerals, n=number)
        # Collect bulk data
        _helper_bulk_data = self.collect_initial_bulk_data(
            list_minerals=list_minerals, _mineral_data=_mineral_data, _helper_composition=_helper_composition)
        # Assign porosity data
        _helper_bulk_data["porosity"] = porosity
        # Collect geophysical data
        _helper_bulk_data = self.collect_geophysical_properties(
            _helper_bulk_data=_helper_bulk_data, rho_f=density_fluid, n=number)
        # Update chemistry data
        _helper_bulk_data = self.update_compositional_bulk_data(
            _helper_elements=_helper_elements, _helper_bulk_data=_helper_bulk_data, _mineral_data=_mineral_data,
            _helper_composition=_helper_composition, _helper_oxides=_helper_oxides,
            _helper_mineral_amounts=_helper_mineral_amounts, n=number)
        # Conversion to a pandas dataframe object
        _bulk_data = pd.DataFrame(_helper_bulk_data)

        return _bulk_data

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = IsotropicRocks(name="Sandstone", random_seed=42).generate_dataset(number=10)