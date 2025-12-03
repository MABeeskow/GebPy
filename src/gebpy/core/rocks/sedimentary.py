#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sedimentary.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		04.12.2025

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

from soupsieve.util import lower
from asteval import Interpreter

# MODULES
from src.gebpy.core.minerals.synthesis import MineralDataGeneration

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

    def _weighted_sum(self, dfs, fractions):
        bulk = None

        for df, frac in zip(dfs, fractions):
            numeric_cols = df.select_dtypes(include=np.number).columns
            weighted = df[numeric_cols]*frac

            if bulk is None:
                bulk = weighted
            else:
                bulk = bulk.add(weighted, fill_value=0)

        return bulk

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
        mineral_amounts = self._sample_bounded_simplex(min_vals=_limits["lower"], max_vals=_limits["upper"])
        _properties = ["rho", "vP", "vS", "K", "G", "GR", "PE"]
        _bulk_data = {}
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

        _bulk_data = self._weighted_sum(dfs=_mineral_data, fractions=mineral_amounts)
        # Density adjustment due to porosity
        _bulk_data["porosity"] = porosity
        _bulk_data["rho_s"] = _bulk_data["rho"]
        _bulk_data["rho"] = (1 - porosity)*_bulk_data["rho"] + porosity*density_fluid
        _bulk_data["rho_f"] = np.ones(number)*density_fluid
        # Value conversion
        val_K = _bulk_data["K"]*1e9
        val_G = _bulk_data["G"]*1e9
        val_rho = _bulk_data["rho"]
        # Seismic velocities
        _bulk_data["vS"] = (val_G/val_rho)**0.5
        _bulk_data["vP"] = ((val_K + 4/3*val_G)/val_rho)**0.5
        _bulk_data["vP/vS"] = _bulk_data["vP"]/_bulk_data["vS"]
        # Elastic moduli
        _bulk_data["E"] = ((9*val_K*val_G)/(3*val_K + val_G))*1e-9
        _bulk_data["poisson"] = (3*val_K - 2*val_G)/(2*(3*val_K + val_G))
        _bulk_data["lame"] = (val_K - (2*val_G)/3)*1e-9
        # Cleanup
        _bulk_data.drop(columns=["nu"], errors="ignore", inplace=True)

        return _bulk_data

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = SedimentaryRocks(name="Sandstone", random_seed=42).generate_dataset(number=10)