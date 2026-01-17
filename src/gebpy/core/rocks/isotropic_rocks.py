#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		isotropic_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.01.2026

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
from ..rocks.common import RockGeneration, CommonRockFunctions
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

    def _extract_mineral_property_data(self, list_minerals, data_mineral, property):
        arrays = [data_mineral[i][property].to_numpy() for i in range(len(list_minerals))]
        return np.vstack(arrays)

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

    def generate_dataset(
            self, number: int = 1, fluid: str = "water", density_fluid=None, element_constraints=None) -> None:
        if element_constraints:
            for el, (lo, hi) in element_constraints.items():
                if not (0.0 <= lo < hi <= 1.0):
                    raise ValueError(
                        f"Element constraint for {el} must be given as fraction (0–1), "
                        f"got ({lo}, {hi})."
                    )

        if density_fluid is None:
            if fluid == "water":
                density_fluid = 1000
            elif fluid == "oil":
                density_fluid = 800
            elif fluid == "natural gas":
                density_fluid = 750

        # YAML processing
        (data_yaml, IsotropicRocks._yaml_cache, IsotropicRocks._mineralogy_cache,
         IsotropicRocks._mineral_groups_cache) = CommonRockFunctions()._load_yaml(
            rock_name=self.name, _yaml_cache=IsotropicRocks._yaml_cache,
            _mineralogy_cache=IsotropicRocks._mineralogy_cache,
            _mineral_groups_cache=IsotropicRocks._mineral_groups_cache, _data_path=self.data_path)
        # Porosity definition
        min_porosity = data_yaml["physical_properties"]["porosity"]["min"]
        max_porosity = data_yaml["physical_properties"]["porosity"]["max"]
        porosity = self.rng.uniform(min_porosity, max_porosity, number)
        # Mineral sampling
        mineral_limits = IsotropicRocks._mineralogy_cache[self.name]
        _limits = {"lower": [], "upper": []}

        for mineral, values in mineral_limits.items():
            _limits["lower"].append(values[0])
            _limits["upper"].append(values[1])

        list_minerals = list(IsotropicRocks._mineralogy_cache[self.name].keys())
        _bulk_data = {}
        # Collect mineralogical composition data
        _helper_composition, _helper_mineral_amounts = CommonRockFunctions()._calculate_chemical_amounts(
            list_minerals=list_minerals, number=number, _limits=_limits, _rng=self.rng,
            _variability=self.variability, _uncertainty=self.uncertainty, element_constraints=element_constraints)
        # Collect mineral data
        if element_constraints != None:
            _mineral_data, _helper_elements, _helper_oxides = CommonRockFunctions().collect_initial_compositional_data(
                list_minerals=list_minerals, n=1, _variability=self.variability, _uncertainty=self.uncertainty)
            _mineral_data = [pd.concat([df]*number, ignore_index=True) for df in _mineral_data]
        else:
            _mineral_data, _helper_elements, _helper_oxides = CommonRockFunctions().collect_initial_compositional_data(
                list_minerals=list_minerals, n=number, _variability=self.variability, _uncertainty=self.uncertainty)
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