#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		anisotropic_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		01.02.2026

#-----------------------------------------------

"""
Module: anisotropic_rocks.py
This module controls the generation of synthetic data for anisotropic rocks.
"""

# PACKAGES
import numpy as np
import pandas as pd
from pathlib import Path

# MODULES
from ..rocks.common import CommonRockFunctions
from ..rocks.common import ElasticCalibration

# Code
BASE_PATH = Path(__file__).resolve().parents[2]
DATA_PATH = BASE_PATH / "data_rocks"

class AnisotropicRocks:
    _yaml_cache = {}
    _mineralogy_cache = {}
    _mineral_groups_cache = {}
    _rocks = {"Shale"}

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
        self.cache = {}
        self.class_commonrockfunctions = CommonRockFunctions()

    def generate_dataset(
            self, number: int = 1, fluid: str = "water", density_fluid=None, element_constraints=None, *,
            porosity=None, mineral_comp=None, additional_assemblage=None, comp_constrained=False,
            mineral_constrained=None, porosity_constrained=None, elastic_transformation=False,
            data_smpl: dict = None, porosity_transformation=True) -> None:
        if additional_assemblage is not None and element_constraints is not None:
            raise ValueError("additional_assemblage is not supported together with element_constraints.")
        if comp_constrained and mineral_comp is not None:
            raise ValueError("mineral_comp (fixed) cannot be used together with comp_constrained (interval sampling).")
        if mineral_comp is not None and element_constraints is not None:
            raise ValueError("element_constraints not allowed with fixed mineral_comp.")
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

        if not comp_constrained:
            # YAML processing
            (data_yaml, AnisotropicRocks._yaml_cache, AnisotropicRocks._mineralogy_cache,
             AnisotropicRocks._mineral_groups_cache) = self.class_commonrockfunctions._load_yaml(
                rock_name=self.name, _yaml_cache=AnisotropicRocks._yaml_cache,
                _mineralogy_cache=AnisotropicRocks._mineralogy_cache,
                _mineral_groups_cache=AnisotropicRocks._mineral_groups_cache, _data_path=self.data_path)
            # Porosity definition
            if porosity is None:
                min_porosity = data_yaml["physical_properties"]["porosity"]["min"]
                max_porosity = data_yaml["physical_properties"]["porosity"]["max"]
                porosity = self.rng.uniform(min_porosity, max_porosity, number)
            else:
                if not (0 < porosity < 1):
                    raise ValueError("Porosity must be between 0 and 1.")
                porosity = np.full(number, porosity)
            # Mineral sampling
            mineral_limits = AnisotropicRocks._mineralogy_cache[self.name]
            list_minerals = list(AnisotropicRocks._mineralogy_cache[self.name].keys())
        else:
            # Porosity definition
            if porosity_constrained is None:
                raise ValueError("porosity_constrained must be provided when comp_constrained=True.")

            min_porosity = porosity_constrained["min"]
            max_porosity = porosity_constrained["max"]

            if not (0.0 <= min_porosity < max_porosity <= 1.0):
                raise ValueError("Invalid porosity_constrained interval.")

            porosity = self.rng.uniform(min_porosity, max_porosity, number)
            # Mineral sampling
            if mineral_constrained is None:
                raise ValueError("mineral_constrained must be provided when comp_constrained=True.")
            for m, (lo, hi) in mineral_constrained.items():
                if not (0.0 <= lo <= hi <= 1.0):
                    raise ValueError(f"Invalid bounds for mineral {m}: ({lo}, {hi})")
            if sum(v[0] for v in mineral_constrained.values()) > 1.0:
                raise ValueError("Sum of lower bounds exceeds 1.")
            if sum(v[1] for v in mineral_constrained.values()) < 1.0:
                raise ValueError("Sum of upper bounds below 1.")

            mineral_limits = mineral_constrained
            list_minerals = list(mineral_constrained.keys())

        _limits = {"lower": [], "upper": []}
        for mineral, values in mineral_limits.items():
            _limits["lower"].append(values[0])
            _limits["upper"].append(values[1])

        # Consider additonal mineral assemblage
        if additional_assemblage is not None:
            _limits, list_minerals = self.class_commonrockfunctions.consider_additional_assemblage_data(
                additional_assemblage=additional_assemblage, _limits=_limits, list_minerals=list_minerals)
        # Collect mineralogical composition data
        if mineral_comp is None:
            _helper_composition, _helper_mineral_amounts = self.class_commonrockfunctions._calculate_chemical_amounts(
                list_minerals=list_minerals, number=number, _limits=_limits, _rng=self.rng,
                _variability=self.variability, _uncertainty=self.uncertainty, element_constraints=element_constraints)
        else:
            if set(mineral_comp.keys()) != set(list_minerals):
                raise ValueError("mineral_comp must define all minerals.")
            if not np.isclose(sum(mineral_comp.values()), 1.0):
                raise ValueError("Mineral fractions must sum to 1.")
            if any(v < 0 for v in mineral_comp.values()):
                raise ValueError("Mineral fractions must be >= 0.")

            list_amounts = []
            _helper_mineral_amounts = {}
            for mineral in list_minerals:
                amount = mineral_comp[mineral]
                list_amounts.append(amount)
                _helper_mineral_amounts[mineral] = np.ones(number)*amount
            _helper_composition = np.ones((number, len(list_minerals)))*list_amounts
        # Collect mineral data
        if element_constraints is not None:
            _mineral_data, _helper_elements, _helper_oxides = (
                self.class_commonrockfunctions.collect_initial_compositional_data(
                    list_minerals=list_minerals, n=1, _variability=self.variability, _uncertainty=self.uncertainty))
            _mineral_data = [pd.concat([df]*number, ignore_index=True) for df in _mineral_data]
        else:
            _mineral_data, _helper_elements, _helper_oxides = (
                self.class_commonrockfunctions.collect_initial_compositional_data(
                    list_minerals=list_minerals, n=number, _variability=self.variability,
                    _uncertainty=self.uncertainty))
        # Collect bulk data
        _helper_bulk_data = self.class_commonrockfunctions.collect_initial_bulk_data(
            list_minerals=list_minerals, _mineral_data=_mineral_data, _helper_composition=_helper_composition)
        # Assign porosity data
        _helper_bulk_data["porosity"] = porosity
        # Collect geophysical data
        _helper_bulk_data = self.class_commonrockfunctions.collect_geophysical_properties(
            _helper_bulk_data=_helper_bulk_data, rho_f=density_fluid, n=number, alpha_K=self.alpha_K,
            alpha_G=self.alpha_G)
        # Consider structural effects
        if elastic_transformation == True:
            if porosity_transformation == True:
                data_rho_phi_trans = ElasticCalibration.transform_bulk_density_and_porosity(
                    rho_smpl=data_smpl["rho"], rho_s_ref=_helper_bulk_data["rho_s"],
                    rho_f_ref=_helper_bulk_data["rho_f"], porosity_ref=_helper_bulk_data["porosity"])
                for param in ["rho", "porosity"]:
                    _helper_bulk_data[param] = data_rho_phi_trans[param + "_opt"]
            data_elast_transf = ElasticCalibration.transform_elastic_moduli(
                rho_smpl=data_smpl["rho"], vp_smpl=data_smpl["vP"], vs_smpl=data_smpl["vS"],
                k_ref=_helper_bulk_data["K"], g_ref=_helper_bulk_data["G"], rho_ref=_helper_bulk_data["rho"])
            for param in ["K", "G", "E", "poisson", "lame", "vP", "vS"]:
                _helper_bulk_data[param] = data_elast_transf[param + "_opt"]
            _helper_bulk_data["vP/vS"] = data_elast_transf["vP_opt"]/data_elast_transf["vS_opt"]
            _helper_bulk_data["SDI"] = data_elast_transf["SDI"]
            _helper_bulk_data["fK"] = data_elast_transf["fK"]
            _helper_bulk_data["fG"] = data_elast_transf["fG"]
        # Update chemistry data
        _helper_bulk_data = self.class_commonrockfunctions.update_compositional_bulk_data(
            _helper_elements=_helper_elements, _helper_bulk_data=_helper_bulk_data, _mineral_data=_mineral_data,
            _helper_composition=_helper_composition, _helper_oxides=_helper_oxides,
            _helper_mineral_amounts=_helper_mineral_amounts, n=number)
        # Conversion to a pandas dataframe object
        _bulk_data = pd.DataFrame(_helper_bulk_data)

        return _bulk_data

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = AnisotropicRocks(name="Shale", random_seed=42).generate_dataset(number=10)