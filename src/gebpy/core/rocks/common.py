#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		common.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.01.2026

#-----------------------------------------------

"""
Module: common.py
This module contains several routines that are commonly used by the different rock-related modules.
"""
import pathlib
# PACKAGES
import re, yaml
import numpy as np
import pandas as pd
from sqlalchemy import BOOLEAN

# MODULES
from ..chemistry.common import PeriodicSystem
from ..minerals.synthesis import MineralDataGeneration
from ..physics.common import Geophysics

class RockGeneration:
    def __init__(self):
        self.elements = {
            "H": PeriodicSystem(name="H").get_data(),
            "C": PeriodicSystem(name="C").get_data(),
            "O": PeriodicSystem(name="O").get_data(),
            "Na": PeriodicSystem(name="Na").get_data(),
            "Mg": PeriodicSystem(name="Mg").get_data(),
            "Al": PeriodicSystem(name="Al").get_data(),
            "Si": PeriodicSystem(name="Si").get_data(),
            "S": PeriodicSystem(name="S").get_data(),
            "Cl": PeriodicSystem(name="Cl").get_data(),
            "K": PeriodicSystem(name="K").get_data(),
            "Ca": PeriodicSystem(name="Ca").get_data(),
            "Mn": PeriodicSystem(name="Mn").get_data(),
            "Fe": PeriodicSystem(name="Fe").get_data(),
            "Ni": PeriodicSystem(name="Ni").get_data(),
            "U": PeriodicSystem(name="U").get_data(),
            "Zn": PeriodicSystem(name="Zn").get_data(),
            "Pb": PeriodicSystem(name="Pb").get_data(),
            "Cu": PeriodicSystem(name="Cu").get_data()}

    def _parse_formula(self, formula: str):
        pattern = r"([A-Z][a-z]?)(\d*)"
        matches = re.findall(pattern, formula)

        composition = {}
        for elem, amount in matches:
            amount = int(amount) if amount else 1
            composition[elem] = composition.get(elem, 0) + amount

        return composition

    def _get_elements_of_compound(self, compound: str) -> str:
        elements = re.findall(r"[A-Z][a-z]?", compound)
        first = elements[0]
        last = elements[-1]

        return  first, last

    def _get_cation_element(self, oxide: str) -> str:
        first = oxide[0]
        if len(oxide) > 1 and oxide[1].islower():
            return oxide[:2]

        return first

    def _get_anion_element(self, compound: str) -> str:
        elements = re.findall(r"[A-Z][a-z]?", compound)
        last = elements[-1]
        return last

    def _determine_oxide_conversion_factors(self):
        list_oxides = [
            "H2O", "CO", "CO2", "Na2O", "MgO", "Al2O3", "SiO2", "Cl2O", "K2O", "CaO", "MnO", "Mn2O3", "MnO2", "MnO3",
            "Mn2O7", "FeO", "Fe2O3", "FeO3", "NiO", "Ni2O3", "TiO2", "Ti2O3", "VO", "V2O3", "VO2", "V2O10", "CrO",
            "Cr2O3", "CrO3", "CoO", "Co2O3", "Cu2O", "CuO", "ZnO", "GeO2", "As2O3", "As2O10", "ZrO2", "Nb2O3", "Nb2O10",
            "MoO", "Mo2O3", "MoO2", "Mo2O10", "MoO3", "CdO", "SnO", "SnO2", "Sb2O3", "Sb2O10", "TeO2", "TeO3", "Ta2O10",
            "WO", "W2O3", "WO2", "W2O10", "WO3", "Au2O", "Au2O3", "Tl2O", "Tl2O3", "PbO", "PbO2", "Bi2O3", "Bi2O10",
            "U2O3", "UO2", "U2O10", "UO3", "Nb2O5", "Ta2O5", "SO", "SO2", "SO3"]
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

class CommonRockFunctions:
    def __init__(self):
        self.geophysics = Geophysics()
        self.rock_gen = RockGeneration()
        self.conversion_factors = self.rock_gen._determine_oxide_conversion_factors()

    # YAML processing
    def _compile_mineralogy(self, rock_name: str, mineralogy_dict: dict, _mineralogy_cache: dict):
        """
        Extracts and compiles all chemistry formulas from the YAML file.
        Stores the compiled ASTs in the global formula cache.
        """
        if rock_name not in _mineralogy_cache:
            _mineralogy_cache[rock_name] = {}

        for element, entry in mineralogy_dict.items():
            mineral = element
            interval = list(entry.values())
            lower_limit = interval[0]
            upper_limit = interval[1]
            compiled = [lower_limit, upper_limit]
            _mineralogy_cache[rock_name][mineral] = compiled

        return _mineralogy_cache

    def _compile_mineral_groups(self, rock_name: str, group_dict: dict, _mineral_groups_cache: dict):
        if rock_name not in _mineral_groups_cache:
            _mineral_groups_cache[rock_name] = {}

        for group, entry in group_dict.items():
            minerals = entry["minerals"]
            min_val = entry["min"]
            max_val = entry["max"]

            _mineral_groups_cache[rock_name][group] = {
                "minerals": minerals,
                "min": min_val,
                "max": max_val
            }

        return _mineral_groups_cache

    def _load_yaml(
            self, rock_name: str, _yaml_cache: dict, _mineralogy_cache: dict, _mineral_groups_cache: dict,
            _data_path: pathlib.WindowsPath) -> dict:
        """
        Extracts and compiles all chemistry formulas from the YAML file. Stores the compiled ASTs in the global formula
        cache.
        >> Parameters
        ----------
            rock_name: str, Name of the rock (YAML filename without extension).
        >> Returns
        -------
            data: dict, Parsed YAML content.
        >>Raises
        ------
            FileNotFoundError: If the YAML file does not exist.
        """
        if rock_name in _yaml_cache:
            return (_yaml_cache[rock_name], _yaml_cache, _mineralogy_cache, _mineral_groups_cache)

        yaml_file = _data_path/f"{rock_name}.yaml"
        if not yaml_file.exists():
            raise FileNotFoundError(f"No YAML file found for {rock_name}.")

        with open(yaml_file, "r") as f:
            data = yaml.safe_load(f)

        if "mineralogy" in data and rock_name not in _mineralogy_cache:
            _mineralogy_cache = self._compile_mineralogy(
                rock_name=rock_name, mineralogy_dict=data["mineralogy"], _mineralogy_cache=_mineralogy_cache)
        if "mineral_groups" in data:
            _mineral_groups_cache = self._compile_mineral_groups(
                rock_name=rock_name, group_dict=data["mineral_groups"], _mineral_groups_cache=_mineral_groups_cache)

        _yaml_cache[rock_name] = data

        return (data, _yaml_cache, _mineralogy_cache, _mineral_groups_cache)

    # Mineral sampling
    def _collect_mineral_data(self, list_minerals, number, _variability, _uncertainty):
        _mineral_data = []
        for index, mineral in enumerate(list_minerals):
            data_init = MineralDataGeneration(
                name=mineral, n_datapoints=number, variability=_variability, uncertainty=_uncertainty)
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
        SULFIDE_OXIDE_MAP = {
            "Py": {"add": ["Fe2O3", "SO3"], "remove": ["FeS2"]},
            "Ccp": {"add": ["Cu2O", "Fe2O3", "SO3"], "remove": ["FeS", "CuS"]},
            "Bn": {"add": ["Cu2O", "Fe2O3", "SO3"], "remove": ["FeS", "Cu2S", "CuS"]},
            "Gn": {"add": ["PbO", "SO3"], "remove": ["PbS"]},
            "Sp": {"add": ["ZnO", "SO3"], "remove": ["ZnS"]}
        }

        seen = set(list_oxides)
        ordered = list(list_oxides)

        for dataset in data_minerals:
            for key in dataset.keys():
                if key.startswith("compounds."):
                    oxide = key[len("compounds."):]
                    if oxide not in seen:
                        seen.add(oxide)
                        ordered.append(oxide)

            mineral = dataset["mineral"][0]
            if mineral in SULFIDE_OXIDE_MAP:
                mapping = SULFIDE_OXIDE_MAP[mineral]
                for oxide in mapping["add"]:
                    if oxide not in seen:
                        seen.add(oxide)
                        ordered.append(oxide)
                for sulfide in mapping["remove"]:
                    if sulfide in seen:
                        seen.remove(sulfide)
                        ordered.remove(sulfide)

        return ordered

    def collect_initial_compositional_data(self, list_minerals, n, _variability, _uncertainty):
        _helper_elements = []
        _helper_oxides = []
        _mineral_data = self._collect_mineral_data(
            list_minerals=list_minerals, number=n, _variability=_variability, _uncertainty=_uncertainty)
        _helper_elements = self._extract_element_data(data_minerals=_mineral_data, list_elements=_helper_elements)
        _helper_oxides = self._extract_oxide_data(data_minerals=_mineral_data, list_oxides=_helper_oxides)

        return _mineral_data, _helper_elements, _helper_oxides

    def _sample_bounded_simplex_batch(self, min_vals, max_vals, number, _rng):
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

            _rng.shuffle(order)
            x = np.zeros(n)

            for j in order[:-1]:
                upper = min(span[j], remaining)
                val = _rng.uniform(0, upper)
                x[j] = val
                remaining -= val

            x[order[-1]] = remaining
            samples[i] = min_vals + x

        return samples

    def _compute_bulk_element(self, element, mineral_data, composition):
        n_samples = composition.shape[0]
        values = np.zeros(n_samples, dtype=float)
        key = "chemistry." + element

        for j, dataset in enumerate(mineral_data):
            if key not in dataset:
                continue

            chem = dataset[key].to_numpy(dtype=float)

            # Broadcast scalar chemistry
            if chem.size == 1:
                chem = np.full(n_samples, chem[0], dtype=float)

            # If chemistry vector but only one sample: take first value deterministically
            elif n_samples == 1 and chem.size > 1:
                chem = np.array([chem[0]], dtype=float)

            elif chem.size != n_samples:
                raise ValueError(
                    f"Shape mismatch for element '{element}': "
                    f"chemistry length {chem.size} vs. composition samples {n_samples}"
                )

            values += composition[:, j]*chem

        return values

    def _calculate_chemical_amounts(
            self, list_minerals, number, _limits, _rng, _variability, _uncertainty, element_constraints=None):
        if element_constraints == None:
            n_minerals = len(list_minerals)
            _helper_composition = np.zeros((number, n_minerals))
            _helper_mineral_amounts = {mineral: np.zeros(number) for mineral in list_minerals}
            _helper_composition = self._sample_bounded_simplex_batch(
                min_vals=_limits["lower"], max_vals=_limits["upper"], number=number, _rng=_rng)
            _helper_mineral_amounts = {mineral: _helper_composition[:, j] for j, mineral in enumerate(list_minerals)}
        elif element_constraints != None:
            valid_comp = []
            attempts = 0
            max_attempts = number*100

            # Mineral data wird für die Elementberechnung benötigt
            _mineral_data, _, _ = self.collect_initial_compositional_data(
                list_minerals=list_minerals, n=1, _variability=_variability, _uncertainty=_uncertainty
            )

            while len(valid_comp) < number:
                comp = self._sample_bounded_simplex_batch(
                    min_vals=_limits["lower"],
                    max_vals=_limits["upper"],
                    number=1, _rng=_rng
                )

                is_valid = True
                if element_constraints:
                    for el, (lo, hi) in element_constraints.items():
                        bulk_val = self._compute_bulk_element(
                            el, _mineral_data, comp
                        )[0]
                        if not (lo <= bulk_val <= hi):
                            is_valid = False
                            break

                if is_valid:
                    valid_comp.append(comp[0])

                attempts += 1
                if attempts%1000 == 0:
                    print(f"Acceptance rate: {len(valid_comp)/attempts:.3f}")
                if attempts > max_attempts:
                    raise RuntimeError(
                        "Element constraints too restrictive for given mineralogy."
                    )

            _helper_composition = np.vstack(valid_comp)
            _helper_mineral_amounts = {
                mineral: _helper_composition[:, j]
                for j, mineral in enumerate(list_minerals)
            }

        return _helper_composition, _helper_mineral_amounts

    # Collecting mineral data
    def _extract_mineral_property_data(self, list_minerals, data_mineral, property):
        arrays = [data_mineral[i][property].to_numpy() for i in range(len(list_minerals))]
        return np.vstack(arrays)

    def collect_initial_bulk_data(self, list_minerals, _mineral_data, _helper_composition):
        _helper_bulk_data = {}
        for property in ["rho", "K", "G", "GR", "PE"]:
            _helper_property = self._extract_mineral_property_data(
                list_minerals=list_minerals, data_mineral=_mineral_data, property=property)
            _helper_bulk_data[property] = np.sum(_helper_composition*_helper_property.T, axis=1)

        return _helper_bulk_data

    # Collecting geophysical bulk data (no anisotropy consideration)
    def collect_geophysical_properties(self, _helper_bulk_data, rho_f, n, alpha_K, alpha_G):
        # Update bulk density data
        (_helper_bulk_data["rho"], _helper_bulk_data["rho_s"],
         _helper_bulk_data["rho_f"]) = self.geophysics.calculate_bulk_density_data(
            v_phi=_helper_bulk_data["porosity"], val_rho=_helper_bulk_data["rho"], val_rho_f=rho_f,
            val_n=n)
        # Update bulk seismic velocity data
        _helper_bulk_data["K"] = _helper_bulk_data["K"]*(1 - _helper_bulk_data["porosity"])**alpha_K
        _helper_bulk_data["G"] = _helper_bulk_data["G"]*(1 - _helper_bulk_data["porosity"])**alpha_G
        (_helper_bulk_data["vP"], _helper_bulk_data["vS"],
         _helper_bulk_data["vP/vS"]) = self.geophysics.calculate_seismic_velocities(
            val_K=_helper_bulk_data["K"], val_G=_helper_bulk_data["G"], val_rho=_helper_bulk_data["rho"])
        # Update elastic parameter data
        (_helper_bulk_data["E"], _helper_bulk_data["poisson"],
         _helper_bulk_data["lame"]) = self.geophysics.calculate_elastic_parameter_data(
            val_K=_helper_bulk_data["K"], val_G=_helper_bulk_data["G"])

        return _helper_bulk_data

    # Update compositional bulk data
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

    def consider_additional_assemblage_data(self, additional_assemblage, _limits, list_minerals):
        extra_fraction = additional_assemblage["volume_fraction"]
        list_extra_minerals = list(additional_assemblage["mineralogy"].keys())
        list_extra_amounts = np.array(list(additional_assemblage["mineralogy"].values()))

        if any(m in list_minerals for m in list_extra_minerals):
            raise ValueError("Additional minerals already present in host mineralogy.")
        if not (0.0 < extra_fraction < 1.0):
            raise ValueError("volume_fraction must be in the interval [0, 1).")
        if not np.isclose(np.sum(list_extra_amounts), 1.0):
            raise ValueError("Mineral amounts of additional mineral assemblage must sum to 1.")
        if any(v < 0 for v in list_extra_amounts):
            raise ValueError("Mineral amounts of additional mineral assemblage must be non-negative.")
        if any(v > 1 for v in list_extra_amounts):
            raise ValueError("Mineral amounts of additional mineral assemblage must be <= 1.")

        list_minerals.extend(list_extra_minerals)
        list_extra_amounts_corrected = list_extra_amounts*extra_fraction
        fraction_correction = 1 - extra_fraction
        for key, values in _limits.items():
            _limits[key] = list(np.array(values)*fraction_correction)
            _limits[key].extend(list_extra_amounts_corrected)

        return _limits, list_minerals