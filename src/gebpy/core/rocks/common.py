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

# MODULES
from ..chemistry.common import PeriodicSystem

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
            "U": PeriodicSystem(name="U").get_data()}

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
        pass

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