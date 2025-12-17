#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		synthesis.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: synthesis.py
Generates the synthetic mineral data based on the settings of the previous configuration.
"""

# PACKAGES
import pandas as pd
from typing import Optional, Union
from functools import cached_property

from pyasn1.type.univ import Boolean

# MODULES
from ..minerals.config import DEFAULT_CONFIG
from ..minerals.carbonates import Carbonates
from ..minerals.oxides import Oxides
from ..minerals.phyllosilicates import Phyllosilicates
from ..minerals.sulfides import Sulfides
from ..minerals.tectosilicates import Tectosilicates

# CODE
class MineralDataGeneration:
    """
    Generation of the synthetic/simulated mineral data.
    - name: mineral name (e.g., 'Olivine')
    - n_datapoints: Number of generated data points (> 0)
    - random_seed: integer seed; if None, defaults to 42
    """

    list_carbonates = {
        "Ankerite", "Aragonite", "Azurite", "Calcite", "Cerrusite", "Dolomite", "Ikaite", "Magnesite", "Malachite",
        "Rhodochrosite", "Siderite", "Smithsonite"}
    list_cyclosilicates = {
        "Benitoite", "Beryl", "Cordierite", "Elbaite", "Liddicoatite", "Schorl", "Sekaninaite"}
    list_halides = {"Carnallite", "Fluorite", "Halite", "Sylvite"}
    list_inosilicates = {
        "Actinolite", "Aegirine", "Arfvedsonite", "Augite", "Ca-Amphibole", "Ca-Pyroxene", "Clinopyroxene",
        "Diopside", "Donpeacorite", "Enstatite", "Ferrosilite", "Glaucophane", "Jadeite", "Mg-Fe-Pyroxene",
        "Na-Amphibole", "Na-Pyroxene", "Orthopyroxene", "Spodumene", "Tremolite", "Wollastonite"}
    list_miscellaneous = {"Organic matter"}
    list_nesosilicates = {
        "Al-Garnet", "Almandine", "Andalusite", "Anhadrite", "Ca-Garnet", "Ca-Olivine", "Fayalite", "Forsterite",
        "Grossular", "Kyanite", "Liebenbergite", "Olivine", "Pyrope", "Sillimanite", "Staurolite", "Tephroite",
        "Thorite", "Titanite", "Topaz", "Uvarovite", "Zircon"}
    list_oxides = Oxides._minerals
    list_phosphates = {"Apatite", "Chloroapatite", "Fluoroapatite", "Hydroxyapatite"}
    list_phospides = {"Allabogdanite"}
    list_phyllosilicates = Phyllosilicates._minerals
    list_sorosilicates = {"Epidote", "Gehlenite", "Zoisite"}
    list_sulfates = {
        "Alunite", "Anglesite", "Anhydrite", "Barite", "Celestite", "Chalcanthite", "Gypsum", "Hanksite",
        "Hexahydrite", "Jarosite", "Kainite", "Kieserite", "Scheelite"}
    list_sulfides = {
        "Acanthite", "Bornite", "Cattierite", "Chalcocite", "Chalcopyrite", "Cinnabar", "Cobaltite", "Covellite",
        "Fahlore", "Galena", "Gallite", "Laforetite", "Lenaite", "Marcasite", "Marmatite", "Millerite",
        "Molybdenite", "Orpiment", "Pentlandite", "Pyrite", "Pyrrhotite", "Realgar", "Roquesite", "Sphalerite",
        "Stibnite", "Vaesite"}
    list_tectosilicates = Tectosilicates._minerals

    mineral_groups = {
        "carbonates": list_carbonates, "cyclosilicates": list_cyclosilicates, "halides": list_halides,
        "inosilicates": list_inosilicates, "miscellaneous": list_miscellaneous, "nesosilicates": list_nesosilicates,
        "oxides": list_oxides, "phosphates": list_phosphates, "phospides": list_phospides,
        "phyllosilicates": list_phyllosilicates, "sorosilicates": list_sorosilicates, "sulfates": list_sulfates,
        "sulfides": list_sulfides, "tectosilicates": list_tectosilicates}

    def __init__(self,
                 name: Optional[str] = None,
                 n_datapoints: Optional[int] = None,
                 random_seed: Optional[int] = None,
                 trace_elements: Optional[list] = None,
                 variability: Optional[Boolean] = False,
                 uncertainty: Optional[float] = 1.0
    ) -> None:
        if name is None and n_datapoints is None and random_seed is None:
            self.name = DEFAULT_CONFIG.name
            self.n_datapoints = DEFAULT_CONFIG.n_datapoints
            self.random_seed = DEFAULT_CONFIG.random_seed
            self.trace_elements = []
            self.variability = variability
            self.uncertainty = uncertainty
        else:
            self.name = name
            self.n_datapoints = n_datapoints
            self.random_seed = 42 if random_seed is None else random_seed
            self.trace_elements = [] if trace_elements is None else trace_elements
            self.variability = variability
            self.uncertainty = uncertainty

    def __repr__(self) -> str:
        return (
        f"<MineralDataGeneration name={self.name!r}, "
        f"n_datapoints={self.n_datapoints}, "
        f"random_seed={self.random_seed}>"
        )

    @cached_property
    def mineral_map(self):
        from gebpy_legacy.modules.halides import Halides
        from gebpy_legacy.modules.silicates import (Cyclosilicates, Inosilicates, Nesosilicates, Sorosilicates)
        from gebpy_legacy.modules.phosphates import Phosphates
        from gebpy_legacy.modules.sulfates import Sulfates
        from gebpy_legacy.modules.organics import Organics
        from gebpy_legacy.modules.phospides import Phospides

        return {
            "carbonates": Carbonates,
            "cyclosilicates": Cyclosilicates,
            "halides": Halides,
            "inosilicates": Inosilicates,
            "miscellaneous": Organics,
            "nesosilicates": Nesosilicates,
            "oxides": Oxides,
            "phosphates": Phosphates,
            "phospides": Phospides,
            "phyllosilicates": Phyllosilicates,
            "sorosilicates": Sorosilicates,
            "sulfates": Sulfates,
            "sulfides": Sulfides,
            "tectosilicates": Tectosilicates
        }

    def get_all_minerals(self):
        return self.mineral_groups

    def assign_mineral_group(self,
            name: Optional[str] = None,
            n_datapoints: Optional[int] = None
    ) -> None:
        if name is None:
            name = self.name
            n_datapoints = self.n_datapoints

        for group, list_group in self.mineral_groups.items():
            if name in list_group:
                cls = self.mineral_map[group]
                if cls.__name__ in ("Phyllosilicates", "Tectosilicates", "Oxides", "Carbonates", "Sulfides"):
                    #print(sorted(cls(name=name, random_seed=self.random_seed)._minerals))
                    return cls(
                        name=name, random_seed=self.random_seed, variability=self.variability,
                        uncertainty=self.uncertainty).generate_dataset(number=n_datapoints)
                else:
                    try:
                        return cls(
                            mineral=name, data_type=True, traces_list=self.trace_elements,
                            random_seed=self.random_seed).generate_dataset(
                            number=n_datapoints)
                    except:
                        return cls(mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(
                            number=n_datapoints)
        raise ValueError(f"Unknown mineral '{name}'. Please check available groups.")

    def generate_data(self, as_dataframe: bool = True) -> Union[pd.DataFrame, dict]:
        data_dict = self.assign_mineral_group(self.name, self.n_datapoints)
        if not as_dataframe:
            return data_dict

        data_mineral = self.dict_to_dataframe_fast(data_dict)
        data_mineral = data_mineral.astype("float32", errors="ignore").round(5)
        return data_mineral

    def dict_to_dataframe_fast(self, data: dict) -> pd.DataFrame:
        """
        Convert nested mineral dataset dictionaries (with lists and subdicts)
        into a flat pandas DataFrame. Handles scalar values gracefully by
        repeating them to match the maximum list length.
        """
        columns = {}

        # Schritt 1: maximale Listenlänge bestimmen
        max_len = 1
        for val in data.values():
            if isinstance(val, dict):
                for subval in val.values():
                    if isinstance(subval, list):
                        max_len = max(max_len, len(subval))
            elif isinstance(val, list):
                max_len = max(max_len, len(val))

        # Schritt 2: flaches Dictionary aufbauen
        for key, val in data.items():
            if isinstance(val, dict):
                for subkey, subval in val.items():
                    col_name = f"{key}.{subkey}"
                    if isinstance(subval, list):
                        columns[col_name] = subval
                    else:
                        columns[col_name] = [subval] * max_len
            elif isinstance(val, list):
                columns[key] = val
            else:
                columns[key] = [val] * max_len

        return pd.DataFrame(columns)

# DEFAULT EXAMPLE
data_config = DEFAULT_CONFIG
data_default = MineralDataGeneration(name=data_config.name, n_datapoints=data_config.n_datapoints)
DEFAULT_DATA = data_default.generate_data()