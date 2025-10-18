#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		synthesis.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		18.10.2025

#-----------------------------------------------

"""
Module: synthesis.py
Generates the synthetic mineral data based on the settings of the previous configuration.
"""

# PACKAGES
from typing import Optional

# MODULES
from src.gebpy.core.minerals.config import DEFAULT_CONFIG

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
    list_oxides = {
        "Al-Spinel", "Anatase", "Arsenolite", "Au(III)-Oxide", "Bismite", "Boehmite", "Brookite", "Brucite",
        "Cassiterite", "Chromite", "Claudetite", "Cochromite", "Coltan", "Columbite", "Corundum", "Cr-Spinel",
        "Crocoite", "Cuprite", "Cuprospinel", "Diaspore", "Fe-Spinel", "Ferberite", "Ferberite-Huebnerite",
        "Franklinite", "Geikielite", "Gibbsite", "Goethite", "Groutite", "Hematite", "Huebnerite", "Ilmenite",
        "Jacobsite", "Litharge", "Magnesiochromite", "Magnesioferrite", "Magnetite", "Manganite", "Manganochromite",
        "Massicot", "Minium", "Nichromite", "Plattnerite", "Pyrolusite", "Pyrophanite", "Quartz", "Rutile",
        "Scrutinyite", "Senarmontite", "Sphaerobismite", "Spinel", "Tantalite", "Trevorite", "Ulvospinel",
        "Uraninite", "Valentinite", "Wolframite", "Wulfenite", "Zincite", "Zincochromite"}
    list_phosphates = {"Apatite", "Chloroapatite", "Fluoroapatite", "Hydroxyapatite"}
    list_phospides = {"Allabogdanite"}
    list_phyllosilicates = {
        "Annite", "Antigorite", "Biotite", "Chamosite", "Chlorite", "Chrysotile", "Clinochlore", "Eastonite",
        "Glauconite", "Illite", "Kaolinite", "Montmorillonite", "Muscovite", "Nimite", "Nontronite", "Pennantite",
        "Phlogopite", "Pyrophyillite", "Saponite", "Siderophyllite", "Talc", "Vermiculite"}
    list_sorosilicates = {"Epidote", "Gehlenite", "Zoisite"}
    list_sulfates = {
        "Alunite", "Anglesite", "Anhydrite", "Barite", "Celestite", "Chalcanthite", "Gypsum", "Hanksite",
        "Hexahydrite", "Jarosite", "Kainite", "Kieserite", "Scheelite"}
    list_sulfides = {
        "Acanthite", "Bornite", "Cattierite", "Chalcocite", "Chalcopyrite", "Cinnabar", "Cobaltite", "Covellite",
        "Fahlore", "Galena", "Gallite", "Laforetite", "Lenaite", "Marcasite", "Marmatite", "Millerite",
        "Molybdenite", "Orpiment", "Pentlandite", "Pyrite", "Pyrrhotite", "Realgar", "Roquesite", "Sphalerite",
        "Stibnite", "Vaesite"}
    list_tectosilicates = {
        "Alkali Feldspar", "Danburite", "Nepheline", "Orthoclase", "Plagioclase", "Scapolite"}

    def __init__(self,
                 name: Optional[str] = None,
                 n_datapoints: Optional[int] = None,
                 random_seed: Optional[int] = None,
                 trace_elements: Optional[list] = None
    ) -> None:
        if name is None and n_datapoints is None and random_seed is None:
            self.name = DEFAULT_CONFIG.name
            self.n_datapoints = DEFAULT_CONFIG.n_datapoints
            self.random_seed = DEFAULT_CONFIG.random_seed
            self.trace_elements = []
        else:
            self.name = name
            self.n_datapoints = n_datapoints
            self.random_seed = 42 if random_seed is None else random_seed
            self.trace_elements = [] if trace_elements is None else trace_elements

    def __repr__(self) -> str:
        return (
        f"<MineralDataGeneration name={self.name!r}, "
        f"n_datapoints={self.n_datapoints}, "
        f"random_seed={self.random_seed}>"
        )

    def assign_mineral_group(self,
            name: Optional[str] = None,
            n_datapoints: Optional[int] = None
    ) -> None:
        if name is None:
            name = self.name
            n_datapoints = self.n_datapoints

        if name in self.list_carbonates:
            from gebpy.modules.carbonates import Carbonates
            data_mineral = Carbonates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_cyclosilicates:
            from gebpy.modules.silicates import Cyclosilicates
            data_mineral = Cyclosilicates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_halides:
            from gebpy.modules.halides import Halides
            data_mineral = Halides(
                mineral=name, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_inosilicates:
            from gebpy.modules.silicates import Inosilicates
            data_mineral = Inosilicates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_miscellaneous:
            from gebpy.modules.organics import Organics
            data_mineral = Organics(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_nesosilicates:
            from gebpy.modules.silicates import Nesosilicates
            data_mineral = Nesosilicates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_oxides:
            from gebpy.modules.oxides import Oxides
            data_mineral = Oxides(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_phosphates:
            from gebpy.modules.phosphates import Phosphates
            data_mineral = Phosphates(
                mineral=name, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_phospides:
            from gebpy.modules.phospides import Phospides
            data_mineral = Phospides(
                mineral=name, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_phyllosilicates:
            from gebpy.modules.silicates import Phyllosilicates
            data_mineral = Phyllosilicates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_sorosilicates:
            from gebpy.modules.silicates import Sorosilicates
            data_mineral = Sorosilicates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_sulfates:
            from gebpy.modules.sulfates import Sulfates
            data_mineral = Sulfates(
                mineral=name, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_sulfides:
            from gebpy.modules.sulfides import Sulfides
            data_mineral = Sulfides(
                mineral=name, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)
        elif name in self.list_tectosilicates:
            from gebpy.modules.silicates import Tectosilicates
            data_mineral = Tectosilicates(
                mineral=name, data_type=True, traces_list=self.trace_elements).generate_dataset(number=n_datapoints)

        return data_mineral

    def generate_data(self):
        data_mineral = self.assign_mineral_group(self.name, self.n_datapoints)

        return data_mineral

# DEFAULT EXAMPLE
data_config = DEFAULT_CONFIG
data_default = MineralDataGeneration(name=data_config.name, n_datapoints=data_config.n_datapoints)
DEFAULT_DATA = data_default.generate_data()