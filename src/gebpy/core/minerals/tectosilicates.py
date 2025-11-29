#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		tectosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		29.11.2025

#-----------------------------------------------

"""
Module: tectosilicates.py
This module controls the generation of the synthetic data of the tectosilicate minerals.
"""

# PACKAGES
import yaml
import numpy as np
from pathlib import Path

from asteval import Interpreter

# MODULES
from modules.chemistry import PeriodicSystem
from modules.geochemistry import MineralChemistry
from modules.geophysics import WellLog as wg
from src.gebpy.core.minerals.common import GeophysicalProperties, CrystallographicProperties, CrystalPhysics
from src.gebpy.core.minerals.common import MineralGeneration as MinGen

# CODE
BASE_PATH = Path(__file__).resolve().parents[2]
DATA_PATH = BASE_PATH / "data"

class Tectosilicates:
    _yaml_cache = {}
    _formula_cache = {}

    def __init__(self, name, random_seed, rounding=3) -> None:
        self.name = name
        self.random_seed = random_seed
        self.rng = np.random.default_rng(random_seed)
        self.current_seed = int(np.round(self.rng.uniform(0, 1000), 0))
        self.data_path = DATA_PATH
        self.rounding = rounding
        self.ae = Interpreter()
        self.cache = {}

        # Chemistry
        self.elements = {
            "H": PeriodicSystem(name="H").get_data(),
            "C": PeriodicSystem(name="C").get_data(),
            "O": PeriodicSystem(name="O").get_data(),
            "Na": PeriodicSystem(name="Na").get_data(),
            "Mg": PeriodicSystem(name="Mg").get_data(),
            "Al": PeriodicSystem(name="Al").get_data(),
            "Si": PeriodicSystem(name="Si").get_data(),
            "Cl": PeriodicSystem(name="Cl").get_data(),
            "K": PeriodicSystem(name="K").get_data(),
            "Ca": PeriodicSystem(name="Ca").get_data(),
            "Mn": PeriodicSystem(name="Mn").get_data(),
            "Fe": PeriodicSystem(name="Fe").get_data(),
            "Ni": PeriodicSystem(name="Ni").get_data(),
        }

        # Geophysics
        self.geophysical_properties = GeophysicalProperties()
        # Crystallography
        self.crystallographic_properties = CrystallographicProperties()

        # Mineral-specific data
        if self.name in [
            "Orthoclase", "Albite", "Microcline", "Anorthite", "Marialite", "Meionite", "Kalsilite", "NaNepheline"]:
            self.yaml_data = self._load_yaml(self.name.lower())
        if self.name == "Alkali feldspar":
            for mineral in ["Albite", "Orthoclase"]:
                self.yaml_data = self._load_yaml(mineral.lower())
        elif self.name == "Plagioclase":
            for mineral in ["Albite", "Anorthite"]:
                self.yaml_data = self._load_yaml(mineral.lower())
        elif self.name == "Scapolite":
            for mineral in ["Marialite", "Meionite"]: #"NaNepheline", "Kalsilite"
                self.yaml_data = self._load_yaml(mineral.lower())
        elif self.name == "Nepheline":
            for mineral in ["NaNepheline", "Kalsilite"]:
                self.yaml_data = self._load_yaml(mineral.lower())

    def _load_yaml(self, mineral_name: str) -> dict:
        # 1) Cache-Hit
        if mineral_name in Tectosilicates._yaml_cache:
            return Tectosilicates._yaml_cache[mineral_name]

        # 2) Laden von Disk
        yaml_file = self.data_path/f"{mineral_name}.yaml"
        if not yaml_file.exists():
            raise FileNotFoundError(f"No YAML file found for {mineral_name}.")

        with open(yaml_file, "r") as f:
            data = yaml.safe_load(f)

        if "chemistry" in data and mineral_name not in Tectosilicates._formula_cache:
            self._compile_chemistry_formulas(mineral_name, data["chemistry"])

        # 3) Cache schreiben
        Tectosilicates._yaml_cache[mineral_name] = data

        return data

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

    def _get_variables(self):
        if "variables" not in self.yaml_data:
            return {}
        vars = {}
        for k, v in self.yaml_data["variables"].items():
            if isinstance(v[0], int):
                vars[k] = self.rng.integers(v[0], v[1])
            else:
                vars[k] = round(self.rng.uniform(v[0], v[1]), 2)
        return vars

    def generate_dataset(self, number: int = 1, as_dataframe=False) -> None:
        fixed = [
            "Albite", "Orthoclase", "Microcline", "Anorthite", "Marialite", "Meionite", "Kalsilite", "NaNepheline"]
        variable = {}
        endmember = ["Alkali feldspar", "Plagioclase", "Scapolite", "Nepheline"]

        generators = {
            **{m: MinGen(
                name=self.name, yaml_data=self.yaml_data, elements=self.elements, cache=self.cache,
                geophysical_properties=self.geophysical_properties, rounding=self.rounding
            ).create_mineral_data_fixed_composition for m in fixed},
            **{m: self.create_mineral_data_variable_composition for m in variable},
            **{m: self.create_mineral_data_endmember_series for m in endmember},
        }

        if self.name not in generators:
            raise ValueError(f"Mineral '{self.name}' not recognized.")

        dataset = {}
        if self.name in fixed:
            dataset = self._evaluate_mineral(index=1, generators=generators, dataset=dataset)
        else:
            for index in range(number):
                dataset = self._evaluate_mineral(index=index, generators=generators, dataset=dataset)

        if as_dataframe:
            import pandas as pd
            return pd.DataFrame(dataset)
        else:
            return dataset

    def _evaluate_mineral(self, index, generators, dataset):
        self.current_seed = np.uint32(self.random_seed + index)
        data_mineral = generators[self.name]()

        for key, value in data_mineral.items():
            if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U",
                       "p"]:
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

    def _compile_chemistry_formulas(self, mineral_name: str, chemistry_dict: dict):
        """
        Extracts and compiles all chemistry formulas from the YAML file.
        Stores the compiled ASTs in the global formula cache.
        """

        compiled = {}

        for element, entry in chemistry_dict.items():

            # Fall A: Einfache Zahl (z.B. 2, 3.5)
            if isinstance(entry, (int, float)):
                expr = str(entry)
            # Fall B: dict mit 'formula'
            elif isinstance(entry, dict) and "formula" in entry:
                expr = str(entry["formula"])
            # Fall C: Alles andere ist invalid
            else:
                raise ValueError(
                    f"Invalid chemistry entry for element '{element}' in {mineral_name}.yaml: {entry}")
            # AST kompilieren
            compiled[element] = self.ae.parse(expr)
        # Cache schreiben
        Tectosilicates._formula_cache[mineral_name] = compiled

    def _evaluate_chemistry(self, chemistry_dict, **variables):
        """
        Evaluates algebraic expressions of element amounts defined in the YAML file.

        Parameters:
            chemistry_dict (dict): Dictionary from YAML containing element formulas as strings.
            **variables: Variable assignments (e.g. x=0.5, y=1, n=8)

        Returns:
            dict: {element: calculated_amount}
        """
        self.ae.symtable.clear()

        # Variablen setzen
        for k, v in variables.items():
            self.ae.symtable[k] = v

        results = {}
        compiled = Tectosilicates._formula_cache[self.name.lower()]

        for el in chemistry_dict:
            results[el] = self.ae.run(compiled[el])

        return results

    def _extract_values_from_yaml(self):
        vals = {}
        # Physical parameters
        for key in ["K", "G", "a_K", "b_K", "a_G", "b_G"]:
            if key in self.yaml_data.get("physical_properties", {}):
                vals[key] = float(self.yaml_data["physical_properties"][key]["value"])
        # Cell parameters
        for key in ["a", "b", "c", "alpha", "beta", "gamma", "Z"]:
            if key in self.yaml_data.get("cell_data", {}):
                vals[key] = float(self.yaml_data["cell_data"][key]["value"])
        # Meta data
        vals["key"] = self.yaml_data["metadata"]["key"]
        vals["crystal_system"] = self.yaml_data["metadata"]["crystal_system"]

        return vals

    def create_mineral_data_variable_composition(self):
        """
        Synthetic mineral data generation for an user-selected mineral.
        All mechanical properties (K, G, E) are stored in Pascals internally.
        For output, they are converted to GPa.
        """
        name_lower = self.name.lower()
        # Chemistry
        val_state = "variable"
        traces_data = []
        # Molar mass, elemental amounts
        majors_data = []
        molar_mass_pure = 0

        amounts_elements = self._evaluate_chemistry(self.yaml_data["chemistry"], **vars)
        for element, amount in amounts_elements.items():
            n_order = int(self.elements[element][1])
            val_amount = float(amount)
            molar_mass = float(self.elements[element][2])
            majors_data.append([element, n_order, val_amount, molar_mass])
            molar_mass_pure += val_amount*molar_mass
        majors_data.sort(key=lambda x: x[1])

        if name_lower not in self.cache:
            vals = self.cache[name_lower]["constants"]
            self.cache[name_lower] = {"constants": vals}
        else:
            vals = self.cache[name_lower]["constants"]

        # Molar mass, element amounts
        molar_mass = 0
        for element, amount in amounts_elements.items():
            molar_mass += amount*float(self.elements[element][2])

        amounts = []
        for element, amount in amounts_elements.items():
            value = amount*float(self.elements[element][2])/molar_mass
            amounts.append([element, self.elements[element][1], value])
        element = [self.elements[name] for name, *_ in amounts]

        # Reading and assigning the mineral-specific information from the YAML file
        val_key = vals["key"]
        val_system = vals["crystal_system"]
        if "K" in vals:
            val_K = vals["K"]
            val_G = vals["G"]
        else:
            val_a_K = float(vals["a_K"])
            val_b_K = float(vals["b_K"])
            val_a_G = float(vals["a_G"])
            val_b_G = float(vals["b_G"])
        val_a = vals["a"]
        val_b = vals["b"]
        val_c = vals["c"]
        val_beta = vals["beta"]
        val_Z = vals["Z"]

        if "alpha" in vals:
            val_alpha = vals["alpha"]
        if "gamma" in vals:
            val_gamma = vals["gamma"]

        # (Molar) Volume
        if "alpha" in vals and "gamma" in vals:
            constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_alpha, val_beta, val_gamma], val_system])
        else:
            constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_beta], val_system])

        constr_minchem = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data)
        V, V_m = self.crystallographic_properties.calculate_molar_volume(
            constr_volume=constr_vol, constr_molar_volume=constr_minchem, cell_z=val_Z)
        # Density
        constr_density = CrystalPhysics([molar_mass, val_Z, V])
        rho = self.crystallographic_properties.calculate_mineral_density(constr_density=constr_density)
        constr_electr_density = wg(amounts=amounts, elements=element, rho_b=rho)
        rho_e = self.crystallographic_properties.calculate_electron_density(
            constr_electron_density=constr_electr_density)

        # Elastic properties
        if "K" not in vals:
            val_K = (val_a_K*rho + val_b_K)*10**9
            val_G = (val_a_G*rho + val_b_G)*10**9
        E, nu = self.geophysical_properties.calculate_elastic_properties(bulk_mod=val_K, shear_mod=val_G)
        # Seismic properties
        vPvS, vP, vS = self.geophysical_properties.calculate_seismic_velocities(
            bulk_mod=val_K, shear_mod=val_G, rho=rho)
        # Radiation properties
        constr_radiation = wg(amounts=amounts, elements=element)
        gamma_ray, pe, U = self.geophysical_properties.calculate_radiation_properties(
            constr_radiation=constr_radiation, rho_electron=rho_e)
        # Electrical resistivity
        p = None
        # Results
        results = {
            "mineral": val_key, "state": val_state, "M": round(molar_mass, self.rounding),
            "chemistry": {name: round(val[1], 6) for name, *val in amounts}, "rho": round(rho, self.rounding),
            "rho_e": round(rho_e, self.rounding), "V": round(V_m, self.rounding), "vP": round(vP, self.rounding),
            "vS": round(vS, self.rounding), "vP/vS": round(vPvS, self.rounding),
            "K": round(val_K*10**(-9), self.rounding), "G": round(val_G*10**(-9), self.rounding),
            "E": round(E*10**(-9), self.rounding), "nu": round(nu, 6), "GR": round(gamma_ray, self.rounding),
            "PE": round(pe, self.rounding), "U": round(U, self.rounding), "p": p}
        return results

    def create_mineral_data_endmember_series(self):
        """
        Synthetic mineral data generation for an user-selected mineral.
        All mechanical properties (K, G, E) are stored in Pascals internally.
        For output, they are converted to GPa.
        """
        val_state = "variable"

        if not hasattr(self, "cache"):
            self.cache = {}

        if self.name == "Alkali feldspar":
            name_lower = self.name.lower()
            val_key = "Afs"
            endmember = ["Albite", "Orthoclase"]
        elif self.name == "Plagioclase":
            name_lower = self.name.lower()
            val_key = "Pl"
            endmember = ["Albite", "Anorthite"]
        elif self.name == "Scapolite":
            name_lower = self.name.lower()
            val_key = "Scp"
            endmember = ["Marialite", "Meionite"]
        elif self.name == "Nepheline":
            name_lower = self.name.lower()
            val_key = "Nph"
            endmember = ["NaNepheline", "Kalsilite"]

        if "endmembers" not in self.cache:
            self.cache["endmembers"] = {}

        endmember_data = {}
        list_elements = []
        for mineral in endmember:
            if mineral not in self.cache["endmembers"]:
                mineral_data = Tectosilicates(name=mineral, random_seed=self.current_seed).generate_dataset(number=1)
                self.cache["endmembers"][mineral] = mineral_data
            endmember_data[mineral] = self.cache["endmembers"][mineral]
            mineral_data = endmember_data[mineral]
            for element in mineral_data["chemistry"]:
                if element not in list_elements:
                    list_elements.append(element)
        weights = self.rng.dirichlet(np.ones(len(endmember)))
        fraction_endmember = dict(zip(endmember, weights))

        if name_lower not in self.cache:
            self.cache[name_lower] = {
                "endmember_data": endmember_data
            }

        properties = ["M", "rho", "rho_e", "V", "K", "G"]
        helper_results = {
            prop: sum(fraction_endmember[m]*endmember_data[m][prop][0] for m in endmember)
            for prop in properties
        }
        # Amounts
        amounts = []
        for element in list_elements:
            amount = sum(fraction_endmember[mineral]*endmember_data[mineral]["chemistry"].get(element, [0])[0]
                         for mineral in endmember)
            amounts.append([element, self.elements[element][1], amount])
        element = [self.elements[name] for name, *_ in amounts]
        # Elastic properties
        val_K = helper_results["K"]*10**9
        val_G = helper_results["G"]*10**9
        rho = helper_results["rho"]
        rho_e = helper_results["rho_e"]
        E, nu = self.geophysical_properties.calculate_elastic_properties(bulk_mod=val_K, shear_mod=val_G)
        # Seismic properties
        vPvS, vP, vS = self.geophysical_properties.calculate_seismic_velocities(
            bulk_mod=val_K, shear_mod=val_G, rho=rho)
        # Radiation properties
        constr_radiation = wg(amounts=amounts, elements=element)
        gamma_ray, pe, U = self.geophysical_properties.calculate_radiation_properties(
            constr_radiation=constr_radiation, rho_electron=rho_e)
        # Electrical resistivity
        p = None
        # Results
        results = {
            "mineral": val_key, "state": val_state, "M": round(helper_results["M"], self.rounding),
            "chemistry": {name: round(val[1], 6) for name, *val in amounts}, "rho": round(rho, self.rounding),
            "rho_e": round(rho_e, self.rounding), "V": round(helper_results["V"], self.rounding),
            "vP": round(vP, self.rounding), "vS": round(vS, self.rounding),
            "vP/vS": round(vPvS, self.rounding), "K": round(val_K*10**(-9), self.rounding),
            "G": round(val_G*10**(-9), self.rounding), "E": round(E*10**(-9), self.rounding), "nu": round(nu, 6),
            "GR": round(gamma_ray, self.rounding), "PE": round(pe, self.rounding), "U": round(U, self.rounding), "p": p}
        return results

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = Tectosilicates(name="Albite", random_seed=42).generate_dataset(number=10)