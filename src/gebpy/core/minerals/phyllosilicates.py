#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		phyllosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		28.11.2025

#-----------------------------------------------

"""
Module: phyllosilicates.py
This module controls the generation of the synthetic data of the phyllosilicate minerals.
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
class Phyllosilicates:
    def __init__(self, name, random_seed, rounding=3) -> None:
        self.name = name
        self.random_seed = random_seed
        self.rng = np.random.default_rng(random_seed)
        self.current_seed = int(np.round(self.rng.uniform(0, 1000), 0))
        self.data_path = Path(__file__).resolve().parents[2] / "data"
        self.rounding = rounding
        self.ae = Interpreter()
        self.cache = {}

        # Chemistry
        self.elements = {
            "H": PeriodicSystem(name="H").get_data(),
            "O": PeriodicSystem(name="O").get_data(),
            "Na": PeriodicSystem(name="Na").get_data(),
            "Mg": PeriodicSystem(name="Mg").get_data(),
            "Al": PeriodicSystem(name="Al").get_data(),
            "Si": PeriodicSystem(name="Si").get_data(),
            "K": PeriodicSystem(name="K").get_data(),
            "Ca": PeriodicSystem(name="Ca").get_data(),
            "Mn": PeriodicSystem(name="Mn").get_data(),
            "Fe": PeriodicSystem(name="Fe").get_data(),
            "Ni": PeriodicSystem(name="Ni").get_data(),
        }

        # Geophysics
        self.geophysical_properties = GeophysicalProperties()

        # Mineral-specific data
        if self.name in [
            "Annite", "Eastonite", "Illite", "Kaolinite", "Phlogopite", "Siderophyllite", "Chamosite", "Clinochlore",
            "Pennantite", "Nimite", "Muscovite", "Talc", "Chrysotile", "Antigorite", "Pyrophyllite", "Montmorillonite",
            "Nontronite", "Saponite", "Glauconite", "Vermiculite", "Chlorite"]:
            self.yaml_data = self._load_yaml(self.name.lower())
        if self.name == "Biotite":
            for mineral in ["Annite", "Phlogopite", "Siderophyllite", "Eastonite"]:
                self.yaml_data = self._load_yaml(mineral.lower())

    def _load_yaml(self, mineral_name: str) -> dict:
        yaml_file = self.data_path / f"{mineral_name}.yaml"
        if not yaml_file.exists():
            raise FileNotFoundError(f"No YAML file found for {mineral_name}.")
        with open(yaml_file, "r") as f:
            return yaml.safe_load(f)

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
        if self.name == "Glauconite":
            x = round(self.rng.uniform(0, 1), 2)
            y1 = round(self.rng.uniform(0, 1), 2)
            y2 = round(self.rng.uniform(0, 1 - y1), 2)
            z = round(self.rng.uniform(0, 1), 2)
            return {"x": x, "y1": y1, "y2": y2, "z": z}
        elif self.name == "Vermiculite":
            x = round(self.rng.uniform(0, 1), 2)
            y = round(self.rng.uniform(0, (1 - x)), 2)
            z = round(self.rng.uniform(0, 1), 2)
            return {"x": x, "y": y, "z": z}
        elif self.name == "Chlorite":
            x = round(self.rng.uniform(0, 1), 2)
            y = round(self.rng.uniform(0, (1 - x)), 2)
            upper = max(0, (1 - x - y))
            z = round(self.rng.uniform(0, upper), 2)
            return {"x": x, "y": y, "z": z}
        else:
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
        fixed = {
            "Annite", "Eastonite", "Illite", "Kaolinite", "Phlogopite", "Siderophyllite",
            "Chamosite", "Clinochlore", "Pennantite", "Nimite", "Muscovite",
            "Talc", "Chrysotile", "Antigorite", "Pyrophyllite"}
        variable = {"Montmorillonite", "Nontronite", "Saponite", "Glauconite", "Vermiculite", "Chlorite"}
        endmember = {"Biotite"}

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
        # Übergib alle Variablen an die Symboltabelle
        for k, v in variables.items():
            self.ae.symtable[k] = v

        # Berechne für jedes Element den Ausdruck aus der YAML
        results = {}
        for el, data in chemistry_dict.items():
            expr = str(data["formula"])
            try:
                results[el] = self.ae(expr)
            except Exception as e:
                raise ValueError(f"Error evaluating formula for element '{el}': {expr} ({e})")

        return results

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
        vars = self._get_variables()

        amounts_elements = self._evaluate_chemistry(self.yaml_data["chemistry"], **vars)
        for element, amount in amounts_elements.items():
            n_order = int(self.elements[element][1])
            val_amount = float(amount)
            molar_mass = float(self.elements[element][2])
            majors_data.append([element, n_order, val_amount, molar_mass])
            molar_mass_pure += val_amount*molar_mass
        majors_data.sort(key=lambda x: x[1])

        if name_lower not in self.cache:
            vals = {}
            for key in ["K", "G", "a_K", "b_K", "a_G", "b_G", "a", "b", "c", "alpha", "beta", "gamma", "Z"]:
                if key in self.yaml_data["cell_data"] or key in self.yaml_data["physical_properties"]:
                    vals[key] = self._get_value(self.yaml_data, ["physical_properties", key]) \
                                if key in ["K", "G", "a_K", "b_K", "a_G", "b_G"] else \
                                self._get_value(self.yaml_data, ["cell_data", key])
            for key in ["key", "crystal_system"]:
                vals[key] = self._get_value(self.yaml_data, ["metadata", key])

            self.cache[name_lower] = {
                "constants": vals,
            }
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

        if self.name != "Chlorite":
            # (Molar) Volume
            if "alpha" in vals and "gamma" in vals:
                constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_alpha, val_beta, val_gamma], val_system])
            else:
                constr_vol = CrystalPhysics([[val_a, val_b, val_c], [val_beta], val_system])

            constr_minchem = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data)
            V, V_m = CrystallographicProperties().calculate_molar_volume(
                constr_volume=constr_vol, constr_molar_volume=constr_minchem, cell_z=val_Z)
            # Density
            constr_density = CrystalPhysics([molar_mass, val_Z, V])
            rho = CrystallographicProperties().calculate_mineral_density(constr_density=constr_density)
            constr_electr_density = wg(amounts=amounts, elements=element, rho_b=rho)
            rho_e = CrystallographicProperties().calculate_electron_density(
                constr_electron_density=constr_electr_density)
        else:
            # Density
            dataV_Fe = CrystalPhysics([[5.373, 9.306, 14.222], [97.88], "monoclinic"])
            V_Fe = dataV_Fe.calculate_volume()
            Z_Fe = 2
            V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
            dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
            rho_Fe = dataRho_Fe.calculate_bulk_density()
            rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()

            dataV_Mg = CrystalPhysics([[5.3, 9.3, 14.3], [97], "monoclinic"])
            V_Mg = dataV_Mg.calculate_volume()
            Z_Mg = 2
            V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
            dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
            rho_Mg = dataRho_Mg.calculate_bulk_density()
            rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()

            dataV_Mn = CrystalPhysics([[5.454, 9.45, 14.4], [97.2], "monoclinic"])
            V_Mn = dataV_Mn.calculate_volume()
            Z_Mn = 2
            V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
            dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
            rho_Mn = dataRho_Mn.calculate_bulk_density()
            rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()

            dataV_Ni = CrystalPhysics([[5.32, 9.214, 14.302], [97.1], "monoclinic"])
            V_Ni = dataV_Ni.calculate_volume()
            Z_Ni = 2
            V_m_Ni = MineralChemistry().calculate_molar_volume(volume_cell=V_Ni, z=Z_Ni)
            dataRho_Ni = CrystalPhysics([molar_mass, Z_Ni, V_Ni])
            rho_Ni = dataRho_Ni.calculate_bulk_density()
            rho_e_Ni = wg(amounts=amounts, elements=element, rho_b=rho_Ni).calculate_electron_density()

            x = vars["x"]
            y = vars["y"]
            z = vars["z"]
            V_m = x*V_m_Fe + y*V_m_Mg + z*V_m_Mn + (1-x-y-z)*V_m_Ni
            rho = x*rho_Fe + y*rho_Mg + z*rho_Mn + (1-x-y-z)*rho_Ni
            rho_e = x*rho_e_Fe + y*rho_e_Mg + z*rho_e_Mn + (1-x-y-z)*rho_e_Ni

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

        if self.name == "Biotite":
            name_lower = self.name.lower()
            val_key = "Bt"
            endmember = ["Annite", "Phlogopite", "Siderophyllite", "Eastonite"]

        if "endmembers" not in self.cache:
            self.cache["endmembers"] = {}

        endmember_data = {}
        list_elements = []
        for mineral in endmember:
            if mineral not in self.cache["endmembers"]:
                mineral_data = Phyllosilicates(name=mineral, random_seed=self.current_seed).generate_dataset(number=1)
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
            "vP": round(vP, self.rounding), "vS": round(vS, self.rounding), "vP/vS": round(vPvS, self.rounding),
            "K": round(val_K*10**(-9), self.rounding), "G": round(val_G*10**(-9), self.rounding),
            "E": round(E*10**(-9), self.rounding), "nu": round(nu, 6), "GR": round(gamma_ray, self.rounding), #
            "PE": round(pe, self.rounding), "U": round(U, self.rounding), "p": p}
        return results

# TEST
if __name__ == "__main__":
    DEFAULT_DATA = Phyllosilicates(name="Annite", random_seed=42).generate_dataset(number=10)