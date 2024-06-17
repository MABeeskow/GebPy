#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		silicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.06.2024

# -----------------------------------------------

# MODULES
import re
import numpy as np
import random as rd
from scipy import stats
from modules.chemistry import PeriodicSystem, DataProcessing, OxideCompounds
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry, TraceElements

# TECTOSILICATES
class Tectosilicates:
    """ Class that generates geophysical and geochemical data of tectosiliacte minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral

        ## Chemistry
        boron = PeriodicSystem(name="B").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        chlorine = PeriodicSystem(name="Cl").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()

        self.boron = ["B", 5, boron[2]]
        self.carbon = ["C", 6, carbon[2]]
        self.oxygen = ["O", 8, oxygen[2]]
        self.sodium = ["Na", 11, sodium[2]]
        self.aluminium = ["Al", 13, aluminium[2]]
        self.silicon = ["Si", 14, silicon[2]]
        self.chlorine = ["Cl", 17, chlorine[2]]
        self.potassium = ["K", 19, potassium[2]]
        self.calcium = ["Ca", 20, calcium[2]]
    #
    def get_data(self, number=1):
        if self.mineral in ["Afs", "Kfs", "Alkali Feldspar"]:
            self.molar_mass_base = round(self.aluminium[2] + 3*self.silicon[2] + 8*self.oxygen[2], 3)
            if number > 1:
                data = [self.create_alkalifeldspar() for n in range(number)]
            else:
                data = self.create_alkalifeldspar()
        elif self.mineral in ["Pl", "Plagioclase"]:
            if number > 1:
                data = [self.create_plagioclase() for n in range(number)]
            else:
                data = self.create_plagioclase()
        elif self.mineral in ["Scp", "Scapolite"]:
            if number > 1:
                data = [self.create_scapolite() for n in range(number)]
            else:
                data = self.create_scapolite()
        elif self.mineral in ["Dnb", "Danburite"]:
            if number > 1:
                data = [self.create_danburite() for n in range(number)]
            else:
                data = self.create_danburite()
        elif self.mineral in ["Nph", "Nepheline"]:
            if number > 1:
                data = [self.create_nepheline() for n in range(number)]
            else:
                data = self.create_nepheline()
        elif self.mineral in ["Fsp", "Orthoclase"]:
            if number > 1:
                data = [self.create_orthoclase() for n in range(number)]
            else:
                data = self.create_orthoclase()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}

        for index in range(number):
            if self.mineral in ["Alkali Feldspar", "Alkaline Feldspar"]:
                data_mineral = self.create_alkalifeldspar()
            elif self.mineral == "Plagioclase":
                data_mineral = self.create_plagioclase()
            elif self.mineral == "Scapolite":
                data_mineral = self.create_scapolite()
            elif self.mineral == "Danburite":
                data_mineral = self.create_danburite()
            elif self.mineral == "Nepheline":
                data_mineral = self.create_nepheline()
            elif self.mineral == "Orthoclase":
                data_mineral = self.create_orthoclase()
            #
            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U",
                           "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif (key in ["mineral", "state", "trace elements", "major elements", "major compounds", "LA-ICP-MS"] and key
                      not in dataset):
                    dataset[key] = value
                elif key in ["chemistry", "compounds"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)
        #
        return dataset
    #
    def create_alkalifeldspar(self, enrichment=None, x_value=None):
        self.enrichment = enrichment
        self.x_value = x_value

        # Major Elements
        majors_name = ["O", "Na", "Al", "Si", "K"]
        # Trace elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["K", "Mg", "Ti", "Fe"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        ## Molar mass
        condition = False
        while condition == False:
            if self.x_value == None:
                if self.enrichment == None:
                    x = round(rd.uniform(0, 1), 4)
                elif self.enrichment == "Na":
                    x = round(rd.uniform(0.75, 1), 4)
                elif self.enrichment == "K":
                    x = round(rd.uniform(0, 0.25), 4)
            else:
                x = self.x_value

            mineral = "Kfs"
            majors_data = np.array(
                [["O", self.oxygen[1], 8, self.oxygen[2]],
                 ["Na", self.sodium[1], round(x, 4), self.sodium[2]],
                 ["Al", self.aluminium[1], 1, self.aluminium[2]],
                 ["Si", self.silicon[1], 3, self.silicon[2]],
                 ["K", self.potassium[1], round((1-x), 4), self.potassium[2]]], dtype=object)

            molar_mass_pure = round(x*self.sodium[2] + (1-x)*self.potassium[2] + self.aluminium[2] + 3*self.silicon[2]
                                    + 8*self.oxygen[2], 3)
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            ## Oxide Composition
            list_oxides = ["Na2O", "K2O", "Al2O3", "SiO2"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)

            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass

        ## Density
        M_Ab = round(self.sodium[2] + self.aluminium[2] + 3*self.silicon[2] + 8*self.oxygen[2], 3)
        dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
        V_Ab = dataV_Ab.calculate_volume()
        V_m_Ab = MineralChemistry().calculate_molar_volume(volume_cell=V_Ab, z=4)
        dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
        rho_Ab = dataRho_Ab.calculate_bulk_density()
        rho_e_Ab = wg(amounts=amounts, elements=element, rho_b=rho_Ab).calculate_electron_density()

        M_Or = round(self.potassium[2] + self.aluminium[2] + 3*self.silicon[2] + 8*self.oxygen[2], 3)
        dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
        V_Or = dataV_Or.calculate_volume()
        V_m_Or = MineralChemistry().calculate_molar_volume(volume_cell=V_Or, z=4)
        dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
        rho_Or = dataRho_Or.calculate_bulk_density()
        rho_e_Or = wg(amounts=amounts, elements=element, rho_b=rho_Or).calculate_electron_density()
        rho = x*rho_Ab + (1 - x)*rho_Or
        rho_e = x*rho_e_Ab + (1 - x)*rho_e_Or
        V = x*V_m_Ab + (1 - x)*V_m_Or
        ## Bulk modulus
        K_Ab = 103.45*10**9
        K_Or = 89.78*10**9
        K = (x*K_Ab + (1-x)*K_Or)
        ## Shear modulus
        G_Ab = 69.0*10**9
        G_Or = 61.52*10**9
        G = (x*G_Ab + (1-x)*G_Or)
        ## Young's modulus
        E = (9*K*G)/(3*K + G)
        ## Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        ## vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        ## P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        ## S-wave velocity
        vS = (G/rho)**0.5
        ## Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        ## Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        ## Electrical resistivity
        p = 5.5*10**11
        ## LA-ICP-MS
        normalized_sensitivity_Si = 44.54
        dict_element_oxides = {}
        for oxide in composition_oxides:
            key = key = re.search("(\D+)(\d*)(\D+)(\d*)", oxide)
            dict_element_oxides[key.group(1)] = oxide
        # RESULTS
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si, "Oxides": dict_element_oxides}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_plagioclase(self, enrichment=None, x_value=None):
        self.enrichment = enrichment
        self.x_value = x_value
        #
        # Major Elements
        majors_name = ["O", "Na", "Al", "Si", "Ca"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["K", "Mg", "Ti", "Fe"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        condition = False
        while condition == False:
            if self.x_value == None:
                if self.enrichment == None:
                    x = round(rd.uniform(0, 1), 4)
                elif self.enrichment == "Na":
                    x = round(rd.uniform(0.75, 1), 4)
                elif self.enrichment == "Ca":
                    x = round(rd.uniform(0, 0.25), 4)
            else:
                x = self.x_value

            mineral = "Pl"
            majors_data = np.array(
                [["O", self.oxygen[1], 8, self.oxygen[2]],
                 ["Na", self.sodium[1], round(x, 4), self.sodium[2]],
                 ["Al", self.aluminium[1], round((2-x), 4), self.aluminium[2]],
                 ["Si", self.silicon[1], round((2+x), 4), self.silicon[2]],
                 ["Ca", self.calcium[1], round((1-x), 4), self.calcium[2]]], dtype=object)
            #
            molar_mass_pure = round(x*self.sodium[2] + (1-x)*self.calcium[2] + (2-x)*self.aluminium[2]
                                    + (2+x)*self.silicon[2] + 8*self.oxygen[2], 3)
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            ## Oxide Composition
            list_oxides = ["Na2O", "CaO", "Al2O3", "SiO2"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)

            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass

        # Density
        M_Ab = round(self.sodium[2] + self.aluminium[2] + 3*self.silicon[2] + 8*self.oxygen[2], 3)
        dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
        V_Ab = dataV_Ab.calculate_volume()
        V_m_Ab = MineralChemistry().calculate_molar_volume(volume_cell=V_Ab, z=4)
        dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
        rho_Ab = dataRho_Ab.calculate_bulk_density()
        rho_e_Ab = wg(amounts=amounts, elements=element, rho_b=rho_Ab).calculate_electron_density()
        #
        M_An = round(self.calcium[2] + 2*self.aluminium[2] + 2*self.silicon[2] + 8*self.oxygen[2], 3)
        dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
        V_An = dataV_An.calculate_volume()
        V_m_An = MineralChemistry().calculate_molar_volume(volume_cell=V_An, z=8)
        dataRho_An = CrystalPhysics([M_An, 8, V_An])
        rho_An = dataRho_An.calculate_bulk_density()
        rho_e_An = wg(amounts=amounts, elements=element, rho_b=rho_An).calculate_electron_density()
        rho = x*rho_Ab + (1-x)*rho_An
        rho_e = x*rho_e_Ab + (1-x)*rho_e_An
        V = x*V_m_Ab + (1-x)*V_m_An
        # Bulk modulus
        K_Ab = 103.45*10**9
        K_An = 114.72*10**9
        K = (x*K_Ab + (1-x)*K_An)
        # Shear modulus
        G_Ab = 69.0*10**9
        G_An = 73.72*10**9
        G = (x*G_Ab + (1-x)*G_An)
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = 5.5*10**11
        # LA-ICP-MS
        normalized_sensitivity_Si = 44.54
        dict_element_oxides = {}
        for oxide in composition_oxides:
            key = key = re.search("(\D+)(\d*)(\D+)(\d*)", oxide)
            dict_element_oxides[key.group(1)] = oxide
        # RESULTS
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si, "Oxides": dict_element_oxides}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_scapolite(self, enrichment=None):
        name = "Scp"
        ## Major Elements
        majors_name = ["C", "O", "Na", "Al", "Si", "Cl", "Ca"]
        ## Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "K", "S", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]

        ## Molar mass
        condition = False
        while condition == False:
            if enrichment in ["Scapolite", "scapolite", "Scp", "scp", None]:
                x = round(rd.uniform(0, 1.0), 4)
            elif enrichment in ["Ca", "Meionite", "meionite", "Mei", "mei"]:
                x = round(rd.uniform(0.75, 1.0), 4)
            elif enrichment in ["Na", "Marialite", "marialite", "Mar", "mar"]:
                x = round(rd.uniform(0, 0.25), 4)

            majors_data = np.array(
                [["C", self.carbon[1], round(x, 4), self.carbon[2]],
                 ["O", self.oxygen[1], round(24 + 3*x, 4), self.oxygen[2]],
                 ["Na", self.sodium[1], round(4*(1 - x), 4), self.sodium[2]],
                 ["Al", self.aluminium[1], round(3*(1 + x), 4), self.aluminium[2]],
                 ["Si", self.silicon[1], round(3*(3 - x), 4), self.silicon[2]],
                 ["Cl", self.chlorine[1], round((1 - x), 4), self.chlorine[2]],
                 ["Ca", self.calcium[1], round(4*x, 4), self.calcium[2]]], dtype=object)

            molar_mass_pure = round(4*(x*self.calcium[2] + (1 - x)*self.sodium[2]) + 3*(1 + x)*self.aluminium[2] \
                              + 3*(3 - x)*self.silicon[2] + 24*self.oxygen[2] \
                              + x*(self.carbon[2] + 3*self.oxygen[2]) + (1 - x)*self.chlorine[2], 3)
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

            majors_data_Mar = np.array(
                [["C", self.carbon[1], 0, self.carbon[2]],
                 ["O", self.oxygen[1], 24, self.oxygen[2]],
                 ["Na", self.sodium[1], 4, self.sodium[2]],
                 ["Al", self.aluminium[1], 3, self.aluminium[2]],
                 ["Si", self.silicon[1], 9, self.silicon[2]],
                 ["Cl", self.chlorine[1], 1, self.chlorine[2]],
                 ["Ca", self.calcium[1], 0, self.calcium[2]]], dtype=object)

            molar_mass_pure_Mar = round(4*self.sodium[2] + 3*self.aluminium[2] + 9*self.silicon[2] + 24*self.oxygen[2]
                                        + self.chlorine[2], 3)
            molar_mass_Mar, amounts_Mar = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Mar, majors=majors_data_Mar).calculate_molar_mass()
            element_Mar = [PeriodicSystem(name=amounts_Mar[i][0]).get_data() for i in range(len(amounts_Mar))]

            majors_data_Mei = np.array(
                [["C", self.carbon[1], 1, self.carbon[2]],
                 ["O", self.oxygen[1], round(24 + 3, 4), self.oxygen[2]],
                 ["Na", self.sodium[1], 0, self.sodium[2]],
                 ["Al", self.aluminium[1], 6, self.aluminium[2]],
                 ["Si", self.silicon[1], 6, self.silicon[2]],
                 ["Cl", self.chlorine[1], 0, self.chlorine[2]],
                 ["Ca", self.calcium[1], 4, self.calcium[2]]], dtype=object)

            molar_mass_pure_Mei = round(4*self.calcium[2] + 6*self.aluminium[2] + 6*self.silicon[2] + 24*self.oxygen[2]
                                        + self.carbon[2] + 3*self.oxygen[2], 3)
            molar_mass_Mei, amounts_Mei = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Mei, majors=majors_data_Mei).calculate_molar_mass()
            element_Mei = [PeriodicSystem(name=amounts_Mei[i][0]).get_data() for i in range(len(amounts_Mei))]
            ## Oxide Composition
            list_oxides = ["CO2", "Na2O", "Al2O3", "SiO2", "CaO"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)

            w_total = np.sum(list(composition_oxides.values()))
            composition_oxides["Cl"] = round(1 - w_total, 6)
            list_oxides = ["CO2", "Na2O", "Al2O3", "SiO2", "Cl", "CaO"]
            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                if enrichment in ["Scapolite", "scapolite", "Scp", "scp", None]:
                    x = round(rd.uniform(0, 1.0), 4)
                elif enrichment in ["Ca", "Meionite", "meionite", "Mei", "mei"]:
                    x = round(rd.uniform(0.75, 1.0), 4)
                elif enrichment in ["Na", "Marialite", "marialite", "Mar", "mar"]:
                    x = round(rd.uniform(0, 0.25), 4)
        # Density
        M_Mar = round(4*self.sodium[2] + 3*self.aluminium[2] + 9*self.silicon[2] + 24*self.oxygen[2] +
                      self.chlorine[2], 3)
        dataV_Mar = CrystalPhysics([[12.075, 7.516], [], "tetragonal"])
        V_Mar = dataV_Mar.calculate_volume()
        V_m_Mar = MineralChemistry().calculate_molar_volume(volume_cell=V_Mar, z=2)
        dataRho_Mar = CrystalPhysics([M_Mar, 2, V_Mar])
        rho_Mar = dataRho_Mar.calculate_bulk_density()
        rho_e_Mar = wg(amounts=amounts_Mar, elements=element_Mar, rho_b=rho_Mar).calculate_electron_density()
        #
        M_Mei = round(4*self.calcium[2] + 6*self.aluminium[2] + 6*self.silicon[2] + 24*self.oxygen[2] + self.carbon[2] +
                      3*self.oxygen[2], 3)
        dataV_Mei = CrystalPhysics([[12.26, 7.61], [], "tetragonal"])
        V_Mei = dataV_Mei.calculate_volume()
        V_m_Mei = MineralChemistry().calculate_molar_volume(volume_cell=V_Mei, z=2)
        dataRho_Mei = CrystalPhysics([M_Mei, 2, V_Mei])
        rho_Mei = dataRho_Mei.calculate_bulk_density()
        rho_e_Mei = wg(amounts=amounts_Mei, elements=element_Mei, rho_b=rho_Mei).calculate_electron_density()
        #
        V = x*V_m_Mar + (1-x)*V_m_Mei
        rho = x*rho_Mar + (1-x)*rho_Mei
        rho_e = x*rho_e_Mar + (1-x)*rho_e_Mei
        # Bulk modulus
        K_Mar = 94.58*10**9
        K_Mei = rd.uniform(0.95, 1.05)*94.58*10**9
        K = x*K_Mar + (1-x)*K_Mei
        # Shear modulus
        G_Mar = 62.83*10**9
        G_Mei = rd.uniform(0.95, 1.05)*62.83*10**9
        G = x*G_Mar + (1-x)*G_Mei
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        # LA-ICP-MS
        normalized_sensitivity_Si = 44.54
        dict_element_oxides = {}
        for oxide in composition_oxides:
            if oxide != "Cl":
                key = key = re.search("(\D+)(\d*)(\D+)(\d*)", oxide)
                dict_element_oxides[key.group(1)] = oxide
            else:
                dict_element_oxides[oxide] = oxide
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si, "Oxides": dict_element_oxides}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_danburite(self):    # Ca B2 (SiO4)2
        name = "Dnb"
        # Major Elements
        majors_name = ["B", "O", "Si", "Ca"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Al", "Mg", "Sr", "Na"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]

        # Molar mass
        condition = False
        while condition == False:
            majors_data = np.array(
                [["B", self.boron[1], 2, self.boron[2]], ["O", self.oxygen[1], 8, self.oxygen[2]],
                 ["Si", self.silicon[1], 2, self.silicon[2]],
                 ["Ca", self.calcium[1], 1, self.calcium[2]]], dtype=object)
            #
            molar_mass_pure = round(self.calcium[2] + 2*self.boron[2] + 2*(self.silicon[2] + 4*self.oxygen[2]), 3)
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            ## Oxide Composition
            list_oxides = ["B2O3", "CaO", "SiO2"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)

            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass
        # Density
        dataV = CrystalPhysics([[8.048, 8.763, 7.731], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 145.62*10**9
        # Shear modulus
        G = 102.17*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        # LA-ICP-MS
        normalized_sensitivity_Si = 44.54
        dict_element_oxides = {}
        for oxide in composition_oxides:
            key = key = re.search("(\D+)(\d*)(\D+)(\d*)", oxide)
            dict_element_oxides[key.group(1)] = oxide
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si, "Oxides": dict_element_oxides}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results

    def create_orthoclase(self):    # K Al Si3 O8
        name = "Kfs"
        # Major Elements
        majors_name = ["O", "Al", "Si", "K"]
        major_compounds = ["Al2O3", "SiO2", "K2O"]
        # Trace elements
        element_traces = {
            "2+": ["Ba", "Ca", "Fe"],
            "1+": ["Na", "Rb"],
            "All": ["Fe", "Ba", "Ca", "Na", "Rb"]}

        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity == "pure"
            var_state = "fixed"

        # Molar mass
        molar_mass_pure = round(self.potassium[2] + self.aluminium[2] + 3*self.silicon[2] + 8*self.oxygen[2], 3)
        condition = False
        while condition == False:
            molar_mass = 0
            amounts = []
            compositon_data = TraceElements(tracer=self.traces_list).calculate_composition_orthoclase()
            for element in compositon_data:
                chem_data = PeriodicSystem(name=element).get_data()
                molar_mass += compositon_data[element]["x"]*chem_data[2]
                amounts.append([chem_data[0], chem_data[1], compositon_data[element]["w"]])

            magic_factor = molar_mass/molar_mass_pure
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            ## Oxide Composition
            list_oxides = ["Al2O3", "SiO2", "K2O"]
            if len(self.traces_list) > 0:
                for key in self.traces_list.keys():
                    if key in ["Ba", "Ca", "Fe"]:
                        list_oxides.append(str(key) + "O")
                    elif key in ["Na", "Rb"]:
                        list_oxides.append(str(key) + "2O")
            composition_oxides = {}
            w_oxide_total = 0
            for index_oxide, var_oxide in enumerate(list_oxides):
                if index_oxide < len(list_oxides) - 1:
                    oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                    if self.impurity == "pure":
                        composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                    else:
                        composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)
                    w_oxide_total += composition_oxides[var_oxide]
                else:
                    if self.impurity == "pure":
                        composition_oxides[var_oxide] = round(1 - w_oxide_total, 4)
                    else:
                        composition_oxides[var_oxide] = round(1 - w_oxide_total, 6)
            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass
        # Density
        dataV = CrystalPhysics([[8.625, 12.996, 7.193], [116.016], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)*magic_factor
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()*magic_factor
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 89.78*10**9*magic_factor
        # Shear modulus
        G = 61.52*10**9*magic_factor
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        # LA-ICP-MS
        normalized_sensitivity_Si = 44.54
        dict_element_oxides = {}
        for oxide in composition_oxides:
            key = key = re.search("(\D+)(\d*)(\D+)(\d*)", oxide)
            dict_element_oxides[key.group(1)] = oxide
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["trace elements"] = element_traces
        results["major elements"] = majors_name
        results["major compounds"] = major_compounds
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si, "Oxides": dict_element_oxides}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_nepheline(self, x_value=None):  # (Na,K) Al SiO4
        name = "Nph"
        # Major Elements
        majors_name = ["O", "Na", "Al", "Si", "K"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mg", "Ca", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]

        # Molar mass
        condition = False
        while condition == False:
            if x_value == None:
                x = round(rd.uniform(0.7, 0.8), 4)
            else:
                x = x_value

            majors_data = np.array([
                ["O", self.oxygen[1], 4, self.oxygen[2]], ["Na", self.sodium[1], x, self.sodium[2]],
                ["Al", self.aluminium[1], 1, self.aluminium[2]], ["Si", self.silicon[1], 1, self.silicon[2]],
                ["K", self.potassium[1], (1 - x), self.potassium[2]]], dtype=object)
            majors_data_Na = np.array([
                ["O", self.oxygen[1], 4, self.oxygen[2]], ["Na", self.sodium[1], 1, self.sodium[2]],
                ["Al", self.aluminium[1], 1, self.aluminium[2]], ["Si", self.silicon[1], 1, self.silicon[2]],
                ["K", self.potassium[1], 0, self.potassium[2]]], dtype=object)
            majors_data_K = np.array([
                ["O", self.oxygen[1], 4, self.oxygen[2]], ["Na", self.sodium[1], 0, self.sodium[2]],
                ["Al", self.aluminium[1], 1, self.aluminium[2]], ["Si", self.silicon[1], 1, self.silicon[2]],
                ["K", self.potassium[1], 1, self.potassium[2]]], dtype=object)

            molar_mass_pure = (x*self.sodium[2] + (1 - x)*self.potassium[2] + self.aluminium[2] + self.silicon[2] +
                               4*self.oxygen[2])
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

            molar_mass_pure_Na = self.sodium[2] + self.aluminium[2] + self.silicon[2] + 4*self.oxygen[2]
            molar_mass_Na, amounts_Na = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Na, majors=majors_data_Na).calculate_molar_mass()
            element_Na = [PeriodicSystem(name=amounts_Na[i][0]).get_data() for i in range(len(amounts_Na))]

            molar_mass_pure_K = self.potassium[2] + self.aluminium[2] + self.silicon[2] + 4*self.oxygen[2]
            molar_mass_K, amounts_K = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_K, majors=majors_data_K).calculate_molar_mass()
            element_K = [PeriodicSystem(name=amounts_K[i][0]).get_data() for i in range(len(amounts_K))]
            ## Oxide Composition
            list_oxides = ["Na2O", "Al2O3", "K2O", "SiO2"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)

            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass
        # Density
        dataV_Na = CrystalPhysics([[10.01, 8.405], [], "hexagonal"])
        V_Na = dataV_Na.calculate_volume()
        Z_Na = 8
        V_m_Na = MineralChemistry().calculate_molar_volume(volume_cell=V_Na, z=Z_Na)
        dataRho_Na = CrystalPhysics([molar_mass_Na, Z_Na, V_Na*10**(6)])
        rho_Na = dataRho_Na.calculate_bulk_density()
        rho_e_Na = wg(amounts=amounts_Na, elements=element_Na, rho_b=rho_Na).calculate_electron_density()

        dataV_K = CrystalPhysics([[5.157, 8.706], [], "trigonal"])
        V_K = dataV_K.calculate_volume()
        Z_K = 2
        V_m_K = MineralChemistry().calculate_molar_volume(volume_cell=V_K, z=Z_K)
        dataRho_K = CrystalPhysics([molar_mass_K, Z_K, V_K])
        rho_K = dataRho_K.calculate_bulk_density()
        rho_e_K = wg(amounts=amounts_K, elements=element_K, rho_b=rho_K).calculate_electron_density()

        x_vector = np.array([x, (1 - x)])

        V_m_vector = np.array([V_m_Na, V_m_K])
        rho_vector = np.array([rho_Na, rho_K])
        rho_e_vector = np.array([rho_e_Na, rho_e_K])
        V_m = np.dot(x_vector, V_m_vector)
        rho = np.dot(x_vector, rho_vector)
        rho_e = np.dot(x_vector, rho_e_vector)
        # Bulk modulus
        K_Na = 92.27*10**9
        K_K = 78.87*10**9
        K_vector = np.array([K_Na, K_K])
        K = np.dot(x_vector, K_vector)
        # Shear modulus
        G_Na = 64.72*10**9
        G_K = 55.98*10**9
        G_vector = np.array([G_Na, G_K])
        G = np.dot(x_vector, G_vector)
        # Other elastic moduli
        denominator_base = np.array([1/(3*K + G), 1/(3*K + G)])
        youngs_poisson_vector = np.array([(9*K*G), (3*K - 2*G)/2])
        youngs_poisson_results = denominator_base*youngs_poisson_vector
        # Young's modulus
        E = youngs_poisson_results[0]
        # Poisson's ratio
        nu = youngs_poisson_results[1]
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe * rho_e * 10 ** (-3)
        # Electrical resistivity
        p = None
        # LA-ICP-MS
        normalized_sensitivity_Si = 44.54
        dict_element_oxides = {}
        for oxide in composition_oxides:
            key = key = re.search("(\D+)(\d*)(\D+)(\d*)", oxide)
            dict_element_oxides[key.group(1)] = oxide
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G * 10 ** (-9), 4)
        results["K"] = round(K * 10 ** (-9), 4)
        results["E"] = round(E * 10 ** (-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si, "Oxides": dict_element_oxides}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
#
class Phyllosilicates:
    """ Class that generates geophysical and geochemical data of phyllosilicate minerals"""

    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral

        self.hydrogen = ["H", 1, 1.008]
        self.carbon = ["C", 6, 12.011]
        self.oxygen = ["O", 8, 15.999]
        self.flourine = ["F", 9, 18.998]
        self.sodium = ["Na", 11, 22.990]
        self.magnesium = ["Mg", 12, 24.304]
        self.aluminium = ["Al", 13, 26.982]
        self.silicon = ["Si", 14, 28.085]
        self.chlorine = ["Cl", 17, 35.450]
        self.potassium = ["K", 19, 39.098]
        self.calcium = ["Ca", 20, 40.078]
        self.iron = ["Fe", 26, 55.845]

    def get_data(self, number=1):
        if self.mineral in ["Ilt", "Illite"]:
            if number > 1:
                data = [self.create_illite() for n in range(number)]
            else:
                data = self.create_illite()
        elif self.mineral in ["Kln", "Kaolinite"]:
            if number > 1:
                data = [self.create_kaolinite() for n in range(number)]
            else:
                data = self.create_kaolinite()
        elif self.mineral in ["Mnt", "Montmorillonite"]:
            if number > 1:
                data = [self.create_montmorillonite() for n in range(number)]
            else:
                data = self.create_montmorillonite()
        elif self.mineral in ["Chl", "Chamosite"]:
            if number > 1:
                data = [self.create_chamosite() for n in range(number)]
            else:
                data = self.create_chamosite()
        elif self.mineral in ["Chl", "Clinochlore"]:
            if number > 1:
                data = [self.create_clinochlore() for n in range(number)]
            else:
                data = self.create_clinochlore()
        elif self.mineral in ["Chl", "Chamosite"]:
            if number > 1:
                data = [self.create_chamosite() for n in range(number)]
            else:
                data = self.create_chamosite()
        elif self.mineral in ["Chl", "Pennantite"]:
            if number > 1:
                data = [self.create_pennantite() for n in range(number)]
            else:
                data = self.create_pennantite()
        elif self.mineral in ["Chl", "Chamosite"]:
            if number > 1:
                data = [self.create_chamosite() for n in range(number)]
            else:
                data = self.create_chamosite()
        elif self.mineral in ["Chl", "Nimite"]:
            if number > 1:
                data = [self.create_nimite() for n in range(number)]
            else:
                data = self.create_nimite()
        elif self.mineral in ["Chl", "Chlorite"]:
            if number > 1:
                data = [self.create_chlorite() for n in range(number)]
            else:
                data = self.create_chlorite()
        elif self.mineral in ["Vrm", "Vermiculute"]:
            if number > 1:
                data = [self.create_vermiculite() for n in range(number)]
            else:
                data = self.create_vermiculite()
        elif self.mineral in ["Ann", "Annite"]:
            if number > 1:
                data = [self.create_annite() for n in range(number)]
            else:
                data = self.create_annite()
        elif self.mineral in ["Phl", "Phlogopite"]:
            if number > 1:
                data = [self.create_phlogopite() for n in range(number)]
            else:
                data = self.create_phlogopite()
        elif self.mineral in ["Eas", "Eastonite"]:
            if number > 1:
                data = [self.create_eastonite() for n in range(number)]
            else:
                data = self.create_eastonite()
        elif self.mineral in ["Sdp", "Siderophyllite"]:
            if number > 1:
                data = [self.create_siderophyllite() for n in range(number)]
            else:
                data = self.create_siderophyllite()
        elif self.mineral in ["Bt", "Biotite"]:
            if number > 1:
                data = [self.create_biotite() for n in range(number)]
            else:
                data = self.create_biotite()
        elif self.mineral in ["Ms", "Muscovite"]:
            if number > 1:
                data = [self.create_muscovite() for n in range(number)]
            else:
                data = self.create_muscovite()
        elif self.mineral in ["Glt", "Glauconite"]:
            if number > 1:
                data = [self.create_glauconite() for n in range(number)]
            else:
                data = self.create_glauconite()
        elif self.mineral in ["Tlc", "Talc"]:
            if number > 1:
                data = [self.create_talc() for n in range(number)]
            else:
                data = self.create_talc()
        elif self.mineral in ["Ctl", "Chrysotile"]:
            if number > 1:
                data = [self.create_chrysotile() for n in range(number)]
            else:
                data = self.create_chrysotile()
        elif self.mineral in ["Ant", "Antigorite"]:
            if number > 1:
                data = [self.create_antigorite() for n in range(number)]
            else:
                data = self.create_antigorite()
        elif self.mineral in ["Prl", "Pyrophyllite"]:
            if number > 1:
                data = [self.create_pyrophyllite() for n in range(number)]
            else:
                data = self.create_pyrophyllite()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        for index in range(number):
            # Clay Minerals
            if self.mineral == "Illite":
                #data_mineral = self.create_illite()
                data_mineral = self.create_illite_simple()
            elif self.mineral == "Kaolinite":
                data_mineral = self.create_kaolinite()
            elif self.mineral == "Montmorillonite":
                data_mineral = self.create_montmorillonite()
            elif self.mineral == "Chamosite":
                data_mineral = self.create_chamosite()
            elif self.mineral == "Clinochlore":
                data_mineral = self.create_clinochlore()
            elif self.mineral == "Pennantite":
                data_mineral = self.create_pennantite()
            elif self.mineral == "Nimite":
                data_mineral = self.create_nimite()
            elif self.mineral == "Chlorite":
                data_mineral = self.create_chlorite()
            elif self.mineral == "Vermiculite":
                data_mineral = self.create_vermiculite()
            elif self.mineral == "Nontronite":
                data_mineral = self.create_nontronite()
            elif self.mineral == "Saponite":
                data_mineral = self.create_saponite()
            # Mica Minerals
            elif self.mineral == "Annite":
                data_mineral = self.create_annite()
            elif self.mineral == "Phlogopite":
                data_mineral = self.create_phlogopite()
            elif self.mineral == "Eastonite":
                data_mineral = self.create_eastonite()
            elif self.mineral == "Siderophyllite":
                data_mineral = self.create_siderophyllite()
            elif self.mineral == "Biotite":
                data_mineral = self.create_biotite()
            elif self.mineral == "Muscovite":
                data_mineral = self.create_muscovite()
            elif self.mineral == "Glauconite":
                data_mineral = self.create_glauconite()
            elif self.mineral == "Talc":
                data_mineral = self.create_talc()
            elif self.mineral == "Chrysotile":
                data_mineral = self.create_chrysotile()
            elif self.mineral == "Antigorite":
                data_mineral = self.create_antigorite()
            elif self.mineral == "Pyrophyllite":
                data_mineral = self.create_pyrophyllite()
            #
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
        #
        return dataset

    def create_illite_simple(self):    # K0.65 Al2 Al0.65 Si3.35 O10 (OH)2
        name = "Ilt"
        # Major Elements
        majors_name = ["H", "O", "Al", "Si", "K"]
        major_compounds = ["H2O", "Al2O3", "SiO2", "K2O"]
        # Trace elements
        element_traces = {
            "2+": ["Mg", "Fe"],
            "All": ["Fe", "Mg"]}

        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        # Molar mass
        molar_mass_pure = round(0.65*self.potassium[2] + (2 + 0.65)*self.aluminium[2] + 3.35*self.silicon[2]
                                + 10*self.oxygen[2] + 2*(self.oxygen[2] + self.hydrogen[2]), 3)
        major_amounts = {"H": 2, "O": 12, "Al": 2.65, "Si": 3.35, "K": 0.65}
        major_masses = {"H": self.hydrogen[2], "O": self.oxygen[2], "Al": self.aluminium[2], "Si": self.silicon[2],
                        "K": self.potassium[2]}
        major_data = {"H": self.hydrogen, "O": self.oxygen, "Al": self.aluminium, "Si": self.silicon,
                      "K": self.potassium}
        condition = False
        while condition == False:
            molar_mass = 0
            amounts = []
            if self.impurity == "impure":
                compositon_data = TraceElements(tracer=self.traces_list).calculate_composition_illite()
                for element in compositon_data:
                    chem_data = PeriodicSystem(name=element).get_data()
                    molar_mass += round(compositon_data[element]["x"]*chem_data[2], 3)
                    amounts.append([chem_data[0], chem_data[1], compositon_data[element]["w"]])
            else:
                for element in majors_name:
                    molar_mass += round(major_amounts[element]*major_masses[element], 3)
                    amounts.append([element, major_data[element][1],
                                    major_amounts[element]*major_masses[element]/molar_mass_pure])
            magic_factor = round(molar_mass/molar_mass_pure, 6)
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            ## Oxide Composition
            list_oxides = ["H2O", "Al2O3", "SiO2", "K2O"]
            if len(self.traces_list) > 0:
                for key in self.traces_list.keys():
                    if key in ["Mg", "Fe"]:
                        list_oxides.append(str(key) + "O")
            composition_oxides = {}
            w_oxide_total = 0
            for index_oxide, var_oxide in enumerate(list_oxides):
                if index_oxide < len(list_oxides) - 1:
                    oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                    if self.impurity == "pure":
                        composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                    else:
                        composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)
                    w_oxide_total += composition_oxides[var_oxide]
                else:
                    if self.impurity == "pure":
                        composition_oxides[var_oxide] = round(1 - w_oxide_total, 4)
                    else:
                        composition_oxides[var_oxide] = round(1 - w_oxide_total, 6)
            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass
        # Density
        dataV = CrystalPhysics([[5.18, 8.98, 10.32], [101.83], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)*magic_factor
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()*magic_factor
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = (35.72 + (62.21-35.72)/(2.706-2.546)*(rho/1000-2.546))*10**9
        # Shear modulus
        G = (17.80 + (25.70-17.80)/(2.706-2.546)*(rho/1000-2.546))*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        results["trace elements"] = element_traces
        results["major elements"] = majors_name
        results["major compounds"] = major_compounds
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_illite(self): # (K,H3O) (Al,Mg,Fe)2 (Si,Al)4 O10 (OH)2
        # Major elements
        majors_name = ["H", "O", "Mg", "Al", "Si", "K", "Fe"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]

        mineral = "Ilt"

        ## Molar mass
        condition = False
        while condition == False:
            a = round(rd.uniform(0.55, 0.65), 4)
            b = round(rd.uniform(0.7, 0.75), 4)
            c = round(rd.uniform(0.875, 0.925), 4)

            mineral = "Ilt"
            majors_data = np.array([
                ["H", self.hydrogen[1], 3*(1 - a) + 2 + 2, self.hydrogen[2]],
                ["O", self.oxygen[1], (1 - a) + 10 + 2 + 1, self.oxygen[2]],
                ["Mg", self.magnesium[1], 2*((1 - b)*0.75), self.magnesium[2]],
                ["Al", self.aluminium[1], 2*b + 4*(1 - c), self.aluminium[2]],
                ["Si", self.silicon[1], 4*c, self.silicon[2]],
                ["K", self.potassium[1], a, self.potassium[2]], ["Fe", self.iron[1], 2*((1 - b)*0.25), self.iron[2]]],
                dtype=object)

            term_01 = a*self.potassium[2] + (1 - a)*(3*self.hydrogen[2] + self.oxygen[2])
            term_02 = 2*(b*self.aluminium[2] + (1 - b)*0.75*self.magnesium[2] + (1 - b)*0.25*self.iron[2])
            term_03 = 4*(c*self.silicon[2] + (1 - c)*self.aluminium[2])
            term_04 = 10*self.oxygen[2] + 2*(self.oxygen[2] + self.hydrogen[2]) + (2*self.hydrogen[2] + self.oxygen[2])

            molar_mass_pure = term_01 + term_02 + term_03 + term_04
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            ## Oxide Composition
            list_oxides = ["H2O", "MgO", "Al2O3", "SiO2", "K2O", "FeO"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()

                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)
            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass
        # Density
        dataV = CrystalPhysics([[5.18, 8.98, 10.32], [101.83], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = (35.72 + (62.21-35.72)/(2.706-2.546)*(rho/1000-2.546))*10**9
        # Shear modulus
        G = (17.80 + (25.70-17.80)/(2.706-2.546)*(rho/1000-2.546))*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_kaolinite(self): # Al2(OH)4Si2O5
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 4, hydrogen[2]], ["O", oxygen[1], 9, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 2, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mg", "Na", "K", "Ti", "Ca"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Kln"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + 4*(oxygen[2]+hydrogen[2]) + 2*silicon[2] + 5*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.13, 8.89, 7.25], [90.0, 104.5, 89.8], "triclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 122.54*10**9
        # Shear modulus
        G = 66.63*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_montmorillonite(self): # (Na,Ca)0.3 (Al,Mg)2 Si4O10 (OH)2 (H2O)10
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "O", "Na", "Mg", "Al", "Si", "Ca"]
        #
        x = round(rd.uniform(0.6, 0.7), 2)
        y = round(rd.uniform(0.9, 1), 2)
        n = rd.randint(8, 12)
        #
        majors_data = np.array([["H", hydrogen[1], 2+2*n, hydrogen[2]], ["O", oxygen[1], 10+2+n, oxygen[2]],
                                ["Na", sodium[1], 0.3*x, sodium[2]], ["Mg", magnesium[1], 2*(1-y), magnesium[2]],
                                ["Al", aluminium[1], 2*y, aluminium[2]], ["Si", silicon[1], 4, silicon[2]],
                                ["Ca", calcium[1], 0.3*(1-x), calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Mnt"
        #
        # Molar mass
        molar_mass_pure = 0.3*(x*sodium[2]+(1-x)*calcium[2]) + 2*(y*aluminium[2]+(1-y)*magnesium[2]) + 4*silicon[2] \
                          + 10*oxygen[2] + 2*(hydrogen[2]+oxygen[2]) + n*(2*hydrogen[2]+oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.17, 8.94, 9.95], [99.9], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 1
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2738, 2788, 3250, 2504, 3182]
        y_K = [37.30, 35.31, 49.46, 29.71, 66.59]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        y_G = [17.00, 20.19, 24.70, 16.30, 27.00]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_nontronite(self): # Na0.3 Fe2 (Si,Al)4 O10 (OH)2 * nH2O
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Na", "Al", "Si", "Fe"]
        #
        x = round(rd.uniform(0.0, 0.5), 2)
        n = rd.randint(1, 10)
        #
        majors_data = np.array([["H", hydrogen[1], 2+2*n, hydrogen[2]], ["O", oxygen[1], 10+2+n, oxygen[2]],
                                ["Na", sodium[1], 0.3, sodium[2]], ["Al", aluminium[1], 4*x, aluminium[2]],
                                ["Si", silicon[1], 4*(1-x), silicon[2]], ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mg", "Ca"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Nntr"
        #
        # Molar mass
        molar_mass_pure = 0.3*sodium[2] + 2*iron[2] + 4*((1-x)*silicon[2] + x*aluminium[2]) + 10*oxygen[2] \
                          + 2*(oxygen[2] + hydrogen[2]) + n*(2*hydrogen[2] + oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.23, 9.11, 15.25], [96.0], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2738, 2788, 3250, 2504, 3182]
        y_K = [37.30, 35.31, 49.46, 29.71, 66.59]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        y_G = [17.00, 20.19, 24.70, 16.30, 27.00]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_saponite(self): # Ca0.25 (Mg,Fe)3 (Si,Al)4 O10 (OH)2 * nH2O
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Mg", "Al", "Si", "Ca", "Fe"]
        #
        x = round(rd.uniform(0.0, 0.5), 2)
        y = round(rd.uniform(0.0, 0.75), 2)
        n = rd.randint(1, 5)
        #
        majors_data = np.array(
            [["H", hydrogen[1], 2+2*n, hydrogen[2]], ["O", oxygen[1], 10+2+n, oxygen[2]],
             ["Mg", magnesium[1], 3*y, magnesium[2]], ["Al", aluminium[1], 4*x, aluminium[2]],
             ["Si", silicon[1], 4*(1-x), silicon[2]], ["Ca", calcium[1], 0.25, calcium[2]],
             ["Fe", iron[1], 3*(1-y), iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Ni", "K", "P"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Sap"
        #
        # Molar mass
        molar_mass_pure = 0.25*calcium[2] + 3*((1-y)*iron[2] + y*magnesium[2]) + 4*((1-x)*silicon[2] + x*aluminium[2]) \
                          + 10*oxygen[2] + 2*(oxygen[2] + hydrogen[2]) + n*(2*hydrogen[2] + oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.3, 9.16, 12.4], [96.5], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2738, 2788, 3250, 2504, 3182]
        y_K = [37.30, 35.31, 49.46, 29.71, 66.59]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        y_G = [17.00, 20.19, 24.70, 16.30, 27.00]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_chamosite(self): # Fe5 Al2 Si3 O10 (OH)8
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Al", "Si", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 8, hydrogen[2]], ["O", oxygen[1], 18, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Fe", iron[1], 5, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Ca", "Na", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Chl"
        #
        # Molar mass
        molar_mass_pure = 5*iron[2] + 2*aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 8*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.373, 9.306, 14.222], [97.88], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_K = [82.6, 83.9, 95.9, 106.5, 123.6, 134.6]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_G = [45.7, 46.8, 46.5, 43.5, 39.3, 36.9]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_clinochlore(self): # Mg5 Al2 Si3 O10 (OH)8
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Mg", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 8, hydrogen[2]], ["O", oxygen[1], 18, oxygen[2]],
                                ["Mg", magnesium[1], 5, magnesium[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Zn", "Ca", "Cr"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Chl"
        #
        # Molar mass
        molar_mass_pure = 5*magnesium[2] + 2*aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 8*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.3, 9.3, 14.3], [97], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_K = [82.6, 83.9, 95.9, 106.5, 123.6, 134.6]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_G = [45.7, 46.8, 46.5, 43.5, 39.3, 36.9]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_pennantite(self): # Mn5 Al2 Si3 O10 (OH)8
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["H", "O", "Al", "Si", "Mn"]
        #
        majors_data = np.array([["H", hydrogen[1], 8, hydrogen[2]], ["O", oxygen[1], 18, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Mn", manganese[1], 5, manganese[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Zn", "Mg", "Ba","H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Chl"
        #
        # Molar mass
        molar_mass_pure = 5*manganese[2] + 2*aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 8*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.454, 9.45, 14.4], [97.2], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_K = [82.6, 83.9, 95.9, 106.5, 123.6, 134.6]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_G = [45.7, 46.8, 46.5, 43.5, 39.3, 36.9]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_nimite(self): # Ni5 Al2 Si3 O10 (OH)8
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["H", "O", "Al", "Si", "Ni"]
        #
        majors_data = np.array([["H", hydrogen[1], 8, hydrogen[2]], ["O", oxygen[1], 18, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ni", nickel[1], 5, nickel[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Cr", "Mn", "Co", "Ca", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Chl"
        #
        # Molar mass
        molar_mass_pure = 5*nickel[2] + 2*aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 8*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.32, 9.214, 14.302], [97.1], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_K = [82.6, 83.9, 95.9, 106.5, 123.6, 134.6]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_G = [45.7, 46.8, 46.5, 43.5, 39.3, 36.9]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_chlorite(self): # (Fe, Mg, Mn,Ni)5 Al2 Si3 O10 (OH)8
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["H", "O", "Mg", "Al", "Si", "Mn", "Fe", "Ni"]
        #
        x = round(rd.uniform(0, 1), 2)
        y = round(rd.uniform(0, (1-x)), 2)
        z = round(rd.uniform(0, (1-x-y)), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 8, hydrogen[2]], ["O", oxygen[1], 18, oxygen[2]],
                                ["Mg", magnesium[1], 5*y, magnesium[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["Mn", manganese[1], 5*z, manganese[2]],
                                ["Fe", iron[1], 5*x, iron[2]], ["Ni", nickel[1], 5*(1-x-y-z), nickel[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Cr", "Mn", "Co", "Ca", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Chl"
        #
        # Molar mass
        molar_mass_pure = 5*(x*iron[2] + y*magnesium[2] + z*manganese[2] + (1-x-y-z)*nickel[2]) + 2*aluminium[2] \
                          + 3*silicon[2] + 10*oxygen[2] + 8*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[5.373, 9.306, 14.222], [97.88], "monoclinic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 2
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Mg = CrystalPhysics([[5.3, 9.3, 14.3], [97], "monoclinic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 2
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[5.454, 9.45, 14.4], [97.2], "monoclinic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 2
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()
        #
        dataV_Ni = CrystalPhysics([[5.32, 9.214, 14.302], [97.1], "monoclinic"])
        V_Ni = dataV_Ni.calculate_volume()
        Z_Ni = 2
        V_m_Ni = MineralChemistry().calculate_molar_volume(volume_cell=V_Ni, z=Z_Ni)
        dataRho_Ni = CrystalPhysics([molar_mass, Z_Ni, V_Ni])
        rho_Ni = dataRho_Ni.calculate_bulk_density()
        rho_e_Ni = wg(amounts=amounts, elements=element, rho_b=rho_Ni).calculate_electron_density()
        #
        V_m = x*V_m_Fe + y*V_m_Mg + z*V_m_Mn + (1-x-y-z)*V_m_Ni
        rho = x*rho_Fe + y*rho_Mg + z*rho_Mn + (1-x-y-z)*rho_Ni
        rho_e = x*rho_e_Fe + y*rho_e_Mg + z*rho_e_Mn + (1-x-y-z)*rho_e_Ni
        # Bulk modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_K = [82.6, 83.9, 95.9, 106.5, 123.6, 134.6]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_G = [45.7, 46.8, 46.5, 43.5, 39.3, 36.9]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_vermiculite(self): # (Mg,Fe,Al)3 (Al,Si)4 O10 (OH)2 4(H2O)
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Mg", "Al", "Si", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        y = round(rd.uniform(0, (1-x)), 2)
        z = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2+4*2, hydrogen[2]], ["O", oxygen[1], 10+2+4, oxygen[2]],
                                ["Mg", magnesium[1], 3*x, magnesium[2]], ["Al", aluminium[1], 3*(1-x-y)+4*z, aluminium[2]],
                                ["Si", silicon[1], 4*(1-z), silicon[2]], ["Fe", iron[1], 3*y, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ca", "Na", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Vrm"
        #
        # Molar mass
        molar_mass_pure = 3*(x*magnesium[2] + y*iron[2] + (1-x-y)*aluminium[2]) + 4*(z*aluminium[2] + (1-z)*silicon[2]) \
                          + 10*oxygen[2] + 2*(oxygen[2]+hydrogen[2]) + 4*(2*hydrogen[2]+oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.26, 9.23, 14.97], [96.82], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_K = [82.6, 83.9, 95.9, 106.5, 123.6, 134.6]
        a_K, b_K, r_value_K, p_value_K, std_err_K = stats.linregress(x_rho, y_K)
        K = (a_K*rho + b_K)*10**9
        # Shear modulus
        x_rho = [2530, 2600, 2670, 2750, 2840, 2930]
        y_G = [45.7, 46.8, 46.5, 43.5, 39.3, 36.9]
        a_G, b_G, r_value_G, p_value_G, std_err_G = stats.linregress(x_rho, y_G)
        G = (a_G*rho + b_G)*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_annite(self): # K Fe3 Al Si3 O10 (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Al", "Si", "K", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 10+2, oxygen[2]],
                                ["Al", aluminium[1], 1, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["K", potassium[1], 1, potassium[2]],
                                ["Fe", iron[1], 3, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Mg", "Ca", "Na", "Cl"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ann"
        #
        # Molar mass K Fe3 Al Si3 O10 (OH)2
        molar_mass_pure = potassium[2] + 3*iron[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(oxygen[2]+hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.39, 9.334, 10.29], [100], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 114.72*10**9
        # Shear modulus
        G = 58.61*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_phlogopite(self): # K Mg3 Al Si3 O10 (F, OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        flourine = PeriodicSystem(name="F").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        majors_name = ["H", "O", "F", "Mg", "Al", "Si", "K"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2*(1-x), hydrogen[2]], ["O", oxygen[1], 10+2*(1-x), oxygen[2]],
                                ["F", flourine[1], 2*x, flourine[2]], ["Mg", magnesium[1], 3, magnesium[2]],
                                ["Al", aluminium[1], 1, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["K", potassium[1], 1, potassium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Ba", "Cr", "Na", "Ti", "Ni", "Zn", "Ca", "Li", "Rb"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Phl"
        #
        # Molar mass
        molar_mass_pure = potassium[2] + 3*magnesium[2] + aluminium[2] + 3*silicon[2] + 10*oxygen[2] + 2*(x*flourine[2] + (1-x)*(oxygen[2]+hydrogen[2]))
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.39, 9.334, 10.29], [100], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 103.78*10**9
        # Shear modulus
        G = 61.69*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_eastonite(self): # K Mg2 Al3 Si2 O10 (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        majors_name = ["H", "O", "Mg", "Al", "Si", "K"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 10+2, oxygen[2]],
                                ["Mg", magnesium[1], 2, magnesium[2]], ["Al", aluminium[1], 3, aluminium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["K", potassium[1], 1, potassium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Eas"
        #
        # Molar mass K Mg2 Al3 Si2 O10 (OH)2
        molar_mass_pure = potassium[2] + 2*magnesium[2] + 3*aluminium[2] + 2*silicon[2] + 10*oxygen[2] + 2*(oxygen[2]+hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.27, 9.13, 10.15], [99.54], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 100.44*10**9
        # Shear modulus
        G = 62.77*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_siderophyllite(self): # K Fe2 Al3 Si2 O10 (F, OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        flourine = PeriodicSystem(name="F").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "F", "Al", "Si", "K", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2*(1-x), hydrogen[2]], ["O", oxygen[1], 10+2*(1-x), oxygen[2]],
                                ["F", flourine[1], 2*x, flourine[2]], ["Al", aluminium[1], 3, aluminium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["K", potassium[1], 1, potassium[2]],
                                ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Mg", "Ca", "Li", "Na", "Rb", "Cs", "Cl"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Sdp"
        #
        # Molar mass
        molar_mass_pure = potassium[2] + 2*iron[2] + 3*aluminium[2] + 2*silicon[2] + 10*oxygen[2] + 2*(x*flourine[2] + (1-x)*(oxygen[2]+hydrogen[2]))
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.369, 9.297, 10.268], [100.06], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 110.91*10**9
        # Shear modulus
        G = 59.61*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results

    def create_biotite(self, enrichment=None, x_value=None, y_value=None):
        # Major Elements
        majors_name = ["H", "O", "Mg", "Al", "Si", "K", "Fe"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = minors = ["Ti", "Mn", "Mg", "Ca", "Li", "Na", "Rb", "Cs", "Cl"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        condition = False
        while condition == False:
            if x_value == None and y_value == None:
                if enrichment == None:
                    x = round(rd.uniform(0, 1), 2)
                    y = round(rd.uniform(0, 1), 2)
                elif enrichment == "Na":
                    x = round(rd.uniform(0, 1), 2)
                    y = round(rd.uniform(0, 1), 2)
                elif enrichment == "Ca":
                    x = round(rd.uniform(0, 1), 2)
                    y = round(rd.uniform(0, 1), 2)
            else:
                x = x_value
                y = y_value

            mineral = "Bt"
            majors_data = np.array(
                [["H", self.hydrogen[1], 2, self.hydrogen[2]], ["O", self.oxygen[1], 12, self.oxygen[2]],
                 ["Mg", self.magnesium[1], (1 - x)*(2 + y), self.magnesium[2]],
                 ["Al", self.aluminium[1], (3 - 2*y), self.aluminium[2]],
                 ["Si", self.silicon[1], (2 + y), self.silicon[2]], ["K", self.potassium[1], 1, self.potassium[2]],
                 ["Fe", self.iron[1], x*(2 + y), self.iron[2]]], dtype=object)
            majors_data_Ann = np.array(
                [["H", self.hydrogen[1], 2, self.hydrogen[2]], ["O", self.oxygen[1], 12, self.oxygen[2]],
                 ["Mg", self.magnesium[1], 0, self.magnesium[2]], ["Al", self.aluminium[1], 1, self.aluminium[2]],
                 ["Si", self.silicon[1], 3, self.silicon[2]], ["K", self.potassium[1], 1, self.potassium[2]],
                 ["Fe", self.iron[1], 3, self.iron[2]]], dtype=object)
            majors_data_Phl = np.array(
                [["H", self.hydrogen[1], 2, self.hydrogen[2]], ["O", self.oxygen[1], 12, self.oxygen[2]],
                 ["Mg", self.magnesium[1], 3, self.magnesium[2]], ["Al", self.aluminium[1], 1, self.aluminium[2]],
                 ["Si", self.silicon[1], 3, self.silicon[2]], ["K", self.potassium[1], 1, self.potassium[2]],
                 ["Fe", self.iron[1], 0, self.iron[2]]], dtype=object)
            majors_data_Sdp = np.array(
                [["H", self.hydrogen[1], 2, self.hydrogen[2]], ["O", self.oxygen[1], 12, self.oxygen[2]],
                 ["Mg", self.magnesium[1], 0, self.magnesium[2]], ["Al", self.aluminium[1], 3, self.aluminium[2]],
                 ["Si", self.silicon[1], 2, self.silicon[2]], ["K", self.potassium[1], 1, self.potassium[2]],
                 ["Fe", self.iron[1], 2, self.iron[2]]], dtype=object)
            majors_data_Eas = np.array(
                [["H", self.hydrogen[1], 2, self.hydrogen[2]], ["O", self.oxygen[1], 12, self.oxygen[2]],
                 ["Mg", self.magnesium[1], 2, self.magnesium[2]], ["Al", self.aluminium[1], 3, self.aluminium[2]],
                 ["Si", self.silicon[1], 2, self.silicon[2]], ["K", self.potassium[1], 1, self.potassium[2]],
                 ["Fe", self.iron[1], 0, self.iron[2]]], dtype=object)

            molar_mass_pure = (self.potassium[2] + (2 + y)*(x*self.iron[2] + (1 - x)*self.magnesium[2]) +
                               (3 - 2*y)*self.aluminium[2] + (2 + y)*self.silicon[2] + 10*self.oxygen[2] +
                               2*(self.oxygen[2] + self.hydrogen[2]))
            molar_mass_pure_Ann = (self.potassium[2] + 3*self.iron[2] + self.aluminium[2] + 3*self.silicon[2] +
                                   10*self.oxygen[2] + 2*(self.oxygen[2] + self.hydrogen[2]))
            molar_mass_pure_Phl = (self.potassium[2] + 3*self.magnesium[2] + self.aluminium[2] + 3*self.silicon[2] +
                                   10*self.oxygen[2] + 2*(self.oxygen[2] + self.hydrogen[2]))
            molar_mass_pure_Sdp = (self.potassium[2] + 2*self.iron[2] + 3*self.aluminium[2] + 2*self.silicon[2] +
                                   10*self.oxygen[2] + 2*(self.oxygen[2] + self.hydrogen[2]))
            molar_mass_pure_Eas = (self.potassium[2] + 2*self.magnesium[2] + 3*self.aluminium[2] + 2*self.silicon[2] +
                                   10*self.oxygen[2] + 2*(self.oxygen[2] + self.hydrogen[2]))

            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
            #
            molar_mass_Ann, amounts_Ann = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Ann, majors=majors_data_Ann).calculate_molar_mass()
            element_Ann = [PeriodicSystem(name=amounts_Ann[i][0]).get_data() for i in range(len(amounts_Ann))]
            #
            molar_mass_Phl, amounts_Phl = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Phl, majors=majors_data_Phl).calculate_molar_mass()
            element_Phl = [PeriodicSystem(name=amounts_Phl[i][0]).get_data() for i in range(len(amounts_Phl))]
            #
            molar_mass_Sdp, amounts_Sdp = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Sdp, majors=majors_data_Sdp).calculate_molar_mass()
            element_Sdp = [PeriodicSystem(name=amounts_Sdp[i][0]).get_data() for i in range(len(amounts_Sdp))]
            #
            molar_mass_Eas, amounts_Eas = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_Eas, majors=majors_data_Eas).calculate_molar_mass()
            element_Eas = [PeriodicSystem(name=amounts_Eas[i][0]).get_data() for i in range(len(amounts_Eas))]

            ## Oxide Composition
            list_oxides = ["Al2O3", "SiO2", "H2O", "MgO", "K2O", "FeO"]
            composition_oxides = {}
            for var_oxide in list_oxides:
                oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                if self.impurity == "pure":
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                else:
                    composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)

            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass

        # Density
        dataV_Ann = CrystalPhysics([[5.39, 9.334, 10.29], [100], "monoclinic"])
        V_Ann = dataV_Ann.calculate_volume()
        Z_Ann = 2
        V_m_Ann = MineralChemistry().calculate_molar_volume(volume_cell=V_Ann, z=Z_Ann)
        dataRho_Ann = CrystalPhysics([molar_mass_Ann, Z_Ann, V_Ann])
        rho_Ann = dataRho_Ann.calculate_bulk_density()
        rho_e_Ann = wg(amounts=amounts_Ann, elements=element_Ann, rho_b=rho_Ann).calculate_electron_density()

        dataV_Phl = CrystalPhysics([[5.39, 9.334, 10.29], [100], "monoclinic"])
        V_Phl = dataV_Phl.calculate_volume()
        Z_Phl = 2
        V_m_Phl = MineralChemistry().calculate_molar_volume(volume_cell=V_Phl, z=Z_Phl)
        dataRho_Phl = CrystalPhysics([molar_mass_Phl, Z_Phl, V_Phl])
        rho_Phl = dataRho_Phl.calculate_bulk_density()
        rho_e_Phl = wg(amounts=amounts_Phl, elements=element_Phl, rho_b=rho_Phl).calculate_electron_density()

        dataV_Eas = CrystalPhysics([[5.27, 9.13, 10.15], [99.54], "monoclinic"])
        V_Eas = dataV_Eas.calculate_volume()
        Z_Eas = 2
        V_m_Eas = MineralChemistry().calculate_molar_volume(volume_cell=V_Eas, z=Z_Eas)
        dataRho_Eas = CrystalPhysics([molar_mass_Eas, Z_Eas, V_Eas])
        rho_Eas = dataRho_Eas.calculate_bulk_density()
        rho_e_Eas = wg(amounts=amounts_Eas, elements=element_Eas, rho_b=rho_Eas).calculate_electron_density()

        dataV_Sdp = CrystalPhysics([[5.369, 9.297, 10.268], [100.06], "monoclinic"])
        V_Sdp = dataV_Sdp.calculate_volume()
        Z_Sdp = 2
        V_m_Sdp = MineralChemistry().calculate_molar_volume(volume_cell=V_Sdp, z=Z_Sdp)
        dataRho_Sdp = CrystalPhysics([molar_mass_Sdp, Z_Sdp, V_Sdp])
        rho_Sdp = dataRho_Sdp.calculate_bulk_density()
        rho_e_Sdp = wg(amounts=amounts_Sdp, elements=element_Sdp, rho_b=rho_Sdp).calculate_electron_density()

        V_m = x*(y*V_m_Ann + (1 - y)*V_m_Sdp) + (1 - x)*(y*V_m_Phl + (1 - y)*V_m_Eas)
        rho = x*(y*rho_Ann + (1 - y)*rho_Sdp) + (1 - x)*(y*rho_Phl + (1 - y)*rho_Eas)
        rho_e = x*(y*rho_e_Ann + (1 - y)*rho_e_Sdp) + (1 - x)*(y*rho_e_Phl + (1 - y)*rho_e_Eas)
        # Bulk modulus
        K_Ann = 114.72*10**9
        K_Phl = 103.78*10**9
        K_Sdp = 110.91*10**9
        K_Eas = 100.44*10**9
        K = x*(y*K_Ann + (1 - y)*K_Sdp) + (1 - x)*(y*K_Phl + (1 - y)*K_Eas)
        # Shear modulus
        G_Ann = 58.61*10**9
        G_Phl = 61.69*10**9
        G_Sdp = 59.61*10**9
        G_Eas = 62.77*10**9
        G = x*(y*G_Ann + (1 - y)*G_Sdp) + (1 - x)*(y*G_Phl + (1 - y)*G_Eas)
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        # RESULTS
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results

    def create_muscovite(self): # K Al3 Si3 O10 (F,OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        fluorine = PeriodicSystem(name="F").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        majors_name = ["H", "O", "F", "Al", "Si", "K"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2*(1-x), hydrogen[2]], ["O", oxygen[1], 10+2*(1-x), oxygen[2]],
                                ["F", fluorine[1], 2*x, fluorine[2]], ["Al", aluminium[1], 3, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["K", potassium[1], 1, potassium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Cr", "Li", "Fe", "V", "Mn", "Na", "Cs", "Rb", "Ca", "Mg", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ms"
        #
        # Molar mass K Al3 Si3 O10 (F,OH)2
        molar_mass_pure = potassium[2] + 3*aluminium[2] + 3*silicon[2] + 10*oxygen[2] \
                          + 2*(x*fluorine[2] + (1-x)*(oxygen[2]+hydrogen[2]))
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.19, 9.03, 20.05], [95.5], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 112.65*10**9
        # Shear modulus
        G = 68.35*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_glauconite(self): # (K,Na) (Fe,Al,Mg)2 (Si,Al)4 O10 (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Na", "Mg", "Al", "Si", "K", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        y1 = round(rd.uniform(0, 1), 2)
        y2 = round(rd.uniform(0, (1-y1)), 2)
        z = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 10+2, oxygen[2]],
                                ["Na", sodium[1], (1-x), sodium[2]], ["Mg", magnesium[1], 2*(1-y1-y2), magnesium[2]],
                                ["Al", aluminium[1], 2*y2+4*(1-z), aluminium[2]], ["Si", silicon[1], 4*z, silicon[2]],
                                ["K", potassium[1], x, potassium[2]], ["Fe", iron[1], 2*y1, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Ca", "P"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Glt"
        #
        # Molar mass (K,Na) (Fe,Al,Mg)2 (Si,Al)4 O10 (OH)2
        molar_mass_pure = (x*potassium[2] + (1-x)*sodium[2]) + 2*(y1*iron[2] + y2*aluminium[2] + (1-y1-y2)*magnesium[2]) \
                          + 4*(z*silicon[2] + (1-z)*aluminium[2]) + 10*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.234, 9.066, 10.16], [100.5], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 112.64*10**9
        # Shear modulus
        G = 68.33*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_talc(self): # Mg3 Si4 O10 (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 12, oxygen[2]],
                                ["Mg", magnesium[1], 3, magnesium[2]], ["Si", silicon[1], 4, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ni", "Fe", "Al", "Ca", "Na"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Tlc"
        #
        # Molar mass
        molar_mass_pure = 3*magnesium[2] + 4*silicon[2] + 10*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.27, 9.12, 18.85], [100.016], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 47*10**9
        # Shear modulus
        G = 25*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_chrysotile(self): # Mg3 Si2 O5 (OH)4
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 4, hydrogen[2]], ["O", oxygen[1], 9, oxygen[2]],
                                ["Mg", magnesium[1], 3, magnesium[2]], ["Si", silicon[1], 2, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ctl"
        #
        # Molar mass
        molar_mass_pure = 3*magnesium[2] + 2*silicon[2] + 5*oxygen[2] + 4*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        #dataV = CrystalPhysics([[5.313, 9.12, 14.637], [93.167], "monoclinic"])
        dataV = CrystalPhysics([[5.371, 14.807], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])    # *10**(6) --> if hexagonal is used
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 108.82*10**9
        # Shear modulus
        G = 58.17*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_antigorite(self): # Mg3 Si2 O5 (OH)4
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 4, hydrogen[2]], ["O", oxygen[1], 9, oxygen[2]],
                                ["Mg", magnesium[1], 3, magnesium[2]], ["Si", silicon[1], 2, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ni", "Al", "Mn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ant"
        #
        # Molar mass
        molar_mass_pure = 3*magnesium[2] + 2*silicon[2] + 5*oxygen[2] + 4*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        #dataV = CrystalPhysics([[5.313, 9.12, 14.637], [93.167], "monoclinic"])
        dataV = CrystalPhysics([[5.371, 14.807], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])    # *10**(6) --> if hexagonal is used
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 108.82*10**9
        # Shear modulus
        G = 58.17*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_pyrophyllite(self): # Al2 Si4 O10 (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 12, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 4, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Prl"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + 4*silicon[2] + 10*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.161, 8.957, 9.351], [91.03, 100.37, 89.75], "triclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 74*10**9
        # Shear modulus
        G = 48*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["state"] = var_state
            results["M"] = round(molar_mass, 3)
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
# NESOSILICATES
class Nesosilicates:
    """ Class that generates geophysical and geochemical data of nesosilicate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1):
        # ["Zircon", "Thorite", "Topaz", "Pyrope", "Olivine", "Ca-Garnet", "Al-Garnet", "Uvarovite", "Andratite",
        # "Grossular", "Almandine", "Liebenbergite", "Ca-Olivine", "Tephroite", "Fayalite", "Forsterite", "Staurolite",
        # "Sillimanite", "Kyanite", "Andalusite", "Spessartine"]
        if self.mineral in ["Zr", "Zircon"]:
            if number > 1:
                data = [self.create_zircon() for n in range(number)]
            else:
                data = self.create_zircon()
        elif self.mineral in ["Thr", "Thorite"]:
            if number > 1:
                data = [self.create_thorite() for n in range(number)]
            else:
                data = self.create_thorite()
        elif self.mineral in ["Tpz", "Topaz"]:
            if number > 1:
                data = [self.create_topaz() for n in range(number)]
            else:
                data = self.create_topaz()
        elif self.mineral in ["Prp", "Pyrope"]:
            if number > 1:
                data = [self.create_pyrope() for n in range(number)]
            else:
                data = self.create_pyrope()
        elif self.mineral in ["Ol", "Olivine"]:
            if number > 1:
                data = [self.create_olivine() for n in range(number)]
            else:
                data = self.create_olivine()
        elif self.mineral in ["Ca-Grt", "Ca-Garnet"]:
            if number > 1:
                data = [self.create_calcium_garnet() for n in range(number)]
            else:
                data = self.create_calcium_garnet()
        elif self.mineral in ["Al-Grt", "Al-Garnet"]:
            if number > 1:
                data = [self.create_aluminium_garnet() for n in range(number)]
            else:
                data = self.create_aluminium_garnet()
        elif self.mineral in ["Uv", "Uvarovite"]:
            if number > 1:
                data = [self.create_uvarovite() for n in range(number)]
            else:
                data = self.create_uvarovite()
        elif self.mineral in ["Adr", "Andratite"]:
            if number > 1:
                data = [self.create_andradite() for n in range(number)]
            else:
                data = self.create_andradite()
        elif self.mineral in ["Grs", "Grossular"]:
            if number > 1:
                data = [self.create_grossular() for n in range(number)]
            else:
                data = self.create_grossular()
        elif self.mineral in ["Alm", "Almandine"]:
            if number > 1:
                data = [self.create_almandine() for n in range(number)]
            else:
                data = self.create_almandine()
        elif self.mineral in ["Lbg", "Liebenbergite"]:
            if number > 1:
                data = [self.create_liebenbergite() for n in range(number)]
            else:
                data = self.create_liebenbergite()
        elif self.mineral in ["Ca-Ol", "Ca-Olivine"]:
            if number > 1:
                data = [self.create_calcio_olivine() for n in range(number)]
            else:
                data = self.create_calcio_olivine()
        elif self.mineral in ["Tep", "Tephroite"]:
            if number > 1:
                data = [self.create_tephroite() for n in range(number)]
            else:
                data = self.create_tephroite()
        elif self.mineral in ["Fa", "Fayalite"]:
            if number > 1:
                data = [self.create_fayalite() for n in range(number)]
            else:
                data = self.create_fayalite()
        elif self.mineral in ["Fo", "Forsterite"]:
            if number > 1:
                data = [self.create_forsterite() for n in range(number)]
            else:
                data = self.create_forsterite()
        elif self.mineral in ["St", "Staurolite"]:
            if number > 1:
                data = [self.create_staurolite() for n in range(number)]
            else:
                data = self.create_staurolite()
        elif self.mineral in ["Sil", "Sillimanite"]:
            if number > 1:
                data = [self.create_sillimanite() for n in range(number)]
            else:
                data = self.create_sillimanite()
        elif self.mineral in ["Ky", "Kyanite"]:
            if number > 1:
                data = [self.create_kyanite() for n in range(number)]
            else:
                data = self.create_kyanite()
        elif self.mineral in ["And", "Andalusite"]:
            if number > 1:
                data = [self.create_andalusite() for n in range(number)]
            else:
                data = self.create_andalusite()
        elif self.mineral in ["Sps", "Spessartine"]:
            if number > 1:
                data = [self.create_spessartine() for n in range(number)]
            else:
                data = self.create_spessartine()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Zircon":
                data_mineral = self.create_zircon()
            elif self.mineral == "Titanite":
                data_mineral = self.create_titanite()
            elif self.mineral == "Thorite":
                data_mineral = self.create_thorite()
            elif self.mineral == "Andalusite":
                data_mineral = self.create_andalusite()
            elif self.mineral == "Kyanite":
                data_mineral = self.create_kyanite()
            elif self.mineral == "Sillimanite":
                data_mineral = self.create_sillimanite()
            elif self.mineral == "Topaz":
                data_mineral = self.create_topaz()
            elif self.mineral == "Staurolite":
                data_mineral = self.create_staurolite()
            elif self.mineral == "Fayalite":
                data_mineral = self.create_fayalite()
            elif self.mineral == "Forsterite":
                data_mineral = self.create_forsterite()
            elif self.mineral == "Tephroite":
                data_mineral = self.create_tephroite()
            elif self.mineral == "Ca-Olivine":
                data_mineral = self.create_calcio_olivine()
            elif self.mineral == "Liebenbergite":
                data_mineral = self.create_liebenbergite()
            elif self.mineral == "Olivine":
                data_mineral = self.create_olivine()
            elif self.mineral == "Pyrope":
                data_mineral = self.create_pyrope()
            elif self.mineral == "Almandine":
                data_mineral = self.create_almandine()
            elif self.mineral == "Grossular":
                data_mineral = self.create_grossular()
            elif self.mineral == "Anhadrite":
                data_mineral = self.create_andradite()
            elif self.mineral == "Uvarovite":
                data_mineral = self.create_uvarovite()
            elif self.mineral == "Al-Garnet":
                data_mineral = self.create_aluminium_garnet()
            elif self.mineral == "Ca-Garnet":
                data_mineral = self.create_calcium_garnet()
            #
            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U",
                           "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif key in ["mineral", "state", "trace elements"] and key not in dataset:
                    dataset[key] = value
                elif key in ["chemistry"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)
        #
        return dataset
    #
    def create_zircon(self): # Zr Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        zirconium = PeriodicSystem(name="Zr").get_data()
        majors_name = ["O", "Si", "Zr"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Zr", zirconium[1], 1, zirconium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Hf", "Th", "U", "REE", "H", "Fe", "Al", "P"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Zr"
        #
        # Molar mass
        molar_mass_pure = zirconium[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[6.604, 5.979], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 203*10**9
        # Shear modulus
        G = 96*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_titanite(self): # Ca Ti Si O5
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        majors_name = ["O", "Si", "Ca", "Ti"]
        #
        majors_data = np.array([["O", oxygen[1], 5, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Ca", calcium[1], 1, calcium[2]], ["Ti", titanium[1], 1, titanium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Y", "Mn", "Al", "Ce", "Sr", "Na", "Nb", "Ta", "Al", "Mg", "V", "F", "Zr", "Sn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ttn"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + titanium[2] + silicon[2] + 5*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[6.56, 8.72, 7.44], [119.716], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 102*10**9
        # Shear modulus
        G = 63*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_thorite(self): # Th Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        thorium = PeriodicSystem(name="Th").get_data()
        majors_name = ["O", "Si", "Th"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Th", thorium[1], 1, thorium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Al", "Fe", "Pb", "Ca", "P", "Ti", "REE", "Y", "Mg", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Thr"
        #
        # Molar mass
        molar_mass_pure = thorium[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.117, 6.295], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 170*10**9
        # Shear modulus
        G = 63*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_andalusite(self): # Al2 Si O5
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Al", "Si"]
        #
        majors_data = np.array([["O", oxygen[1], 5, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 1, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "And"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + silicon[2] + 5*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.79, 7.9, 5.56], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 142*10**9
        # Shear modulus
        G = 89*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_kyanite(self): # Al2 Si O5
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Al", "Si"]
        #
        majors_data = np.array([["O", oxygen[1], 5, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 1, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ky"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + silicon[2] + 5*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.112, 7.844, 5.574], [90.12, 101.1, 105.9], "triclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 191*10**9
        # Shear modulus
        G = 120.85*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_sillimanite(self): # Al2 Si O5
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Al", "Si"]
        #
        majors_data = np.array([["O", oxygen[1], 5, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 1, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Sil"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + silicon[2] + 5*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.484, 7.672, 5.77], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 162*10**9
        # Shear modulus
        G = 101.5*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_topaz(self): # Al2 Si O4 (F,OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        flourine = PeriodicSystem(name="F").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "F", "Al", "Si"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2*(1-x), hydrogen[2]], ["O", oxygen[1], 4+2*(1-x), oxygen[2]],
                                ["F", flourine[1], 2*x, flourine[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 1, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Tpz"
        #
        # Molar mass
        molar_mass_pure = 2*aluminium[2] + silicon[2] + 4*oxygen[2] + 2*(x*flourine[2] + (1-x)*(oxygen[2]+hydrogen[2]))
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.6499, 8.796, 8.3909], [], "orthorhombic"])   # Ref. Handbook of Mineralogy
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K_F = 170.59*10**9
        K_OH = 176.53*10**9
        K = x*K_F + (1-x)*K_OH
        # Shear modulus
        G_F = 97.71*10**9
        G_OH = 95.28*10**9
        G = x*G_F + (1-x)*G_OH
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_staurolite(self): # Fe2 Al9 Si4 O23 (OH)
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Al", "Si", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 23+1, oxygen[2]],
                                ["Al", aluminium[1], 9, aluminium[2]], ["Si", silicon[1], 4, silicon[2]],
                                ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Cr", "Mn", "Co", "Zn", "Li"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "St"
        #
        # Molar mass
        molar_mass_pure = 2*iron[2] + 9*aluminium[2] + 4*silicon[2] + 23*oxygen[2] + (oxygen[2]+hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.88, 16.62, 5.66], [90], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 181.2*10**9
        # Shear modulus
        G = 100.77*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_fayalite(self): # Fe2 Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Si", "Fe"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Fa"
        #
        # Molar mass
        molar_mass_pure = 2*iron[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.76, 10.2, 5.98], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 143.76*10**9
        # Shear modulus
        G = 65.66*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_forsterite(self): # Mg2 Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Mg", "Si"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], 2, magnesium[2]],
                                ["Si", silicon[1], 1, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Fo"
        #
        # Molar mass
        molar_mass_pure = 2*magnesium[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.756, 10.195, 5.981], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 119*10**9
        # Shear modulus
        G = 74*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_tephroite(self): # Mn2 Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Si", "Mn"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Mn", manganese[1], 2, manganese[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Zn", "Ca", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Tep"
        #
        # Molar mass
        molar_mass_pure = 2*manganese[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.76, 10.2, 5.98], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 116*10**9
        # Shear modulus
        G = 50*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_calcio_olivine(self): # Ca2 Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Si", "Ca"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Ca", calcium[1], 2, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "CaOl"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.074, 11.211, 6.753], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 90*10**9
        # Shear modulus
        G = 49*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_liebenbergite(self): # Ni2 Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["O", "Si", "Ni"]
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Ni", nickel[1], 2, nickel[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Lbg"
        #
        # Molar mass
        molar_mass_pure = 2*nickel[2] + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.725, 10.118, 5.908], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 152*10**9
        # Shear modulus
        G = 66*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_olivine(self): # (Fe,Mg,Mn)2 Si O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mg", "Si", "Mn", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        y = round(rd.uniform(0, (1-x)), 2)
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], 2*y, magnesium[2]],
                                ["Si", silicon[1], 1, silicon[2]], ["Mn", manganese[1], 2*(1-x-y), manganese[2]],
                                ["Fe", iron[1], 2*x, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ca", "Ni"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ol"
        #
        # Molar mass
        molar_mass_pure = 2*(x*iron[2] + y*magnesium[2] + (1-x-y)*manganese[2]) + silicon[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[4.76, 10.2, 5.98], [], "orthorhombic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 4
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Mg = CrystalPhysics([[4.756, 10.195, 5.981], [], "orthorhombic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 4
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[4.76, 10.2, 5.98], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 4
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()
        #
        V_m = x*V_m_Fe + y*V_m_Mg + (1-x-y)*V_m_Mn
        rho = x*rho_Fe + y*rho_Mg + (1-x-y)*rho_Mn
        rho_e = x*rho_e_Fe + y*rho_e_Mg + (1-x-y)*rho_e_Mn
        # Bulk modulus
        K_Fe = 143.76*10**9
        K_Mg = 119*10**9
        K_Mn = 116*10**9
        K = x*K_Fe + y*K_Mg + (1-x-y)*K_Mn
        # Shear modulus
        G_Fe = 65.66*10**9
        G_Mg = 74*10**9
        G_Mn = 50*10**9
        G = x*G_Fe + y*G_Mg + (1-x-y)*G_Mn
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_pyrope(self): # Mg3 Al2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Mg", "Al", "Si"]
        #
        majors_data = np.array([["O", oxygen[1], 3*4, oxygen[2]], ["Mg", magnesium[1], 3, magnesium[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 3, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Ca"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Prp"
        #
        # Molar mass
        molar_mass_pure = 3*magnesium[2] + 2*aluminium[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.459], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 160*10**9
        # Shear modulus
        G = 85*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_almandine(self): # Fe3 Al2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Al", "Si", "Fe"]
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["Fe", iron[1], 3, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Alm"
        #
        # Molar mass
        molar_mass_pure = 3*iron[2] + 2*aluminium[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.526], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 164.92*10**9
        # Shear modulus
        G = 80.07*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_spessartine(self): # Mn3 Al2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Al", "Si", "Mn"]
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["Mn", manganese[1], 3, manganese[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Fe", "Mg", "Ca", "H", "Y"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Sps"
        #
        # Molar mass
        molar_mass_pure = 3*manganese[2] + 2*aluminium[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.621], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 175.88*10**9
        # Shear modulus
        G = 93.58*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_grossular(self): # Ca3 Al2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Al", "Si", "Ca"]
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["Ca", calcium[1], 3, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Cr"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Grs"
        #
        # Molar mass
        molar_mass_pure = 3*calcium[2] + 2*aluminium[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.851], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 154*10**9
        # Shear modulus
        G = 97*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_andradite(self): # Ca3 Fe2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Si", "Ca", "Fe"]
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ca", calcium[1], 3, calcium[2]], ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Cr", "Al", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Adr"
        #
        # Molar mass
        molar_mass_pure = 3*calcium[2] + 2*iron[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[12.05], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 128.57*10**9
        # Shear modulus
        G = 72.90*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_uvarovite(self): # Ca3 Cr2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        majors_name = ["O", "Si", "Ca", "Cr"]
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ca", calcium[1], 3, calcium[2]], ["Cr", chromium[1], 2, chromium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Al", "Fe", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Uv"
        #
        # Molar mass
        molar_mass_pure = 3*calcium[2] + 2*chromium[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[12.0], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 136.17*10**9
        # Shear modulus
        G = 82.58*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_aluminium_garnet(self): # (Fe,Mg,Ca,Mn)3 Al2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mg", "Al", "Si", "Ca", "Mn", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        y = round(rd.uniform(0, (1-x)), 2)
        z = round(rd.uniform(0, (1-x-y)), 2)
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Mg", magnesium[1], 3*y, magnesium[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ca", calcium[1], 3*z, calcium[2]], ["Mn", manganese[1], 3*(1-x-y-z), manganese[2]],
                                ["Fe", iron[1], 3*x, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Grt"
        #
        # Molar mass
        molar_mass_pure = 3*(x*iron[2] + y*magnesium[2] + z*calcium[2] + (1-x-y-z)*manganese[2]) + 2*aluminium[2] + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[11.526], [], "cubic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 8
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Mg = CrystalPhysics([[11.459], [], "cubic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 8
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Ca = CrystalPhysics([[11.851], [], "cubic"])
        V_Ca = dataV_Ca.calculate_volume()
        Z_Ca = 8
        V_m_Ca = MineralChemistry().calculate_molar_volume(volume_cell=V_Ca, z=Z_Ca)
        dataRho_Ca = CrystalPhysics([molar_mass, Z_Ca, V_Ca])
        rho_Ca = dataRho_Ca.calculate_bulk_density()
        rho_e_Ca = wg(amounts=amounts, elements=element, rho_b=rho_Ca).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[11.621], [], "cubic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 8
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()
        #
        V_m = x*V_m_Fe + y*V_m_Mg + z*V_m_Ca + (1-x-y-z)*V_m_Mn
        rho = x*rho_Fe + y*rho_Mg + z*rho_Ca + (1-x-y-z)*rho_Mn
        rho_e = x*rho_e_Fe + y*rho_e_Mg + z*rho_e_Ca + (1-x-y-z)*rho_e_Mn
        # Bulk modulus
        K_Fe = 164.92*10**9
        K_Mg = 160*10**9
        K_Ca = 154*10**9
        K_Mn = 175.88*10**9
        K = x*K_Fe + y*K_Mg + z*K_Ca + (1-x-y-z)*K_Mn
        # Shear modulus
        G_Fe = 80.07*10**9
        G_Mg = 85*10**9
        G_Ca = 97*10**9
        G_Mn = 93.58*10**9
        G = x*G_Fe + y*G_Mg + z*G_Ca + (1-x-y-z)*G_Mn
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_calcium_garnet(self): # Ca3 (Al,Fe,Cr)2 (SiO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Al", "Si", "Ca", "Cr", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        y = round(rd.uniform(0, (1-x)), 2)
        #
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["Al", aluminium[1], 2*x, aluminium[2]],
                                ["Si", silicon[1], 3, silicon[2]], ["Ca", calcium[1], 3, calcium[2]],
                                ["Cr", chromium[1], 2*(1-x-y), chromium[2]], ["Fe", iron[1], 2*y, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Grt"
        #
        # Molar mass
        molar_mass_pure = 3*calcium[2] + 2*(x*aluminium[2] + y*iron[2] + (1-x-y)*chromium[2]) + 3*(silicon[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Al = CrystalPhysics([[11.526], [], "cubic"])
        V_Al = dataV_Al.calculate_volume()
        Z_Al = 8
        V_m_Al = MineralChemistry().calculate_molar_volume(volume_cell=V_Al, z=Z_Al)
        dataRho_Al = CrystalPhysics([molar_mass, Z_Al, V_Al])
        rho_Al = dataRho_Al.calculate_bulk_density()
        rho_e_Al = wg(amounts=amounts, elements=element, rho_b=rho_Al).calculate_electron_density()
        #
        dataV_Fe = CrystalPhysics([[11.459], [], "cubic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 8
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Cr = CrystalPhysics([[11.851], [], "cubic"])
        V_Cr = dataV_Cr.calculate_volume()
        Z_Cr = 8
        V_m_Cr = MineralChemistry().calculate_molar_volume(volume_cell=V_Cr, z=Z_Cr)
        dataRho_Cr = CrystalPhysics([molar_mass, Z_Cr, V_Cr])
        rho_Cr = dataRho_Cr.calculate_bulk_density()
        rho_e_Cr = wg(amounts=amounts, elements=element, rho_b=rho_Cr).calculate_electron_density()
        #
        V_m = x*V_m_Al + y*V_m_Fe + (1-x-y)*V_m_Cr
        rho = x*rho_Al + y*rho_Fe + (1-x-y)*rho_Cr
        rho_e = x*rho_e_Al + y*rho_e_Fe + (1-x-y)*rho_e_Cr
        # Bulk modulus
        K_Al = 154*10**9
        K_Fe = 128.57*10**9
        K_Cr = 136.17*10**9
        K = x*K_Al + y*K_Fe + (1-x-y)*K_Cr
        # Shear modulus
        G_Al = 97*10**9
        G_Fe = 72.90*10**9
        G_Cr = 82.58*10**9
        G = x*G_Al + y*G_Fe + (1-x-y)*G_Cr
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
#
# SOROSILICATES
class Sorosilicates:
    """ Class that generates geophysical and geochemical data of sorosilicate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1):
        if self.mineral in ["Ep", "Epidote"]:
            if number > 1:
                data = [self.create_epidote() for n in range(number)]
            else:
                data = self.create_epidote()
        elif self.mineral in ["Zo", "Zoisite"]:
            if number > 1:
                data = [self.create_zoisite() for n in range(number)]
            else:
                data = self.create_zoisite()
        elif self.mineral in ["Gh", "Gehlenite"]:
            if number > 1:
                data = [self.create_gehlenite() for n in range(number)]
            else:
                data = self.create_gehlenite()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Epidote":
                data_mineral = self.create_epidote()
            elif self.mineral == "Zoisite":
                data_mineral = self.create_zoisite()
            elif self.mineral == "Gehlenite":
                data_mineral = self.create_gehlenite()
            #
            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U",
                           "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif key in ["mineral", "state", "trace elements"] and key not in dataset:
                    dataset[key] = value
                elif key in ["chemistry"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)
        #
        return dataset
    #
    def create_epidote(self): # Ca2 Al2 Fe Si3 H O13
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Al", "Si", "Ca", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 13, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ca", calcium[1], 2, calcium[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mg", "Mn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Ep"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + 2*aluminium[2] + iron[2] + 3*silicon[2] + hydrogen[2] + 13*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.98, 5.64, 10.22], [115.4], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 138.35*10**9
        # Shear modulus
        G = 78.35*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_zoisite(self): # Ca2 Al3 Si3 O12 (OH)
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "O", "Al", "Si", "Ca"]
        #
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 13, oxygen[2]],
                                ["Al", aluminium[1], 3, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ca", calcium[1], 2, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Mg", "Cr", "Ti", "Ca", "Na", "V", "Sr"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Zo"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + 3*aluminium[2] + 3*silicon[2] + 12*oxygen[2] + (hydrogen[2]+oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[16.24, 5.58, 10.1], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 145.56*10**9
        # Shear modulus
        G = 89.97*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_gehlenite(self): # Ca2 Al2 Si O7
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Al", "Si", "Ca"]
        #
        majors_data = np.array([["O", oxygen[1], 7, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
                                ["Si", silicon[1], 1, silicon[2]], ["Ca", calcium[1], 2, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Fe", "Mg", "Mn", "Na", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Gh"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + 2*aluminium[2] + silicon[2] + 7*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.69, 5.067], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 102*10**9
        # Shear modulus
        G = 53*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
#
# INOSILICATES
class Inosilicates:
    """ Class that generates geophysical and geochemical data of inosilicate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def get_data(self, number=1): # ["Enstatite", "Diopside", "Augite", "Jadeite", "Aegirine", "Actinolite",
        # "Arfvedsonite", "Ca-Amphibole", "Ca-Pyroxene", "Ferrosilite", "Glaucophane", "Mg-Fe-Pyroxene", "Riebeckite",
        # "Na-Amphibole", "Spodumene", "Tremolite", "Wollastonite"]
        if self.mineral in ["En", "Enstatite"]:
            if number > 1:
                data = [self.create_enstatite() for n in range(number)]
            else:
                data = self.create_enstatite()
        elif self.mineral in ["Di", "Diopside"]:
            if number > 1:
                data = [self.create_diopside() for n in range(number)]
            else:
                data = self.create_diopside()
        elif self.mineral in ["Aug", "Augite"]:
            if number > 1:
                data = [self.create_augite() for n in range(number)]
            else:
                data = self.create_augite()
        elif self.mineral in ["Jd", "Jadeite"]:
            if number > 1:
                data = [self.create_jadeite() for n in range(number)]
            else:
                data = self.create_jadeite()
        elif self.mineral in ["Aeg", "Aegirine"]:
            if number > 1:
                data = [self.create_aegirine() for n in range(number)]
            else:
                data = self.create_aegirine()
        elif self.mineral in ["Act", "Actinolite"]:
            if number > 1:
                data = [self.create_actinolite() for n in range(number)]
            else:
                data = self.create_actinolite()
        elif self.mineral in ["Arf", "Arfvedsonite"]:
            if number > 1:
                data = [self.create_arfvedsonite() for n in range(number)]
            else:
                data = self.create_arfvedsonite()
        elif self.mineral in ["Amph", "Ca-Amphibole"]:
            if number > 1:
                data = [self.create_calcium_amphibole() for n in range(number)]
            else:
                data = self.create_calcium_amphibole()
        elif self.mineral in ["Px", "Ca-Pyroxene"]:
            if number > 1:
                data = [self.create_calium_pyroxene() for n in range(number)]
            else:
                data = self.create_calium_pyroxene()
        elif self.mineral in ["Fs", "Ferrosilite"]:
            if number > 1:
                data = [self.create_ferrosilite() for n in range(number)]
            else:
                data = self.create_ferrosilite()
        elif self.mineral in ["Dnp", "Donpeacorite"]:
            if number > 1:
                data = [self.create_donpeacorite() for n in range(number)]
            else:
                data = self.create_donpeacorite()
        elif self.mineral in ["Gln", "Glaucophane"]:
            if number > 1:
                data = [self.create_glaucophane() for n in range(number)]
            else:
                data = self.create_glaucophane()
        elif self.mineral in ["Px", "Mg-Fe-Pyroxene"]:
            if number > 1:
                data = [self.create_mg_fe_pyroxene() for n in range(number)]
            else:
                data = self.create_mg_fe_pyroxene()
        elif self.mineral in ["Rbk", "Riebeckite"]:
            if number > 1:
                data = [self.create_riebeckite() for n in range(number)]
            else:
                data = self.create_riebeckite()
        elif self.mineral in ["Amph", "Na-Amphibole"]:
            if number > 1:
                data = [self.create_sodium_amphibole() for n in range(number)]
            else:
                data = self.create_sodium_amphibole()
        elif self.mineral in ["Spd", "Spodumene"]:
            if number > 1:
                data = [self.create_spodumene() for n in range(number)]
            else:
                data = self.create_spodumene()
        elif self.mineral in ["Tr", "Tremolite"]:
            if number > 1:
                data = [self.create_tremolite() for n in range(number)]
            else:
                data = self.create_tremolite()
        elif self.mineral in ["Wo", "Wollastonite"]:
            if number > 1:
                data = [self.create_wollastonite() for n in range(number)]
            else:
                data = self.create_wollastonite()
        #
        return data
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Enstatite":
                data_mineral = self.create_enstatite()
            elif self.mineral == "Ferrosilite":
                data_mineral = self.create_ferrosilite()
            elif self.mineral == "Diopside":
                data_mineral = self.create_diopside()
            #
            # Na-Pyroxenes
            elif self.mineral == "Jadeite":
                data_mineral = self.create_jadeite()
            elif self.mineral == "Aegirine":
                data_mineral = self.create_aegirine()
            elif self.mineral == "Kosmochlor":
                data_mineral = self.create_kosmochlor()
            #
            elif self.mineral == "Spodumene":
                data_mineral = self.create_spodumene()
            elif self.mineral == "Wollastonite":
                data_mineral = self.create_wollastonite()
            elif self.mineral == "Tremolite":
                data_mineral = self.create_tremolite()
            elif self.mineral == "Actinolite":
                data_mineral = self.create_actinolite()
            elif self.mineral == "Glaucophane":
                data_mineral = self.create_glaucophane()
            elif self.mineral == "Augite":
                data_mineral = self.create_augite()
            elif self.mineral == "Riebeckite":
                data_mineral = self.create_riebeckite()
            elif self.mineral == "Arfvedsonite":
                data_mineral = self.create_arfvedsonite()
            elif self.mineral == "Ca-Amphibole":
                data_mineral = self.create_calcium_amphibole()
            elif self.mineral == "Na-Amphibole":
                data_mineral = self.create_sodium_amphibole()
            elif self.mineral == "Mg-Fe-Pyroxene":
                data_mineral = self.create_mg_fe_pyroxene()
            elif self.mineral == "Ca-Pyroxene":
                data_mineral = self.create_calium_pyroxene()
            elif self.mineral == "Na-Pyroxene":
                data_mineral = self.create_sodium_pyroxene()
            elif self.mineral == "Donpeacorite":
                data_mineral = self.create_donpeacorite()
            elif self.mineral == "Orthopyroxene":
                data_mineral = self.create_orthopyroxene()
            elif self.mineral == "Clinopyroxene":
                data_mineral = self.create_clinopyroxene()
            #
            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U",
                           "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif key in ["mineral", "state", "trace elements"] and key not in dataset:
                    dataset[key] = value
                elif key in ["chemistry"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)
        #
        return dataset
    #
    def create_enstatite(self): # Mg (Mg,Fe) Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mg", "Si", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], 1+x, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Fe", iron[1], (1-x), iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ca", "Al", "Co", "Ni", "Mn", "Ti", "Cr", "Na", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "En"
        #
        # Molar mass
        molar_mass_pure = magnesium[2] + x*magnesium[2] + (1-x)*iron[2] + 2*silicon[2] +6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[18.228, 8.805, 5.185], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K_Mg = 136.17*10**9
        K_Fe = 146.44*10**9
        K = x*K_Mg + (1-x)*K_Fe
        # Shear modulus
        G_Mg = 80.79*10**9
        G_Fe = 78.29*10**9
        G = x*G_Mg + (1-x)*G_Fe
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_ferrosilite(self): # Fe Si O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Si", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Si", silicon[1], 2, silicon[2]],
                                ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ca", "Na", "K", "Al", "Co", "Ni", "Mn", "Ti", "Cr"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Fs"
        #
        # Molar mass
        molar_mass_pure = 2*iron[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[18.418, 9.078, 5.237], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 141.77*10**9
        # Shear modulus
        G = 70.44*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_donpeacorite(self):  # (Mg,Mn)  Mg Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Mg", "Si", "Mn"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], 1 + x, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Mn", manganese[1], (1 - x), manganese[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Al", "Fe", "Ca", "Na"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Dnp"
        #
        # Molar mass
        molar_mass_pure = (1 + x)*magnesium[2] + (1 - x)*manganese[2] +  2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Mg = CrystalPhysics([[18.228, 8.805, 5.185], [], "orthorhombic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 8
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[18.834, 8.878, 5.226], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 8
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()
        #
        V_m = x*V_m_Mg + (1 - x)*V_m_Mn
        rho = x*rho_Mg + (1 - x)*rho_Mn
        rho_e = x*rho_e_Mg + (1 - x)*rho_e_Mn
        # Bulk modulus
        K_Mg = 138.51*10**9
        K_Mn = 144.5*10**9
        K = x*K_Mg + (1 - x)*K_Mn
        # Shear modulus
        G_Mg = 81.89*10**9
        G_Mn = 78.86*10**9
        G = x*G_Mg + (1 - x)*G_Mn
        # Young's modulus
        E = (9 * K * G) / (3 * K + G)
        # Poisson's ratio
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
        # vP/vS
        vPvS = ((K + 4 / 3 * G) / G) ** 0.5
        # P-wave velocity
        vP = ((K + 4 / 3 * G) / rho) ** 0.5
        # S-wave velocity
        vS = (G / rho) ** 0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe * rho_e * 10 ** (-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append(
                [round(K * 10 ** (-9), 2), round(G * 10 ** (-9), 2), round(E * 10 ** (-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G * 10 ** (-9), 4)
            results["K"] = round(K * 10 ** (-9), 4)
            results["E"] = round(E * 10 ** (-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_orthopyroxene(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mg", "Si", "Mn", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        y = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], x*y, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Mn", manganese[1], x*(1 - y), manganese[2]],
                                ["Fe", iron[1], (1 - x)*2, iron[2]]], dtype=object)
        majors_data_Mg = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], 2, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Mn", manganese[1], 0, manganese[2]],
                                ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_Mn = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Mn", manganese[1], 1, manganese[2]],
                                ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_Fe = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], 0, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Mn", manganese[1], 0, manganese[2]],
                                ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Al", "Ca", "Na", "K", "Co", "Ni", "Ti", "Cr"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Opx"
        #
        # Molar mass
        molar_mass_pure = 2*(x*(y*magnesium[2] + (1 - y)*manganese[2]) + (1 - x)*iron[2]) + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        molar_mass_pure_Mg = 2*magnesium[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass_pure_Mn = magnesium[2] + manganese[2] + 2 * silicon[2] + 6 * oxygen[2]
        molar_mass_pure_Fe = 2 * iron[2] + 2 * silicon[2] + 6 * oxygen[2]
        molar_mass_Mg, amounts_Mg = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Mg, majors=majors_data_Mg).calculate_molar_mass()
        element_Mg = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        molar_mass_Mn, amounts_Mn = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Mn, majors=majors_data_Mn).calculate_molar_mass()
        element_Mn = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        molar_mass_Fe, amounts_Fe = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Fe, majors=majors_data_Fe).calculate_molar_mass()
        element_Fe = [PeriodicSystem(name=amounts_Fe[i][0]).get_data() for i in range(len(amounts_Fe))]
        # Density
        dataV_Mg = CrystalPhysics([[18.228, 8.805, 5.185], [], "orthorhombic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 8
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass_Mg, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts_Mg, elements=element_Mg, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[18.834, 8.878, 5.226], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 8
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass_Mn, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts_Mn, elements=element_Mn, rho_b=rho_Mn).calculate_electron_density()
        #
        dataV_Fe = CrystalPhysics([[18.418, 9.078, 5.237], [], "orthorhombic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 8
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass_Fe, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts_Fe, elements=element_Fe, rho_b=rho_Fe).calculate_electron_density()
        #
        V_m = x*(y*V_m_Mg + (1 - y)*V_m_Mn) + (1 - x)*V_m_Fe
        rho = x*(y*rho_Mg + (1 - y)*rho_Mn) + (1 - x)*rho_Fe
        rho_e = x*(y*rho_e_Mg + (1 - y)*rho_e_Mn) + (1 - x)*rho_e_Fe
        # Bulk modulus
        K_Mg = 138.51 * 10 ** 9
        K_Mn = 144.5 * 10 ** 9
        K_Fe = 141.77 * 10 ** 9
        K = x*(y*K_Mg + (1 - y)*K_Mn) + (1 - x)*K_Fe
        # Shear modulus
        G_Mg = 81.89 * 10 ** 9
        G_Mn = 78.86 * 10 ** 9
        G_Fe = 70.44 * 10 ** 9
        G = x*(y*G_Mg + (1 - y)*G_Mn) + (1 - x)*G_Fe
        # Young's modulus
        E = (9 * K * G) / (3 * K + G)
        # Poisson's ratio
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
        # vP/vS
        vPvS = ((K + 4 / 3 * G) / G) ** 0.5
        # P-wave velocity
        vP = ((K + 4 / 3 * G) / rho) ** 0.5
        # S-wave velocity
        vS = (G / rho) ** 0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe * rho_e * 10 ** (-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append(
                [round(K * 10 ** (-9), 2), round(G * 10 ** (-9), 2), round(E * 10 ** (-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G * 10 ** (-9), 4)
            results["K"] = round(K * 10 ** (-9), 4)
            results["E"] = round(E * 10 ** (-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_sodium_pyroxene(self, enrichment=None):
        mineralogy = {"x": self.create_jadeite(), "y": self.create_kosmochlor(), "z": self.create_aegirine()}
        #
        list_elements = []
        chemistry = {}
        for key, dataset in mineralogy.items():
            for element, value in dataset["chemistry"].items():
                if element not in list_elements:
                    list_elements.append(element)
                    chemistry[element] = 0
        #
        old_index = list_elements.index("O")
        list_elements += [list_elements.pop(old_index)]
        #
        x = round(rd.uniform(0, 1), 4)
        y = round(rd.uniform(0, 1 - x), 4)
        z = round(1 - x - y, 4)
        amounts_components = {"x": x, "y": y, "z": z}
        #
        mineral = "Na-Cpx"
        var_state = "variable"
        #
        molar_mass = 0
        rho = 0
        rho_e = 0
        V_m = 0
        K = 0
        G = 0
        #
        w_total = 0
        for key, value in mineralogy.items():
            molar_mass += amounts_components[key]*mineralogy[key]["M"]
            rho += amounts_components[key]*mineralogy[key]["rho"]
            rho_e += amounts_components[key]*mineralogy[key]["rho_e"]
            V_m += amounts_components[key]*mineralogy[key]["V"]
            K += amounts_components[key]*mineralogy[key]["K"]
            G += amounts_components[key]*mineralogy[key]["G"]
            #
            for element in list_elements:
                if element in mineralogy[key]["chemistry"]:
                    if element != "O":
                        value = amounts_components[key]*mineralogy[key]["chemistry"][element]
                        chemistry[element] += value
                        w_total += value
        #
        chemistry["O"] = 1 - round(w_total, 6)
        for key, value in chemistry.items():
            chemistry[key] = round(value, 6)
        #
        # Molar mass
        molar_mass = molar_mass
        # Density
        rho = rho
        # Electron density
        rho_e = rho_e
        # Molar volume
        V_m = V_m
        # Bulk modulus
        K = K
        # Shear modulus
        G = G
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        amounts = []
        for element in list_elements:
            amounts.append([element, PeriodicSystem(name="O").get_data()[0], chemistry[element]])
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G, 4)
        results["K"] = round(K, 4)
        results["E"] = round(E, 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_clinopyroxene(self, enrichment=None):
        mineralogy = {"x": self.create_calium_pyroxene(), "y": self.create_sodium_pyroxene()}
        #
        list_elements = []
        chemistry = {}
        for key, dataset in mineralogy.items():
            for element, value in dataset["chemistry"].items():
                if element not in list_elements:
                    list_elements.append(element)
                    chemistry[element] = 0
        #
        old_index = list_elements.index("O")
        list_elements += [list_elements.pop(old_index)]
        #
        x = round(rd.uniform(0, 1), 4)
        y = round(1 - x, 4)
        amounts_components = {"x": x, "y": y}
        #
        mineral = "Cpx"
        var_state = "variable"
        #
        molar_mass = 0
        rho = 0
        rho_e = 0
        V_m = 0
        K = 0
        G = 0
        #
        w_total = 0
        for key, value in mineralogy.items():
            molar_mass += amounts_components[key]*mineralogy[key]["M"]
            rho += amounts_components[key]*mineralogy[key]["rho"]
            rho_e += amounts_components[key]*mineralogy[key]["rho_e"]
            V_m += amounts_components[key]*mineralogy[key]["V"]
            K += amounts_components[key]*mineralogy[key]["K"]
            G += amounts_components[key]*mineralogy[key]["G"]
            #
            for element in list_elements:
                if element in mineralogy[key]["chemistry"]:
                    if element != "O":
                        value = amounts_components[key]*mineralogy[key]["chemistry"][element]
                        chemistry[element] += value
                        w_total += value
        #
        chemistry["O"] = 1 - round(w_total, 6)
        for key, value in chemistry.items():
            chemistry[key] = round(value, 6)
        #
        # Molar mass
        molar_mass = molar_mass
        # Density
        rho = rho
        # Electron density
        rho_e = rho_e
        # Molar volume
        V_m = V_m
        # Bulk modulus
        K = K
        # Shear modulus
        G = G
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        amounts = []
        for element in list_elements:
            amounts.append([element, PeriodicSystem(name="O").get_data()[0], chemistry[element]])
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G, 4)
        results["K"] = round(K, 4)
        results["E"] = round(E, 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_diopside(self): # Ca Mg Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Mg", "Si", "Ca"]
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "V", "Cr", "Mn", "Zn", "Al", "Ti", "Na", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Di"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + magnesium[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.761, 8.926, 5.258], [105.8], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 102*10**9
        # Shear modulus
        G = 64*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_jadeite(self): # Na (Al,Fe) Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Na", "Al", "Si", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Na", sodium[1], 1, sodium[2]],
                                ["Al", aluminium[1], x, aluminium[2]], ["Si", silicon[1], 2, silicon[2]],
                                ["Fe", iron[1], (1-x), iron[2]]], dtype=object)
        majors_data_al = np.array([["O", oxygen[1], 6, oxygen[2]], ["Na", sodium[1], 1, sodium[2]],
                                ["Al", aluminium[1], 1, aluminium[2]], ["Si", silicon[1], 2, silicon[2]],
                                ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_fe = np.array([["O", oxygen[1], 6, oxygen[2]], ["Na", sodium[1], 1, sodium[2]],
                                ["Al", aluminium[1], 0, aluminium[2]], ["Si", silicon[1], 2, silicon[2]],
                                ["Fe", iron[1], 1, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Mg", "Ca", "K", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Jd"
        #
        # Molar mass
        molar_mass_pure = sodium[2] + (x*aluminium[2] + (1-x)*iron[2]) + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        molar_mass_pure_al = sodium[2] + aluminium[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass_al, amounts_al = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_al, majors=majors_data_al).calculate_molar_mass()
        element_al = [PeriodicSystem(name=amounts_al[i][0]).get_data() for i in range(len(amounts_al))]
        #
        molar_mass_pure_fe = sodium[2] + iron[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass_fe, amounts_fe = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_fe, majors=majors_data_fe).calculate_molar_mass()
        element_fe = [PeriodicSystem(name=amounts_fe[i][0]).get_data() for i in range(len(amounts_fe))]
        # Density
        dataV = CrystalPhysics([[9.418, 8.562, 5.219], [107.56], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        dataV_al = CrystalPhysics([[9.418, 8.562, 5.219], [107.56], "monoclinic"])
        V_al = dataV_al.calculate_volume()
        Z_al = 4
        V_m_al = MineralChemistry().calculate_molar_volume(volume_cell=V_al, z=Z_al)
        dataRho_al = CrystalPhysics([molar_mass_al, Z_al, V_al])
        rho_al = dataRho_al.calculate_bulk_density()
        rho_e_al = wg(amounts=amounts_al, elements=element_al, rho_b=rho_al).calculate_electron_density()
        #
        dataV_fe = CrystalPhysics([[9.65, 8.79, 5.29], [107.5], "monoclinic"])
        V_fe = dataV_fe.calculate_volume()
        Z_fe = 4
        V_m_fe = MineralChemistry().calculate_molar_volume(volume_cell=V_fe, z=Z_fe)
        dataRho_fe = CrystalPhysics([molar_mass_fe, Z_fe, V_fe])
        rho_fe = dataRho_fe.calculate_bulk_density()
        rho_e_fe = wg(amounts=amounts_fe, elements=element_fe, rho_b=rho_fe).calculate_electron_density()
        #
        V_m = x*V_m_al + (1-x)*V_m_fe
        rho = x*rho_al + (1-x)*rho_fe
        rho_e = x*rho_e_al + (1-x)*rho_e_fe
        # Bulk modulus
        K_Al = 74.38*10**9
        K_Fe = 124.44*10**9
        K = x*K_Al + (1-x)*K_Fe
        # Shear modulus
        G_Al = 49.63*10**9
        G_Fe = 74.16*10**9
        G = x*G_Al + (1-x)*G_Fe
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_aegirine(self): # Na Fe Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Na", "Si", "Fe"]
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Na", sodium[1], 1, sodium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Al", "Ti", "V", "Mn", "Mg", "Ca", "K", "Zr", "Ce"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Aeg"
        #
        # Molar mass
        molar_mass_pure = sodium[2] + iron[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.65, 8.79, 5.29], [107.5], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 124.44*10**9
        # Shear modulus
        G = 74.16*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_kosmochlor(self):  # Na Cr Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        chromium = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Na", "Si", "Cr"]
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Na", sodium[1], 1, sodium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Cr", chromium[1], 1, chromium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Al", "Fe", "Mn", "Mg", "Ca", "K", "P"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Kos"
        #
        # Molar mass
        molar_mass_pure = sodium[2] + chromium[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.574, 8.712, 5.265], [107.49], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 115.0*10**9
        # Shear modulus
        G = 72.0*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_spodumene(self): # Li Al Si2 O6
        # Major elements
        lithium = PeriodicSystem(name="Li").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["Li", "O", "Al", "Si"]
        #
        majors_data = np.array([["Li", lithium[1], 1, lithium[2]], ["O", oxygen[1], 6, oxygen[2]],
                                ["Al", aluminium[1], 1, aluminium[2]], ["Si", silicon[1], 2, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Mg", "Ca", "Na", "K", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Spd"
        #
        # Molar mass
        molar_mass_pure = lithium[2] + aluminium[2] + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.46, 8.39, 5.22], [110.17], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 144.81*10**9
        # Shear modulus
        G = 103.66*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_wollastonite(self): # Ca Si O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Si", "Ca"]
        #
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Si", silicon[1], 1, silicon[2]],
                                ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = [None]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Wo"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + silicon[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.94, 7.32, 7.07], [90.033, 95.367, 103.433], "triclinic"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 93*10**9
        # Shear modulus
        G = 47*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_tremolite(self): # Ca2 Mg5 (Si8 O22) (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "O", "Mg", "Si", "Ca"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Mg", magnesium[1], 5, magnesium[2]], ["Si", silicon[1], 8, silicon[2]],
                                ["Ca", calcium[1], 2, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Al", "Na", "K", "F", "Cl", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Tr"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + 5*magnesium[2] + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.8385, 18.0554, 5.2778], [104.751], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 124.09*10**9
        # Shear modulus
        G = 72.90*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_actinolite(self): # Ca2 (Mg,Fe)5 (Si8 O22) (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Mg", "Si", "Ca", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Mg", magnesium[1], 5*x, magnesium[2]], ["Si", silicon[1], 8, silicon[2]],
                                ["Ca", calcium[1], 2, calcium[2]], ["Fe", iron[1], 5*(1-x), iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Al", "Na", "K", "Ti"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Act"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + 5*(x*magnesium[2] + (1-x)*iron[2]) + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.84, 18.1, 5.28], [104.7], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K_Mg = 124.09*10**9
        K_Fe = 85*10**9
        K = x*K_Mg + (1-x)*K_Fe
        # Shear modulus
        G_Mg = 72.90*10**9
        G_Fe = 47*10**9
        G = x*G_Mg + (1-x)*G_Fe
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_glaucophane(self): # Na2 Mg3 Al2 (Si8 O22) (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["H", "O", "Na", "Mg", "Al", "Si"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Na", sodium[1], 2, sodium[2]], ["Mg", magnesium[1], 3, magnesium[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 8, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Li", "Ti", "Cr", "Mn", "Ca", "K", "F", "Cl"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Gln"
        #
        # Molar mass
        molar_mass_pure = 2*sodium[2] + 3*magnesium[2] + 2*aluminium[2] + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.541, 17.74, 5.295], [103.67], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 96*10**9
        # Shear modulus
        G = 64*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_augite(self):   # (Ca,Mg,Fe) (Mg,Fe) Si2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mg", "Si", "Ca", "Fe"]
        #
        x_Ca = round(rd.uniform(0.4, 0.9), 2)
        x_Mg1 = round(rd.uniform(0.0, 1-x_Ca), 2)
        x_Fe1 = 1-x_Ca-x_Mg1
        x_Mg2 = round(rd.uniform(0.0, 1.0), 2)
        x_Fe2 = 1-x_Mg2
        #
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], x_Mg1+x_Mg2, magnesium[2]],
                                ["Si", silicon[1], 2, silicon[2]], ["Ca", calcium[1], x_Ca, calcium[2]],
                                ["Fe", iron[1], x_Fe1+x_Fe2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Cr", "Na", "Mn", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Aug"
        #
        # Molar mass
        molar_mass_pure = (x_Ca*calcium[2]+x_Mg1*magnesium[2]+x_Fe1*iron[2]) + (x_Mg2*magnesium[2]+x_Fe2*iron[2]) + 2*silicon[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.8, 9.0, 5.25], [105], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = (95.72 + (106.97-95.72)/(3420-3320)*(rho-3320))*10**9
        # Shear modulus
        G = (58.01 + (57.21-58.01)/(3420-3320)*(rho-3320))*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_riebeckite(self): # Na2 Fe5 (Si8 O22) (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Na", "Si", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Na", sodium[1], 2, sodium[2]], ["Si", silicon[1], 8, silicon[2]],
                                ["Fe", iron[1], 5, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Ti", "Mg", "Al", "Mn"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Rbk"
        #
        # Molar mass
        molar_mass_pure = 2*sodium[2] + 5*iron[2] + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.769, 18.048, 5.335], [103.6], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 101.1*10**9
        # Shear modulus
        G = 43.7*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_arfvedsonite(self): # Na3 Fe5 (Si8 O22) (OH)2
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Na", "Si", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Na", sodium[1], 3, sodium[2]], ["Si", silicon[1], 8, silicon[2]],
                                ["Fe", iron[1], 5, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Ca", "Al", "K", "F"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Arf"
        #
        # Molar mass
        molar_mass_pure = 3*sodium[2] + 5*iron[2] + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.9, 18.0, 5.3], [104], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 101.1*10**9
        # Shear modulus
        G = 43.7*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_calcium_amphibole(self): # Tr + Act
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Mg", "Si", "Ca", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Mg", magnesium[1], 5*x, magnesium[2]], ["Si", silicon[1], 8, silicon[2]],
                                ["Ca", calcium[1], 2, calcium[2]], ["Fe", iron[1], 5*(1-x), iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Al", "Na", "K", "Ti", "F", "Cl", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Amph"
        #
        # Molar mass
        molar_mass_pure = 2*calcium[2] + 5*(x*magnesium[2] + (1-x)*iron[2]) + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.84, 18.1, 5.28], [104.7], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K_Mg = 124.09*10**9
        K_Fe = 85*10**9
        K = x*K_Mg + (1-x)*K_Fe
        # Shear modulus
        G_Mg = 72.90*10**9
        G_Fe = 47*10**9
        G = x*G_Mg + (1-x)*G_Fe
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_sodium_amphibole(self): # Rbk + Arf
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Na", "Si", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 22+2, oxygen[2]],
                                ["Na", sodium[1], 2+x, sodium[2]], ["Si", silicon[1], 8, silicon[2]],
                                ["Fe", iron[1], 5, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Ca", "Al", "K", "F", "Fe", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Amph"
        #
        # Molar mass
        molar_mass_pure = (2+x)*sodium[2] + 5*iron[2] + 8*silicon[2] + 22*oxygen[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Rbk = CrystalPhysics([[9.769, 18.048, 5.335], [103.6], "monoclinic"])
        V_Rbk = dataV_Rbk.calculate_volume()
        Z_Rbk = 2
        V_m_Rbk = MineralChemistry().calculate_molar_volume(volume_cell=V_Rbk, z=Z_Rbk)
        dataRho_Rbk = CrystalPhysics([molar_mass, Z_Rbk, V_Rbk])
        rho_Rbk = dataRho_Rbk.calculate_bulk_density()
        rho_e_Rbk = wg(amounts=amounts, elements=element, rho_b=rho_Rbk).calculate_electron_density()
        #
        dataV_Arf = CrystalPhysics([[9.9, 18.0, 5.3], [104], "monoclinic"])
        V_Arf = dataV_Arf.calculate_volume()
        Z_Arf = 2
        V_m_Arf = MineralChemistry().calculate_molar_volume(volume_cell=V_Arf, z=Z_Arf)
        dataRho_Arf = CrystalPhysics([molar_mass, Z_Arf, V_Arf])
        rho_Arf = dataRho_Arf.calculate_bulk_density()
        rho_e_Arf = wg(amounts=amounts, elements=element, rho_b=rho_Arf).calculate_electron_density()
        #
        V_m = x*V_m_Rbk + (1-x)*V_m_Arf
        rho = x*rho_Rbk + (1-x)*rho_Arf
        rho_e = x*rho_e_Rbk + (1-x)*rho_e_Arf
        # Bulk modulus
        K = 101.1*10**9
        # Shear modulus
        G = 43.7*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            element_list = np.array(amounts)[:, 0]
            results["chemistry"] = {}
            for index, element in enumerate(element_list, start=0):
                results["chemistry"][element] = amounts[index][2]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_mg_fe_pyroxene(self): # En + Fs
        # Minerals
        enstatite = Inosilicates(data_type=True).create_enstatite()
        ferrosilite = Inosilicates(data_type=True).create_ferrosilite()
        majors_name = ["O", "Mg", "Si", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Ca", "Al", "K", "F", "Fe", "Mg"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Px"
        #
        # Molar mass
        molar_mass = x*enstatite["M"] + (1-x)*ferrosilite["M"]
        V_m = x*enstatite["V"] + (1-x)*ferrosilite["V"]
        amounts = {}
        for element in majors_name:
            if element != "Mg":
                amounts[element] = round(x*enstatite["chemistry"][element] + (1-x)*ferrosilite["chemistry"][element], 6)
            elif element == "Mg":
                amounts[element] = round(x*enstatite["chemistry"][element], 6)
        sum_w = 0
        for element in majors_name:
            sum_w += amounts[element]
        amounts["O"] = round(amounts["O"]-(sum_w-1), 6)
        # Density
        rho = x*enstatite["rho"] + (1-x)*ferrosilite["rho"]
        rho_e = x*enstatite["rho_e"] + (1-x)*ferrosilite["rho_e"]
        # Bulk modulus
        K = (x*enstatite["K"] + (1-x)*ferrosilite["K"])*10**9
        # Shear modulus
        G = (x*enstatite["G"] + (1-x)*ferrosilite["G"])*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = x*enstatite["GR"] + (1-x)*ferrosilite["GR"]
        # Photoelectricity
        pe = x*enstatite["PE"] + (1-x)*ferrosilite["PE"]
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            results["chemistry"] = {}
            for index, element in enumerate(majors_name, start=0):
                results["chemistry"][element] = amounts[element]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
    #
    def create_calium_pyroxene(self): # Aug + Di
        # Minerals
        augite = Inosilicates(data_type=True).create_augite()
        diopside = Inosilicates(data_type=True).create_diopside()
        majors_name = ["O", "Mg", "Si", "Ca", "Fe"]
        #
        x = round(rd.uniform(0, 1), 2)
        #
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "variable"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Cr", "Na", "Mn", "K", "Fe", "V", "Zn", "Al"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Px"
        #
        # Molar mass
        molar_mass = x*augite["M"] + (1-x)*diopside["M"]
        V_m = x*augite["V"] + (1-x)*diopside["V"]
        amounts = {}
        for element in majors_name:
            if element != "Fe":
                amounts[element] = round(x*augite["chemistry"][element] + (1-x)*diopside["chemistry"][element], 6)
            elif element == "Fe":
                amounts[element] = round(x*augite["chemistry"][element], 6)
        sum_w = 0
        for element in majors_name:
            sum_w += amounts[element]
        amounts["O"] = round(amounts["O"]-(sum_w-1), 6)
        # Density
        rho = x*augite["rho"] + (1-x)*diopside["rho"]
        rho_e = x*augite["rho_e"] + (1-x)*diopside["rho_e"]
        # Bulk modulus
        K = (x*augite["K"] + (1-x)*diopside["K"])*10**9
        # Shear modulus
        G = (x*augite["G"] + (1-x)*diopside["G"])*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = x*augite["GR"] + (1-x)*diopside["GR"]
        # Photoelectricity
        pe = x*augite["PE"] + (1-x)*diopside["PE"]
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        if self.data_type == False:
            data = []
            data.append(mineral)
            data.append(round(molar_mass, 3))
            data.append(round(rho, 2))
            data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
            data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
            data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
            data.append(amounts)
            #
            return data
        else:
            #
            results = {}
            results["mineral"] = mineral
            results["M"] = round(molar_mass, 3)
            results["state"] = var_state
            results["chemistry"] = {}
            for index, element in enumerate(majors_name, start=0):
                results["chemistry"][element] = amounts[element]
            results["rho"] = round(rho, 4)
            results["rho_e"] = round(rho_e, 4)
            results["V"] = round(V_m, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vPvS, 4)
            results["G"] = round(G*10**(-9), 4)
            results["K"] = round(K*10**(-9), 4)
            results["E"] = round(E*10**(-9), 4)
            results["nu"] = round(nu, 4)
            results["GR"] = round(gamma_ray, 4)
            results["PE"] = round(pe, 4)
            results["U"] = round(U, 4)
            if p != None:
                results["p"] = round(p, 4)
            else:
                results["p"] = p
            #
            return results
#
# CYCLOSILICATES
class Cyclosilicates:
    """ Class that generates geophysical and geochemical data of cyclosilicate minerals"""
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral

        # Chemistry
        self.hydrogen = ["H", 1, 1.008]
        self.lithium = ["Li", 3, 6.938]
        self.boron = ["B", 5, 10.806]
        self.carbon = ["C", 6, 12.011]
        self.oxygen = ["O", 8, 15.999]
        self.flourine = ["F", 9, 18.998]
        self.sodium = ["Na", 11, 22.990]
        self.magnesium = ["Mg", 12, 24.304]
        self.aluminium = ["Al", 13, 26.982]
        self.silicon = ["Si", 14, 28.085]
        self.chlorine = ["Cl", 17, 35.450]
        self.potassium = ["K", 19, 39.098]
        self.calcium = ["Ca", 20, 40.078]
        self.iron = ["Fe", 26, 55.845]

    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Beryl":
                data_mineral = self.create_beryl()
            elif self.mineral == "Benitoite":
                data_mineral = self.create_benitoite()
            elif self.mineral == "Cordierite":
                data_mineral = self.create_cordierite()
            elif self.mineral == "Sekaninaite":
                data_mineral = self.create_sekaninaite()
            elif self.mineral == "Schorl":
                data_mineral = self.create_schorl()
            elif self.mineral == "Elbaite":
                data_mineral = self.create_elbaite()
            elif self.mineral == "Liddicoatite":
                data_mineral = self.create_liddicoatite()
            #
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
        #
        return dataset
    #
    def create_beryl(self): # Be3 Al2 (Si6 O18)
        # Major elements
        beryllium = PeriodicSystem(name="Be").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["Be", "O", "Al", "Si"]
        #
        majors_data = np.array([["Be", beryllium[1], 3, beryllium[2]], ["O", oxygen[1], 18, oxygen[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Si", silicon[1], 6, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Mg", "Ca", "Cr", "Na", "Li", "Cs", "H", "K", "Rb"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Brl"
        #
        # Molar mass
        molar_mass_pure = 3*beryllium[2] + 2*aluminium[2] + 6*silicon[2] + 18*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.215, 9.192], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 131.29*10**9
        # Shear modulus
        G = 84.61*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_benitoite(self): # Ba Ti Si3 O9
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        barium = PeriodicSystem(name="Ba").get_data()
        majors_name = ["O", "Si", "Ti", "Ba"]
        #
        majors_data = np.array([["O", oxygen[1], 9, oxygen[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["Ti", titanium[1], 1, titanium[2]], ["Ba", barium[1], 1, barium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Na"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Bnt"
        #
        # Molar mass
        molar_mass_pure = barium[2] + titanium[2] + 3*silicon[2] + 9*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[6.6, 9.71], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 116.01*10**9
        # Shear modulus
        G = 71.49*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_cordierite(self): # Mg2 Al4 Si5 O18
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Mg", "Al", "Si"]
        #
        majors_data = np.array([["O", oxygen[1], 18, oxygen[2]], ["Mg", magnesium[1], 2, magnesium[2]],
                                ["Al", aluminium[1], 4, aluminium[2]], ["Si", silicon[1], 5, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Fe", "Ti", "Ca", "Na", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Crd"
        #
        # Molar mass
        molar_mass_pure = 2*magnesium[2] + 4*aluminium[2] + 5*silicon[2] + 18*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[17.13, 9.8, 9.35], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 122*10**9
        # Shear modulus
        G = 50*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_sekaninaite(self): # Fe2 Al4 Si5 O18
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Al", "Si", "Fe"]
        #
        majors_data = np.array([["O", oxygen[1], 18, oxygen[2]], ["Al", aluminium[1], 4, aluminium[2]],
                                ["Si", silicon[1], 5, silicon[2]], ["Fe", iron[1], 2, iron[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Ti", "Mn", "Ca", "Na", "K", "H"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Skt"
        #
        # Molar mass
        molar_mass_pure = 2*iron[2] + 4*aluminium[2] + 5*silicon[2] + 18*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[17.7488, 9.827, 9.288], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 109.45*10**9
        # Shear modulus
        G = 60.78*10**9
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_schorl(self):   # Na Fe3 Al6 (B O3)3 Si6 O18 (O H)4
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        boron = PeriodicSystem(name="B").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "B", "O", "Na", "Al", "Si", "Fe"]
        #
        majors_data = np.array([["H", hydrogen[1], 4, hydrogen[2]], ["B", boron[1], 3, boron[2]],
                                ["O", oxygen[1], 3*3+18+4, oxygen[2]], ["Na", sodium[1], 1, sodium[2]],
                                ["Al", aluminium[1], 6, aluminium[2]], ["Si", silicon[1], 6, silicon[2]],
                                ["Fe", iron[1], 3, iron[2]]], dtype=object)
        #
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Mg", "Ca", "Li", "Cr", "Ti", "F", "K"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Tour"
        #
        # Molar mass
        molar_mass_pure = sodium[2] + 3*iron[2] + 6*aluminium[2] + (6*silicon[2] + 18*oxygen[2]) \
                          + 3*(boron[2] + 3*oxygen[2]) + 4*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[15.99, 7.195], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 3
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_elbaite(self, x_value=None, enrichment=None):   # Na Li2.5 Al6.5 (BO3)3 Si6 O18 (OH)4
        # Major elements
        majors_name = ["H", "Li", "B", "O", "Na", "Al", "Si"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Fe", "Mn", "Cu", "Ti", "Ca", "F"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]

        mineral = "Tour"

        ## Molar mass
        condition = False
        while condition == False:
            if x_value == None:
                if enrichment == None:
                    #x = round(rd.uniform(0, 1), 4)
                    x = round(rd.uniform(5/6 - 1/24, 5/6 + 1/24), 4)
                elif enrichment == "Li":
                    x = round(rd.uniform(0.80, 0.85), 4)
                elif enrichment == "Al":
                    x = round(rd.uniform(0.15, 0.20), 4)
            else:
                x = x_value

            majors_data = np.array([
                ["H", self.hydrogen[1], 4, self.hydrogen[2]], ["Li", self.lithium[1], 3*x, self.lithium[2]],
                ["B", self.boron[1], 3, self.boron[2]], ["O", self.oxygen[1], 3*3 + 18 + 4, self.oxygen[2]],
                ["Na", self.sodium[1], 1, self.sodium[2]], ["Al", self.aluminium[1], 6 + 3*(1 - x), self.aluminium[2]],
                ["Si", self.silicon[1], 6, self.silicon[2]]], dtype=object)

            molar_mass_pure = (self.sodium[2] + 3*(x*self.lithium[2] + (1 - x)*self.aluminium[2]) + 6*self.aluminium[2] +
                               (6*self.silicon[2] + 18*self.oxygen[2]) + 3*(self.boron[2] + 3*self.oxygen[2]) +
                               4*(self.oxygen[2] + self.hydrogen[2]))
            molar_mass, amounts = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

            majors_data_ideal = np.array([
                ["H", self.hydrogen[1], 4, self.hydrogen[2]], ["Li", self.lithium[1], 2.5, self.lithium[2]],
                ["B", self.boron[1], 3, self.boron[2]], ["O", self.oxygen[1], 3*3 + 18 + 4, self.oxygen[2]],
                ["Na", self.sodium[1], 1, self.sodium[2]], ["Al", self.aluminium[1], 6 + 0.5, self.aluminium[2]],
                ["Si", self.silicon[1], 6, self.silicon[2]]], dtype=object)

            molar_mass_pure_ideal = (
                    self.sodium[2] + 2.5*self.lithium[2] + 0.5*self.aluminium[2] + 6*self.aluminium[2] +
                    (6*self.silicon[2] + 18*self.oxygen[2]) + 3*(self.boron[2] + 3*self.oxygen[2]) +
                    4*(self.oxygen[2] + self.hydrogen[2]))
            molar_mass_ideal, amounts_ideal = MineralChemistry(
                w_traces=traces_data, molar_mass_pure=molar_mass_pure_ideal,
                majors=majors_data_ideal).calculate_molar_mass()
            element_ideal = [PeriodicSystem(name=amounts_ideal[i][0]).get_data() for i in range(len(amounts_ideal))]
            ## Oxide Composition
            list_oxides = ["H2O", "Li2O", "B2O3", "Na2O", "Al2O3", "SiO2"]
            composition_oxides = {}
            w_oxide_total = 0
            for index, var_oxide in enumerate(list_oxides):
                if index < len(list_oxides) - 1:
                    oxide_data = OxideCompounds(var_compound=var_oxide, var_amounts=amounts).get_composition()
                    if self.impurity == "pure":
                        composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 4)
                    else:
                        composition_oxides[var_oxide] = round(oxide_data["Oxide"][1], 6)
                    w_oxide_total += composition_oxides[var_oxide]
                else:
                    if self.impurity == "pure":
                        composition_oxides[var_oxide] = round(1 - w_oxide_total, 4)
                    else:
                        composition_oxides[var_oxide] = round(1 - w_oxide_total, 6)

            if np.isclose(np.sum(list(composition_oxides.values())), 1.0000) == True:
                condition = True
            else:
                pass
        # Density
        dataV = CrystalPhysics([[15.838, 7.103], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 3
        V_m_pre = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()

        dataRho_ideal = CrystalPhysics([molar_mass_ideal, Z, V])
        rho_ideal = dataRho_ideal.calculate_bulk_density()
        rho_e_ideal = wg(amounts=amounts_ideal, elements=element_ideal, rho_b=rho_ideal).calculate_electron_density()
        V_m = (rho/rho_ideal)*V_m_pre
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        # Results
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["compounds"] = {}
        for index, oxide in enumerate(list_oxides, start=0):
            results["compounds"][oxide] = composition_oxides[oxide]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_liddicoatite(self):   # Ca Li2 Al7 (BO3)3 Si6 O18 (OH)4
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        lithium = PeriodicSystem(name="Li").get_data()
        boron = PeriodicSystem(name="B").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "Li", "B", "O", "Al", "Si", "Ca"]
        #
        majors_data = np.array([["H", hydrogen[1], 4, hydrogen[2]], ["Li", lithium[1], 2, lithium[2]],
                                ["B", boron[1], 3, boron[2]], ["O", oxygen[1], 3*3+18+4, oxygen[2]],
                                ["Al", aluminium[1], 7, aluminium[2]], ["Si", silicon[1], 6, silicon[2]],
                                ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
        #
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Mn", "Fe", "Ti", "Mg", "Na"]
                n = rd.randint(1, len(minors))
                while len(self.traces_list) < n:
                    selection = rd.choice(minors)
                    if selection not in self.traces_list and selection not in majors_name:
                        self.traces_list.append(selection)
                    else:
                        continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        mineral = "Tour"
        #
        # Molar mass
        molar_mass_pure = calcium[2] + 2*lithium[2] + 7*aluminium[2] + (6*silicon[2] + 18*oxygen[2]) \
                          + 3*(boron[2] + 3*oxygen[2]) + 4*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[15.875, 7.126], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 3
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 99.47*10**9    # estimated
        # Shear modulus
        G = 79.85*10**9    # estimated
        # Young's modulus
        E = (9*K*G)/(3*K + G)
        # Poisson's ratio
        nu = (3*K - 2*G)/(2*(3*K + G))
        # vP/vS
        vPvS = ((K + 4/3*G)/G)**0.5
        # P-wave velocity
        vP = ((K + 4/3*G)/rho)**0.5
        # S-wave velocity
        vS = (G/rho)**0.5
        # Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        # Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        results = {}
        results["mineral"] = mineral
        results["M"] = round(molar_mass, 3)
        results["state"] = var_state
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        results["rho"] = round(rho, 4)
        results["rho_e"] = round(rho_e, 4)
        results["V"] = round(V_m, 4)
        results["vP"] = round(vP, 4)
        results["vS"] = round(vS, 4)
        results["vP/vS"] = round(vPvS, 4)
        results["G"] = round(G*10**(-9), 4)
        results["K"] = round(K*10**(-9), 4)
        results["E"] = round(E*10**(-9), 4)
        results["nu"] = round(nu, 4)
        results["GR"] = round(gamma_ray, 4)
        results["PE"] = round(pe, 4)
        results["U"] = round(U, 4)
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results