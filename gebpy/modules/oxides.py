#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		oxides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.06.2024

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry, TraceElements

# OXIDES
class Oxides():
    """ Class that generates geophysical and geochemical data of oxide minerals"""
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral

        # Chemistry
        boron = PeriodicSystem(name="B").get_data()
        carbon = PeriodicSystem(name="C").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        chlorine = PeriodicSystem(name="Cl").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        antimony = PeriodicSystem(name="Sb").get_data()
        tungsten = PeriodicSystem(name="W").get_data()

        self.boron = ["B", 5, boron[2]]
        self.carbon = ["C", 6, carbon[2]]
        self.oxygen = ["O", 8, oxygen[2]]
        self.sodium = ["Na", 11, sodium[2]]
        self.aluminium = ["Al", 13, aluminium[2]]
        self.silicon = ["Si", 14, silicon[2]]
        self.chlorine = ["Cl", 17, chlorine[2]]
        self.potassium = ["K", 19, potassium[2]]
        self.calcium = ["Ca", 20, calcium[2]]
        self.manganese = ["Mn", 25, manganese[2]]
        self.iron = ["Fe", 26, iron[2]]
        self.antimony = ["Sb", 51, antimony[2]]
        self.tungsten = ["W", 74, tungsten[2]]

    def get_data(self, number=1):
        if self.mineral == "Quartz":
            data = self.create_quartz()
        # Rutile Group
        elif self.mineral == "Argutite":
            data = self.create_argutite()
        elif self.mineral == "Cassiterite":
            data = self.create_cassiterite()
        elif self.mineral == "Paratellurite":
            data = self.create_paratellurite()
        elif self.mineral == "Plattnerite":
            data = self.create_plattnerite()
        elif self.mineral == "Pyrolusite":
            data = self.create_pyrolusite()
        elif self.mineral == "Rutile":
            data = self.create_rutile()
        elif self.mineral == "Stishovite":
            data = self.create_stishovite()
        # Boehmite Group
        elif self.mineral == "Boehmite":
            data = self.create_boehmite()
        # Periclase Group
        elif self.mineral == "Periclase":
            data = self.create_periclase()
        elif self.mineral in ["Wüstite", "Wustite"]:
            data = self.create_wustite()
        elif self.mineral == "Manganosite":
            data = self.create_manganosite()
        elif self.mineral == "Bunsenite":
            data = self.create_bunsenite()
        elif self.mineral == "Monteponite":
            data = self.create_monteponite()
        elif self.mineral == "Lime":
            data = self.create_lime()
        # Scheelite Group
        elif self.mineral == "Wulfenite":
            data = self.create_wulfenite()
        elif self.mineral == "Stolzite":
            data = self.create_stolzite()
        elif self.mineral == "Scheelite":
            data = self.create_scheelite()
        elif self.mineral == "Powellite":
            data = self.create_powellite()
        # Hematite Group
        elif self.mineral == "Hematite":
            data = self.create_hematite()
        elif self.mineral == "Corundum":
            data = self.create_corundum()
        elif self.mineral == "Eskolaite":
            data = self.create_eskolaite()
        elif self.mineral == "Karelianite":
            data = self.create_karelianite()
        elif self.mineral == "Tistarite":
            data = self.create_tistarite()
        # Chromite Group
        elif self.mineral == "Chromite":
            data = self.create_chromite()
        elif self.mineral == "Magnesiochromite":
            data = self.create_magnesiochromite()
        elif self.mineral == "Zincochromite":
            data = self.create_zincochromite()
        elif self.mineral == "Manganochromite":
            data = self.create_manganochromite()
        elif self.mineral == "Niochromite":
            data = self.create_nichromite()
        elif self.mineral == "Cochromite":
            data = self.create_cochromite()
        #
        elif self.mineral == "Crocoite":
            data = self.create_crocoite()
        elif self.mineral == "Uraninite":
            data = self.create_uraninite()
        elif self.mineral == "Magnetite":
            data = self.create_magnetite()
        elif self.mineral == "Spinel":
            data = self.create_spinel()
        elif self.mineral == "Boehmite":
            data = self.create_boehmite()
        elif self.mineral == "Diaspore":
            data = self.create_diaspore()
        elif self.mineral == "Gibbsite":
            data = self.create_gibbsite()
        elif self.mineral == "Cuprite":
            data = self.create_cuprite()
        elif self.mineral == "Goethite":
            data = self.create_goethite()
        elif self.mineral == "Ilmenite":
            data = self.create_ilmenite()
        elif self.mineral == "Brookite":
            data = self.create_brookite()
        elif self.mineral == "Anatase":
            data = self.create_anatase()
        elif self.mineral == "Manganite":
            data = self.create_manganite()
        elif self.mineral == "Groutite":
            data = self.create_groutite()
        elif self.mineral == "Pyrophanite":
            data = self.create_pyrophanite()
        elif self.mineral == "Geikielite":
            data = self.create_geikielite()
        elif self.mineral == "Claudetite":
            data = self.create_claudetite()
        elif self.mineral == "Arsenolite":
            data = self.create_arsenolite()
        elif self.mineral == "Senarmontite":
            data = self.create_senarmontite()
        elif self.mineral == "Valentinite":
            data = self.create_valentinite()
        elif self.mineral == "Bismite":
            data = self.create_bismite()
        elif self.mineral == "Sphaerobismite":
            data = self.create_sphaerobismite()
        elif self.mineral == "Brucite":
            data = self.create_brucite()
        elif self.mineral == "Ulvöspinel":
            data = self.create_ulvoespinel()
        elif self.mineral == "Al-Spinel":
            data = self.create_aluminium_spinel()
        elif self.mineral == "Cr-Spinel":
            data = self.create_chromium_spinel()
        elif self.mineral == "Cuprospinel":
            data = self.create_cuprospinel()
        elif self.mineral == "Jacobsite":
            data = self.create_jacobsite()
        elif self.mineral == "Magnesioferrite":
            data = self.create_magnesioferrite()
        elif self.mineral == "Trevorite":
            data = self.create_trevorite()
        elif self.mineral == "Franklinite":
            data = self.create_franklinite()
        elif self.mineral == "Fe-Spinel":
            data = self.create_iron_spinel()
        elif self.mineral == "Litharge":
            data = self.create_litharge()
        elif self.mineral == "Massicot":
            data = self.create_massicot()
        elif self.mineral == "Minium":
            data = self.create_minium()
        elif self.mineral == "Scrutinyite":
            data = self.create_scrutinyite()
        elif self.mineral == "Zincite":
            data = self.create_zincite()
        elif self.mineral == "Columbite":
            data = self.create_columbite()
        elif self.mineral == "Tantalite":
            data = self.create_tantalite()
        elif self.mineral == "Au2O3":
            data = self.create_gold3oxide()
        elif self.mineral == "Coltan":
            if number > 1:
                data = [self.create_coltan() for n in range(number)]
            else:
                data = self.create_coltan()
        else:
            data = "Nothing found!"
        #
        return data
    #
    def generate_dataset(self, number):
        pure_minerals = ["Quartz", "Hematite", "Corundum"]
        solid_solutions = ["Al-Spinel", "Cr-Spinel", "Fe-Spinel"]
        dataset = {}
        for index in range(number):
            if self.mineral == "Quartz":
                data_mineral = self.create_quartz()
            # Hematite Group
            elif self.mineral == "Hematite":
                data_mineral = self.create_hematite()
            elif self.mineral == "Corundum":
                data_mineral = self.create_corundum()
            elif self.mineral == "Eskolaite":
                data_mineral = self.create_eskolaite()
            elif self.mineral == "Karelianite":
                data_mineral = self.create_karelianite()
            elif self.mineral == "Tistarite":
                data_mineral = self.create_tistarite()
            # Scheelite Group
            elif self.mineral == "Wulfenite":
                data_mineral = self.create_wulfenite()
            elif self.mineral == "Stolzite":
                data_mineral = self.create_stolzite()
            elif self.mineral == "Scheelite":
                data_mineral = self.create_scheelite()
            elif self.mineral == "Powellite":
                data_mineral = self.create_powellite()
            # Chromite Group
            elif self.mineral == "Chromite":
                data_mineral = self.create_chromite()
            elif self.mineral == "Magnesiochromite":
                data_mineral = self.create_magnesiochromite()
            elif self.mineral == "Zincochromite":
                data_mineral = self.create_zincochromite()
            elif self.mineral == "Manganochromite":
                data_mineral = self.create_manganochromite()
            elif self.mineral == "Niochromite":
                data_mineral = self.create_nichromite()
            elif self.mineral == "Cochromite":
                data_mineral = self.create_cochromite()
            #
            elif self.mineral == "Magnetite":
                data_mineral = self.create_magnetite()
            elif self.mineral == "Ilmenite":
                data_mineral = self.create_ilmenite()
            elif self.mineral == "Cassiterite":
                data_mineral = self.create_cassiterite()
            elif self.mineral == "Corundum":
                data_mineral = self.create_corundum()
            elif self.mineral == "Rutile":
                data_mineral = self.create_rutile()
            # Boehmite Group
            elif self.mineral == "Boehmite":
                data_mineral = self.create_boehmite()
            elif self.mineral == "Diaspore":
                data_mineral = self.create_diaspore()
            elif self.mineral == "Pyrolusite":
                data_mineral = self.create_pyrolusite()
            elif self.mineral == "Al-Spinel":
                data_mineral = self.create_aluminium_spinel()
            elif self.mineral == "Cr-Spinel":
                data_mineral = self.create_chromium_spinel()
            elif self.mineral == "Cuprospinel":
                data_mineral = self.create_cuprospinel()
            elif self.mineral == "Spinel":
                data_mineral = self.create_spinel()
            elif self.mineral == "Brucite":
                data_mineral = self.create_brucite()
            elif self.mineral == "Jacobsite":
                data_mineral = self.create_jacobsite()
            elif self.mineral == "Magnesioferrite":
                data_mineral = self.create_magnesioferrite()
            elif self.mineral == "Trevorite":
                data_mineral = self.create_trevorite()
            elif self.mineral == "Franklinite":
                data_mineral = self.create_franklinite()
            elif self.mineral == "Ulvospinel":
                data_mineral = self.create_ulvoespinel()
            elif self.mineral == "Fe-Spinel":
                data_mineral = self.create_iron_spinel()
            elif self.mineral == "Litharge":
                data_mineral = self.create_litharge()
            elif self.mineral == "Massicot":
                data_mineral = self.create_massicot()
            elif self.mineral == "Minium":
                data_mineral = self.create_minium()
            elif self.mineral == "Plattnerite":
                data_mineral = self.create_plattnerite()
            elif self.mineral == "Scrutinyite":
                data_mineral = self.create_scrutinyite()
            elif self.mineral == "Zincite":
                data_mineral = self.create_zincite()
            elif self.mineral == "Columbite":
                data_mineral = self.create_columbite()
            elif self.mineral == "Tantalite":
                data_mineral = self.create_tantalite()
            elif self.mineral == "Coltan":
                data_mineral = self.create_coltan()
            elif self.mineral == "Crocoite":
                data_mineral = self.create_crocoite()
            elif self.mineral == "Wulfenite":
                data_mineral = self.create_wulfenite()
            elif self.mineral == "Goethite":
                data_mineral = self.create_goethite()
            elif self.mineral == "Uraninite":
                data_mineral = self.create_uraninite()
            elif self.mineral == "Wolframite":
                data_mineral = self.create_wolframite()
            elif self.mineral == "Huebnerite":
                data_mineral = self.create_huebnerite()
            elif self.mineral == "Ferberite":
                data_mineral = self.create_ferberite()
            elif self.mineral == "Ferberite-Huebnerite":
                data_mineral = self.create_ferberite_huebnerite()
            elif self.mineral == "Gibbsite":
                data_mineral = self.create_gibbsite()
            elif self.mineral in ["Au2O3", "Au(III)-Oxide"]:
                data_mineral = self.create_gold3oxide()
            elif self.mineral == "Brookite":
                data_mineral = self.create_brookite()
            elif self.mineral == "Anatase":
                data_mineral = self.create_anatase()
            elif self.mineral == "Manganite":
                data_mineral = self.create_manganite()
            elif self.mineral == "Groutite":
                data_mineral = self.create_groutite()
            elif self.mineral == "Pyrophanite":
                data_mineral = self.create_pyrophanite()
            elif self.mineral == "Geikielite":
                data_mineral = self.create_geikielite()
            elif self.mineral == "Claudetite":
                data_mineral = self.create_claudetite()
            elif self.mineral == "Arsenolite":
                data_mineral = self.create_arsenolite()
            elif self.mineral == "Bismite":
                data_mineral = self.create_bismite()
            elif self.mineral == "Sphaerobismite":
                data_mineral = self.create_sphaerobismite()
            # Hydroxides
            elif self.mineral == "Valentinite":
                data_mineral = self.create_valentinite()
            elif self.mineral == "Senarmontite":
                data_mineral = self.create_senarmontite()

            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U", "p"]:
                    if key not in dataset:
                        dataset[key] = [value]
                    else:
                        dataset[key].append(value)
                elif key in ["mineral", "state", "trace elements", "LA-ICP-MS"] and key not in dataset:
                    dataset[key] = value
                elif key in ["chemistry"]:
                    if key not in dataset:
                        dataset[key] = {}
                        for key_2, value_2 in value.items():
                            dataset[key][key_2] = [value_2]
                    else:
                        for key_2, value_2 in value.items():
                            dataset[key][key_2].append(value_2)

            if len(self.traces_list) == 0 and self.mineral in pure_minerals:
                break

        return dataset

    def create_quartz(self, var_T=298.15):
        ## General Information
        name = "Qz"
        ## Trace elements
        element_traces = {
            "4+": ["Ti", "Ge", "Sn", "C"],
            "3+": ["Al", "Fe", "Ga", "As", "B", "P"],
            "2+": ["Mg", "Cu", "Be", "Mn"],
            "1+": ["H", "Li", "Na", "Ag", "K"],
            "All": ["Ti", "Ge", "Sn", "C", "Al", "Fe", "Ga", "As", "B", "P", "H", "Li", "Na", "Ag", "K", "Mg", "Cu",
                    "Be", "Mn"]}

        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity == "pure"
            var_state = "fixed"

        composition_data = TraceElements(tracer=self.traces_list).calculate_composition_quartz()

        ## Molar mass
        molar_mass_pure = self.silicon[2] + 2*self.oxygen[2]
        molar_mass = 0
        amounts = []

        for element in composition_data:
            chem_data = PeriodicSystem(name=element).get_data()
            molar_mass += composition_data[element]["x"]*chem_data[2]
            amounts.append([chem_data[0], chem_data[1], composition_data[element]["w"]])

        magic_factor = molar_mass/molar_mass_pure
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

        ## Density
        dataV = CrystalPhysics([[4.9135, 5.4050], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 3
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)*magic_factor
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()*magic_factor
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()

        ## ELASTICITY
        x = rd.uniform(0, 1)
        ## Bulk modulus
        K_raw = 1.7*x + 36.5
        #K_raw = 38
        K = K_raw*10**9*magic_factor
        ## Shear modulus
        G_raw = -2.3*x + 45.6
        #G_raw = 44
        G = G_raw*10**9*magic_factor
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
        p = 2*10**14*magic_factor
        ## Thermodynamics
        dGf0 = -856281  # J/mol
        dHf0 = -910700  # J/mol
        dS0 = 41.439    # J/(mol K)
        Cp0 = 44.59     # J/(mol K)
        CpT = 39.62 + 44.78*10**(-3)*var_T - 7.45*10**5*var_T**(-2)
        ## LA-ICP-MS
        normalized_sensitivity_Si = 44.54

        ## Data Export
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        results["thermodynamics"] = {"dGf0": dGf0, "dHf0": dHf0, "dS0": dS0, "Cp0": Cp0, "CpT": CpT}
        results["LA-ICP-MS"] = {"Si": normalized_sensitivity_Si}
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_wolframite(self):  # (Fe,Mn) WO4
        name = "Wf"
        # Major elements
        majors_name = ["O", "Mn", "Fe", "W"]

        x = round(rd.uniform(0, 1), 4)

        majors_data = np.array(
            [["O", self.oxygen[1], 4, self.oxygen[2]], ["Mn", self.manganese[1], 1 - x, self.manganese[2]],
             ["Fe", self.iron[1], x, self.iron[2]], ["W", self.tungsten[1], 1, self.tungsten[2]]], dtype=object)

        majors_data_fe = np.array(
            [["O", self.oxygen[1], 4, self.oxygen[2]], ["Mn", self.manganese[1], 0, self.manganese[2]],
             ["Fe", self.iron[1], 1, self.iron[2]], ["W", self.tungsten[1], 1, self.tungsten[2]]], dtype=object)

        majors_data_mn = np.array(
            [["O", self.oxygen[1], 4, self.oxygen[2]], ["Mn", self.manganese[1], 1, self.manganese[2]],
             ["Fe", self.iron[1], 0, self.iron[2]], ["W", self.tungsten[1], 1, self.tungsten[2]]], dtype=object)
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
                minors = ["Nb", "Ta", "Sc", "Sn"]
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
        molar_mass_pure = x*self.iron[2] + (1 - x)*self.manganese[2] + (self.tungsten[2] + 4*self.oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

        molar_mass_pure_fe = self.iron[2] + (self.tungsten[2] + 4*self.oxygen[2])
        molar_mass_fe, amounts_fe = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_fe, majors=majors_data_fe).calculate_molar_mass()
        element_fe = [PeriodicSystem(name=amounts_fe[i][0]).get_data() for i in range(len(amounts_fe))]

        molar_mass_pure_mn = self.manganese[2] + (self.tungsten[2] + 4*self.oxygen[2])
        molar_mass_mn, amounts_mn = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_mn, majors=majors_data_mn).calculate_molar_mass()
        element_mn = [PeriodicSystem(name=amounts_mn[i][0]).get_data() for i in range(len(amounts_mn))]
        # Density
        dataV_fe = CrystalPhysics([[4.76, 5.68, 4.92], [90.016], "monoclinic"])
        V_fe = dataV_fe.calculate_volume()
        Z_fe = 2
        V_m_fe = MineralChemistry().calculate_molar_volume(volume_cell=V_fe, z=Z_fe)
        dataRho_fe = CrystalPhysics([molar_mass_fe, Z_fe, V_fe])
        rho_fe = dataRho_fe.calculate_bulk_density()
        rho_e_fe = wg(amounts=amounts_fe, elements=element_fe, rho_b=rho_fe).calculate_electron_density()

        dataV_mn = CrystalPhysics([[4.8238, 5.7504, 4.9901], [91.18], "monoclinic"])
        V_mn = dataV_mn.calculate_volume()
        Z_mn = 2
        V_m_mn = MineralChemistry().calculate_molar_volume(volume_cell=V_mn, z=Z_mn)
        dataRho_mn = CrystalPhysics([molar_mass_mn, Z_mn, V_mn])
        rho_mn = dataRho_mn.calculate_bulk_density()
        rho_e_mn = wg(amounts=amounts_mn, elements=element_mn, rho_b=rho_mn).calculate_electron_density()

        V_m = x*V_m_fe + (1 - x)*V_m_mn
        rho = x*rho_fe + (1 - x)*rho_mn
        rho_e = x*rho_e_fe + (1 - x)*rho_e_mn
        # Bulk modulus
        K_fe = 141 * 10 ** 9
        K_mn = 127 * 10 ** 9
        K = x * K_fe + (1 - x) * K_mn
        # Shear modulus
        G_fe = 54 * 10 ** 9
        G_mn = 55 * 10 ** 9
        G = x * G_fe + (1 - x) * G_mn
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
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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

        return results

    def create_huebnerite(self):  # Mn WO4
        #
        name = "Hbr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        tungsten = PeriodicSystem(name="W").get_data()
        majors_name = ["O", "Mn", "W"]
        #
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["Mn", manganese[1], 1, manganese[2]],
             ["W", tungsten[1], 1, tungsten[2]]], dtype=object)
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
                minors = ["Nb", "Ta", "Sc", "Sn"]
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
        molar_mass_pure = manganese[2] + (tungsten[2] + 4 * oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        # Density
        dataV = CrystalPhysics([[4.8238, 5.7504, 4.9901], [91.18], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        # Bulk modulus
        K = 127 * 10 ** 9
        # Shear modulus
        G = 55 * 10 ** 9
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
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
    def create_ferberite(self):  # Fe WO4
        #
        name = "Feb"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        tungsten = PeriodicSystem(name="W").get_data()
        majors_name = ["O", "Fe", "W"]
        #
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["Fe", iron[1], 1, iron[2]],
             ["W", tungsten[1], 1, tungsten[2]]], dtype=object)
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
                minors = ["Nb", "Ta", "Sc", "Sn"]
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
        molar_mass_pure = iron[2] + (tungsten[2] + 4 * oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        # Density
        dataV = CrystalPhysics([[4.76, 5.68, 4.92], [90.016], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        # Bulk modulus
        K = 141 * 10 ** 9
        # Shear modulus
        G = 54 * 10 ** 9
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
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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

    def create_ferberite_huebnerite(self):  # (Fe,Mn) WO4
        name = "Feb-Hbr"
        # Major elements
        majors_name = ["O", "Mn", "Fe", "W"]
        x = rd.uniform(0, 1)
        majors_data = np.array([
            ["O", self.oxygen[1], 4,  self.oxygen[2]], ["Mn",  self.manganese[1], (1 - x),  self.manganese[2]],
            ["Fe",  self.iron[1], x,  self.iron[2]], ["W",  self.tungsten[1], 1,  self.tungsten[2]]], dtype=object)

        majors_data_Fe = np.array([
            ["O", self.oxygen[1], 4, self.oxygen[2]], ["Mn", self.manganese[1], 0, self.manganese[2]],
            ["Fe", self.iron[1], 1, self.iron[2]], ["W", self.tungsten[1], 1, self.tungsten[2]]], dtype=object)

        majors_data_Mn = np.array([
            ["O", self.oxygen[1], 4, self.oxygen[2]], ["Mn", self.manganese[1], 1, self.manganese[2]],
            ["Fe", self.iron[1], 0, self.iron[2]], ["W", self.tungsten[1], 1, self.tungsten[2]]], dtype=object)
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
                minors = ["Nb", "Ta", "Sc", "Sn"]
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
        molar_mass_pure = x*self.iron[2] + (1 - x)*self.manganese[2] + (self.tungsten[2] + 4*self.oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]

        molar_mass_pure_Fe = self.iron[2] + (self.tungsten[2] + 4*self.oxygen[2])
        molar_mass_Fe, amounts_Fe = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Fe, majors=majors_data_Fe).calculate_molar_mass()
        element_Fe = [PeriodicSystem(name=amounts_Fe[i][0]).get_data() for i in range(len(amounts_Fe))]

        molar_mass_pure_Mn = self.manganese[2] + (self.tungsten[2] + 4*self.oxygen[2])
        molar_mass_Mn, amounts_Mn = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_Mn, majors=majors_data_Mn).calculate_molar_mass()
        element_Mn = [PeriodicSystem(name=amounts_Mn[i][0]).get_data() for i in range(len(amounts_Mn))]
        # Density
        dataV_Fe = CrystalPhysics([[4.76, 5.68, 4.92], [90.016], "monoclinic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 2
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass_Fe, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts_Fe, elements=element_Fe, rho_b=rho_Fe).calculate_electron_density()

        dataV_Mn = CrystalPhysics([[4.8238, 5.7504, 4.9901], [91.18], "monoclinic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 2
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass_Mn, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts_Mn, elements=element_Mn, rho_b=rho_Mn).calculate_electron_density()

        V_m = x*V_m_Fe + (1 - x)*V_m_Mn
        rho = x*rho_Fe + (1 - x)*rho_Mn
        rho_e = x*rho_e_Fe + (1 - x)*rho_e_Mn
        # Bulk modulus
        K_Fe = 141*10**9
        K_Mn = 127*10**9
        K = x*K_Fe + (1 - x)*K_Mn
        # Shear modulus
        G_Fe = 54*10**9
        G_Mn = 55*10**9
        G = x*G_Fe + (1 - x)*G_Mn
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
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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

        return results

    def create_uraninite(self): # U O2
        #
        name = "Urn"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        uranium = PeriodicSystem(name="U").get_data()
        majors_name = ["O", "Si"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["U", uranium[1], 1, uranium[2]]], dtype=object)
        # Trace elements
        elements_traces = ["Th", "Zr", "Pb", "Ra", "Ac", "Po", "Ce", "Y", "Er", "La"]
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            if self.impurity == "random":
                self.traces_list = []
                minors = ["Th", "Zr", "Pb", "Ra", "Ac", "Po", "Ce", "Y", "Er", "La"]
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
        molar_mass_pure = uranium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.4682], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 200*10**9
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
        p = 28*10**(-8)
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = elements_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_magnetite(self):
        ## General Information
        name = "Mag"
        oxides = ["Fe2O3", "FeO"]
        elements_list = ["Fe", "O"]
        #
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        #
        molar_mass_ideal = 3*iron[2] + 4*oxygen[2]
        molar_mass_fe2o3 = 2*iron[2] + 3*oxygen[2]
        molar_mass_feo = iron[2] + oxygen[2]
        amounts_oxides = {
            "Fe2O3": round(molar_mass_fe2o3/molar_mass_ideal, 6), "FeO": round(molar_mass_feo/molar_mass_ideal, 6)}
        #
        ## Trace elements
        composition_oxides = {}
        for oxide in oxides:
            composition_oxides[oxide] = int(amounts_oxides[oxide]*10**6)
        #
        element_traces = {
            "4+": ["Ti", "V"],
            "3+": ["Cr", "Al"],
            "2+": ["Mg", "Zn", "Mn", "Ni"],
            "All": ["Mg", "Zn", "Mn", "Ni", "Cr", "Ti", "V", "Al"]}
        #
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
            #
            for trace_element, value in self.traces_list.items():
                if trace_element in ["Ti", "V"]:
                    compound = trace_element + str("O2")
                    oxides.append(compound)
                    elements_list.append(trace_element)
                    data_element = PeriodicSystem(name=trace_element).get_data()
                    amount_element = data_element[2]/(data_element[2] + 2*oxygen[2])
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max) / 2
                    sigma = (mean - val_min) / 3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0]/amount_element)
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_oxides[compound] = amount_ppm
                    composition_oxides["Fe2O3"] -= int(amounts_oxides["Fe2O3"]*amount_ppm)
                    composition_oxides["FeO"] -= amount_ppm - int(amounts_oxides["Fe2O3"]*amount_ppm)
                    #
                elif trace_element in ["Cr", "Al"]:
                    compound = trace_element + str("2O3")
                    oxides.append(compound)
                    elements_list.append(trace_element)
                    data_element = PeriodicSystem(name=trace_element).get_data()
                    amount_element = 2*data_element[2]/(2*data_element[2] + 3*oxygen[2])
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max) / 2
                    sigma = (mean - val_min) / 3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0]/amount_element)
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_oxides[compound] = amount_ppm
                    composition_oxides["Fe2O3"] -= amount_ppm
                    #
                elif trace_element in ["Mg", "Zn", "Mn", "Ni"]:
                    compound = trace_element + str("O")
                    oxides.append(compound)
                    elements_list.append(trace_element)
                    data_element = PeriodicSystem(name=trace_element).get_data()
                    amount_element = data_element[2]/(data_element[2] + oxygen[2])
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0]/amount_element)
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_oxides[compound] = amount_ppm
                    composition_oxides["FeO"] -= amount_ppm
            #
        else:
            self.impurity == "pure"
            var_state = "fixed"
        #
        composition_data = TraceElements(
            tracer=self.traces_list).calculate_composition_oxides(
            var_oxides=oxides, var_composition=composition_oxides, var_mineral="Magnetite", var_elements=elements_list)
        #
        ## Molar mass
        molar_mass_pure = molar_mass_ideal
        molar_mass = 0
        amounts = []
        #
        for element in composition_data:
            chem_data = PeriodicSystem(name=element).get_data()
            molar_mass += composition_data[element]["x"] * chem_data[2]
            amounts.append([chem_data[0], chem_data[1], composition_data[element]["w"]])
        #
        magic_factor = molar_mass / molar_mass_pure
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        ## Density
        dataV = CrystalPhysics([[8.396], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z) * magic_factor
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density() * magic_factor
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        ## Bulk modulus
        K = 176*10**9 * magic_factor
        ## Shear modulus
        G = 64*10**9 * magic_factor
        ## Young's modulus
        E = (9 * K * G) / (3 * K + G)
        ## Poisson's ratio
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
        ## vP/vS
        vPvS = ((K + 4 / 3 * G) / G) ** 0.5
        ## P-wave velocity
        vP = ((K + 4 / 3 * G) / rho) ** 0.5
        ## S-wave velocity
        vS = (G / rho) ** 0.5
        ## Gamma ray
        gamma_ray = wg(amounts=amounts, elements=element).calculate_gr()
        ## Photoelectricity
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        U = pe * rho_e * 10 ** (-3)
        ## Electrical resistivity
        p = 2850 * magic_factor
        #
        ## Data Export
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_hematite(self):  # Fe2O3
        #
        results = {}
        results["mineral"] = "Hem"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Fe"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Fe", iron[1], 2, iron[2]]], dtype=object)
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
                minors = ["H", "Ti", "Al", "Mn"]
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
        name = "Hem"
        #
        # Molar mass
        molar_mass_pure = 2*iron[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        #
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.038, 13.772], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 143.58*10**9
        # Shear modulus
        G = 53.43*10**9
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
        p = 5*10**6
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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

        return results

    def create_corundum(self):   # Al2O3
        #
        name = "Crn"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        majors_name = ["O", "Al"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]]], dtype=object)
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
                minors = ["Cr", "Fe", "V", "Ti"]
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
        molar_mass_pure = 2*aluminium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.75, 12.982], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 232*10**9
        # Shear modulus
        G = 147*10**9
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
        p = 5*10**6
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]

        return results
    #
    def create_tistarite(self):   # Ti2O3
        #
        name = "Tta"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        majors_name = ["O", "Ti"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Ti", titanium[1], 2, titanium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = 2*titanium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.158, 13.611], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 216*10**9
        # Shear modulus
        G = 87*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_hematite_group(self):   # X2O3
        #
        name = "Hem-Group"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        vanadium = PeriodicSystem(name="V").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        #
        probabilities = list(np.around(np.random.dirichlet(np.ones(5), size=1)[0], 4))
        for index, value in enumerate(probabilities):
            if index == 0:
                v = value
            elif index == 1:
                w = value
            elif index == 2:
                x = value
            elif index == 3:
                y = value
            elif index == 4:
                z = value
        majors_name = ["O", "Al", "Ti", "V", "Cr", "Fe"]
        majors_data = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 2*v, aluminium[2]],
             ["Ti", titanium[1], 2*w, titanium[2]], ["V", vanadium[1], 2*x, vanadium[2]],
             ["Cr", chromium[1], 2*y, chromium[2]], ["Fe", iron[1], 2*z, iron[2]]], dtype=object)
        majors_data_al = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 2, aluminium[2]],
             ["Ti", titanium[1], 0, titanium[2]], ["V", vanadium[1], 0, vanadium[2]],
             ["Cr", chromium[1], 0, chromium[2]], ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_ti = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 0, aluminium[2]],
             ["Ti", titanium[1], 2, titanium[2]], ["V", vanadium[1], 0, vanadium[2]],
             ["Cr", chromium[1], 0, chromium[2]], ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_v = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 0, aluminium[2]],
             ["Ti", titanium[1], 0, titanium[2]], ["V", vanadium[1], 2, vanadium[2]],
             ["Cr", chromium[1], 0, chromium[2]], ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_cr = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 0, aluminium[2]],
             ["Ti", titanium[1], 0, titanium[2]], ["V", vanadium[1], 0, vanadium[2]],
             ["Cr", chromium[1], 2, chromium[2]], ["Fe", iron[1], 0, iron[2]]], dtype=object)
        majors_data_fe = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Al", aluminium[1], 0, aluminium[2]],
             ["Ti", titanium[1], 0, titanium[2]], ["V", vanadium[1], 0, vanadium[2]],
             ["Cr", chromium[1], 0, chromium[2]], ["Fe", iron[1], 2, iron[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = 2*(v*aluminium[2] + w*titanium[2] + x*vanadium[2] + y*chromium[2] + z*iron[2]) + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        molar_mass_pure_al = 2*aluminium[2] + 3*oxygen[2]
        molar_mass_al, amounts_al = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_al, majors=majors_data_al).calculate_molar_mass()
        element_al = [PeriodicSystem(name=amounts_al[i][0]).get_data() for i in range(len(amounts_al))]
        #
        molar_mass_pure_ti = 2*titanium[2] + 3*oxygen[2]
        molar_mass_ti, amounts_ti = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_ti, majors=majors_data_ti).calculate_molar_mass()
        element_ti = [PeriodicSystem(name=amounts_ti[i][0]).get_data() for i in range(len(amounts_ti))]
        #
        molar_mass_pure_v = 2*vanadium[2] + 3*oxygen[2]
        molar_mass_v, amounts_v = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_v, majors=majors_data_v).calculate_molar_mass()
        element_v = [PeriodicSystem(name=amounts_v[i][0]).get_data() for i in range(len(amounts_v))]
        #
        molar_mass_pure_cr = 2*chromium[2] + 3*oxygen[2]
        molar_mass_cr, amounts_cr = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_cr, majors=majors_data_cr).calculate_molar_mass()
        element_cr = [PeriodicSystem(name=amounts_cr[i][0]).get_data() for i in range(len(amounts_cr))]
        #
        molar_mass_pure_fe = 2*iron[2] + 3*oxygen[2]
        molar_mass_fe, amounts_fe = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure_fe, majors=majors_data_fe).calculate_molar_mass()
        element_fe = [PeriodicSystem(name=amounts_fe[i][0]).get_data() for i in range(len(amounts_fe))]
        #
        # Density
        dataV = CrystalPhysics([[5.158, 13.611], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        dataV_al = CrystalPhysics([[4.751, 12.97], [], "trigonal"])
        V_al = dataV_al.calculate_volume()
        Z_al = 6
        V_m_al = MineralChemistry().calculate_molar_volume(volume_cell=V_al, z=Z_al)
        dataRho_al = CrystalPhysics([molar_mass, Z_al, V_al])
        rho_al = dataRho_al.calculate_bulk_density()
        rho_e_al = wg(amounts=amounts_al, elements=element_al, rho_b=rho_al).calculate_electron_density()
        #
        dataV_ti = CrystalPhysics([[5.158, 13.611], [], "trigonal"])
        V_ti = dataV_ti.calculate_volume()
        Z_ti = 4
        V_m_ti = MineralChemistry().calculate_molar_volume(volume_cell=V_ti, z=Z_ti)
        dataRho_ti = CrystalPhysics([molar_mass, Z_ti, V_ti])
        rho_ti = dataRho_ti.calculate_bulk_density()
        rho_e_ti = wg(amounts=amounts_ti, elements=element_ti, rho_b=rho_ti).calculate_electron_density()
        #
        dataV_v = CrystalPhysics([[4.99, 13.98], [], "trigonal"])
        V_v = dataV_v.calculate_volume()
        Z_v = 6
        V_m_v = MineralChemistry().calculate_molar_volume(volume_cell=V_v, z=Z_v)
        dataRho_v = CrystalPhysics([molar_mass, Z_v, V_v])
        rho_v = dataRho_v.calculate_bulk_density()
        rho_e_v = wg(amounts=amounts_v, elements=element_v, rho_b=rho_v).calculate_electron_density()
        #
        dataV_cr = CrystalPhysics([[4.958, 13.6], [], "trigonal"])
        V_cr = dataV_cr.calculate_volume()
        Z_cr = 6
        V_m_cr = MineralChemistry().calculate_molar_volume(volume_cell=V_cr, z=Z_cr)
        dataRho_cr = CrystalPhysics([molar_mass, Z_cr, V_cr])
        rho_cr = dataRho_cr.calculate_bulk_density()
        rho_e_cr = wg(amounts=amounts_cr, elements=element_cr, rho_b=rho_cr).calculate_electron_density()
        #
        dataV_fe = CrystalPhysics([[5.0317, 13.737], [], "trigonal"])
        V_fe = dataV_fe.calculate_volume()
        Z_fe = 6
        V_m_fe = MineralChemistry().calculate_molar_volume(volume_cell=V_fe, z=Z_fe)
        dataRho_fe = CrystalPhysics([molar_mass, Z_fe, V_fe])
        rho_fe = dataRho_fe.calculate_bulk_density()
        rho_e_fe = wg(amounts=amounts_fe, elements=element_fe, rho_b=rho_fe).calculate_electron_density()
        #
        V_m = v*V_m_al + w*V_m_ti + x*V_m_v + y*V_m_cr + z*V_m_fe
        rho = v*rho_al + w*rho_ti + x*rho_v + y*rho_cr + z*rho_fe
        rho_e = v*rho_e_al + w*rho_e_ti + x*rho_e_v + y*rho_e_cr + z*rho_e_fe
        #
        # Bulk modulus
        K_al = 232*10**9
        K_ti = 216*10**9
        K_v = 202*10**9
        K_cr = 203*10**9
        K_fe = 143.58*10**9
        K = v*K_al + w*K_ti + x*K_v + y*K_cr + z*K_fe
        # Shear modulus
        G_al = 147*10**9
        G_ti = 87*10**9
        G_v = 80*10**9
        G_cr = 113*10**9
        G_fe = 53.43*10**9
        G = v*G_al + w*G_ti + x*G_v + y*G_cr + z*G_fe
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_chromite(self):   # FeCr2O4
        name = "Chr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Cr", "Fe"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Cr", chromium[1], 2, chromium[2]],
                                ["Fe", iron[1], 1, iron[2]]], dtype=object)
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
                minors = ["Mg", "Mn", "Zn", "Al", "Ti"]
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
        molar_mass_pure = 1*iron[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.344], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 147.21*10**9
        # Shear modulus
        G = 55.77*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]

        return results

    def create_manganochromite(self):   # MnCr2O4
        name = "Chr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Cr", "Mn"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Cr", chromium[1], 2, chromium[2]],
                                ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
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
                minors = ["Mg", "Fe", "Ni", "Co", "Zn"]
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
        molar_mass_pure = 1*manganese[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.47], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 147.47*10**9    # estimated
        # Shear modulus
        G = 56.54*10**9     # estimated
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]

        return results

    def create_nichromite(self):   # NiCr2O4
        name = "Chr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["O", "Cr", "Ni"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Cr", chromium[1], 2, chromium[2]],
                                ["Ni", nickel[1], 1, nickel[2]]], dtype=object)
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
                minors = ["Mg", "Fe", "Mn", "Co", "Zn"]
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
        molar_mass_pure = 1*nickel[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.316], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 148.00*10**9    # estimated
        # Shear modulus
        G = 55.72*10**9     # estimated
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]

        return results

    def create_cochromite(self):   # CoCr2O4
        name = "Chr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        cobalt = PeriodicSystem(name="Co").get_data()
        majors_name = ["O", "Cr", "Co"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Cr", chromium[1], 2, chromium[2]],
                                ["Co", cobalt[1], 1, cobalt[2]]], dtype=object)
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
                minors = ["Mg", "Fe", "Mn", "Ni", "Zn"]
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
        molar_mass_pure = 1*cobalt[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.292], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 148.03*10**9    # estimated
        # Shear modulus
        G = 55.67*10**9     # estimated
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]

        return results
    #
    def create_spinel(self):   # MgAl2O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        majors_name = ["O", "Mg", "Al"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]],
                                ["Al", aluminium[1], 2, aluminium[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "4+": ["Ti"],
            "3+": ["Fe",],
            "2+": ["Mn", "Zn", "Ca"],
            "All": ["Ti", "Fe", "Zn", "Mn", "Ca"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["Ti", "Fe", "Zn", "Mn", "Ca"]
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
        mineral = "Spl"
        #
        # Molar mass
        molar_mass_pure = 1*magnesium[2] + 2*aluminium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.0898], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 180*10**9
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
        # Results
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_boehmite(self):  # AlO(OH)
        ## General Information
        name = "Bhm"
        oxides = ["H2O", "Al2O3"]
        elements_list = ["H", "O", "Al"]
        #
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        #
        molar_mass_ideal = aluminium[2] + oxygen[2] + (oxygen[2] + hydrogen[2])
        molar_mass_al2o3 = 2*aluminium[2] + 3*oxygen[2]
        molar_mass_h2o = 2*hydrogen[2] + oxygen[2]
        amounts_oxides = {
            "Al2O3": round(0.5*molar_mass_al2o3/molar_mass_ideal, 6),
            "H2O": round(0.5*molar_mass_h2o/molar_mass_ideal, 6)}
        #
        ## Trace elements
        composition_oxides = {}
        for oxide in oxides:
            composition_oxides[oxide] = int(amounts_oxides[oxide]*10**6)
        #
        element_traces = {
            "4+": ["Si"],
            "3+": ["Fe", "Cr"],
            "2+": ["Mn"],
            "All": ["Fe", "Mn", "Cr", "Si"]}
        #
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
            #
            for trace_element, value in self.traces_list.items():
                if trace_element in ["Si"]:
                    compound = trace_element + str("O2")
                    oxides.append(compound)
                    elements_list.append(trace_element)
                    data_element = PeriodicSystem(name=trace_element).get_data()
                    amount_element = data_element[2]/(data_element[2] + 2*oxygen[2])
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0]/amount_element)
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_oxides[compound] = amount_ppm
                    composition_oxides["Al2O3"] -= amount_ppm
                    #
                elif trace_element in ["Fe", "Cr", "Al"]:
                    compound = trace_element + str("2O3")
                    oxides.append(compound)
                    elements_list.append(trace_element)
                    data_element = PeriodicSystem(name=trace_element).get_data()
                    amount_element = 2*data_element[2]/(2*data_element[2] + 3*oxygen[2])
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0]/amount_element)
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_oxides[compound] = amount_ppm
                    composition_oxides["Al2O3"] -= amount_ppm
                    #
                elif trace_element in ["Mn"]:
                    compound = trace_element + str("O")
                    oxides.append(compound)
                    elements_list.append(trace_element)
                    data_element = PeriodicSystem(name=trace_element).get_data()
                    amount_element = data_element[2]/(data_element[2] + oxygen[2])
                    #
                    val_min = self.traces_list[trace_element]["Min"]
                    val_max = self.traces_list[trace_element]["Max"]
                    mean = (val_min + val_max)/2
                    sigma = (mean - val_min)/3
                    #
                    condition = False
                    while condition == False:
                        amount_ppm = int(np.random.normal(loc=mean, scale=sigma, size=1)[0]/amount_element)
                        if amount_ppm >= 0 and val_min <= amount_ppm <= val_max:
                            condition = True
                    #
                    composition_oxides[compound] = amount_ppm
                    composition_oxides["Al2O3"] -= amount_ppm
            #
        else:
            self.impurity == "pure"
            var_state = "fixed"
        #
        composition_data = TraceElements(
            tracer=self.traces_list).calculate_composition_oxides(
            var_oxides=oxides, var_composition=composition_oxides, var_mineral="Boehmite", var_elements=elements_list)
        #
        ## Molar mass
        molar_mass_pure = molar_mass_ideal
        molar_mass = 0
        amounts = []
        #
        for element in composition_data:
            chem_data = PeriodicSystem(name=element).get_data()
            molar_mass += composition_data[element]["x"] * chem_data[2]
            amounts.append([chem_data[0], chem_data[1], composition_data[element]["w"]])
        #
        magic_factor = molar_mass / molar_mass_pure
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        #
        ## Density
        dataV = CrystalPhysics([[2.868, 12.227, 3.7], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)*magic_factor
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()*magic_factor
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        #
        ## Bulk modulus
        K = 114*10**9*magic_factor
        ## Shear modulus
        G = 82*10**9*magic_factor
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
        p = None
        #
        ## Data Export
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p
        #
        return results
    #
    def create_diaspore(self):  # AlO(OH)
        mineral = "Dsp"
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        majors_name = ["H", "O", "Al"]
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 2, oxygen[2]],
                                ["Al", aluminium[1], 1, aluminium[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "4+": ["Si"],
            "3+": ["Fe", "Cr"],
            "2+": ["Mn"],
            "All": ["Fe", "Mn", "Cr", "Si"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Mn", "Cr", "Si"]
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
        molar_mass_pure = 1*aluminium[2] + 1*oxygen[2] + (oxygen[2]+hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.397, 9.421, 2.8439], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 174.02*10**9
        # Shear modulus
        G = 96.08*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_gibbsite(self): # Al(OH)3
        #
        name = "Gbs"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        majors_name = ["H", "O", "Al"]
        majors_data = np.array([["H", hydrogen[1], 3, hydrogen[2]], ["O", oxygen[1], 3, oxygen[2]],
                                ["Al", aluminium[1], 1, aluminium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        if self.impurity == "random":
            self.traces_list = []
            var_state = "variable"
            minors = ["Fe", "Ga"]
            n = rd.randint(1, len(minors))
            while len(self.traces_list) < n:
                selection = rd.choice(minors)
                if selection not in self.traces_list and selection not in majors_name:
                    self.traces_list.append(selection)
                else:
                    continue
        else:
            var_state = "fixed"
        #
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        for i in range(len(self.traces_list)):
            traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
        if len(traces_data) > 0:
            traces_data = np.array(traces_data, dtype=object)
            traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        # Molar mass
        molar_mass_pure = 1*aluminium[2] + 3*(oxygen[2]+hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.641, 5.07, 9.719], [94.566], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 58.2094*10**9
        # Shear modulus
        G = 29.0147*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results

    def create_cuprite(self):  # Cu2O
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["O", "Cu"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Cu", copper[1], 2, copper[2]]], dtype=object)
        # Minor elements
        element_traces = {"All": []}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

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

        mineral = "Cup"
        #
        # Molar mass
        molar_mass_pure = 2*copper[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.2696], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 111*10**9
        # Shear modulus
        G = 8*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_ilmenite(self):  # FeTiO3
        #
        name = "Ilm"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Ti", "Fe"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Ti", titanium[1], 1, titanium[2]],
                                ["Fe", iron[1], 1, iron[2]]], dtype=object)
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
                minors = ["Mn", "Mg", "V"]
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
        molar_mass_pure = iron[2] + titanium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.093, 14.06], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 166*10**9
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
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_rutile(self):  # TiO2
        #
        results = {}
        results["mineral"] = "Rt"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        majors_name = ["O", "Ti"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Ti", titanium[1], 1, titanium[2]]], dtype=object)
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
                minors = ["Fe", "Ta", "Nb", "Cr", "V", "Sn", "W", "Sb"]
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
        name = "Rt"
        #
        # Molar mass
        molar_mass_pure = titanium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        #
        results["M"] = molar_mass
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.594, 2.958], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 209*10**9
        # Shear modulus
        G = 110*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
    def create_brookite(self):  # TiO2
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        majors_name = ["O", "Ti"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Ti", titanium[1], 1, titanium[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "5+": ["Ta", "Nb"],
            "3+": ["Fe"],
            "All": ["Fe", "Ta", "Nb"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Ta", "Nb"]
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

        mineral = "Brk"
        #
        # Molar mass
        molar_mass_pure = titanium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.456, 9.182, 5.143], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 191*10**9
        # Shear modulus
        G = 78*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_anatase(self):  # TiO2
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        majors_name = ["O", "Ti"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Ti", titanium[1], 1, titanium[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "5+": ["Ta", "Nb"],
            "4+": ["Sn", "V"],
            "3+": ["Fe"],
            "All": ["Fe", "Sn", "V", "Nb"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Sn", "V", "Nb"]
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

        mineral = "Anat"
        #
        # Molar mass
        molar_mass_pure = titanium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.793, 9.51], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 179*10**9
        # Shear modulus
        G = 55*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_manganite(self):  # MnO(OH)
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["H", "O", "Mn"]
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 2, oxygen[2]],
                                ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "3+": ["Fe", "Al"],
            "2+": ["Pb", "Ba", "Cu", "Ca"],
            "All": ["Fe", "Ba", "Pb", "Cu", "Al", "Ca"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Ba", "Pb", "Cu", "Al", "Ca"]
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

        mineral = "Man"
        #
        # Molar mass
        molar_mass_pure = manganese[2] + oxygen[2] + (oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.94, 5.28, 5.74], [90], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 133.29*10**9
        # Shear modulus
        G = 49.80*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_groutite(self):  # MnO(OH)
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["H", "O", "Mn"]
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 2, oxygen[2]],
                                ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "All": []}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

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

        mineral = "Grou"
        #
        # Molar mass
        molar_mass_pure = manganese[2] + oxygen[2] + (oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.56, 10.7, 2.85], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 127.15*10**9
        # Shear modulus
        G = 47.53*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_pyrophanite(self):  # MnTiO3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Ti", "Mn"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Ti", titanium[1], 1, titanium[2]],
                                ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "3+": ["Fe"],
            "2+": ["Zn"],
            "All": ["Fe", "Zn"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Zn"]
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

        mineral = "Pyrph"
        #
        # Molar mass
        molar_mass_pure = manganese[2] + titanium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.137, 14.29], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 159*10**9
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
        # Results
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_geikielite(self):  # MgTiO3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        majors_name = ["O", "Mg", "Ti"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]],
                                ["Ti", titanium[1], 1, titanium[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "3+": ["Fe", "Cr"],
            "2+": ["Mn", "Ca"],
            "All": ["Fe", "Cr", "Mn", "Ca"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Cr", "Mn", "Ca"]
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

        mineral = "Gkl"
        #
        # Molar mass
        molar_mass_pure = magnesium[2] + titanium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.086, 14.093], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 163*10**9
        # Shear modulus
        G = 84*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_eskolaite(self):  # Cr2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        majors_name = ["O", "Cr"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Cr", chromium[1], 2, chromium[2]]], dtype=object)
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
                minors = ["V"]
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
        mineral = "Esk"
        #
        # Molar mass
        molar_mass_pure = 2*chromium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.958, 13.6], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 203*10**9
        # Shear modulus
        G = 113*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_karelianite(self):  # V2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        vanadium = PeriodicSystem(name="V").get_data()
        majors_name = ["O", "V"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["V", vanadium[1], 2, vanadium[2]]], dtype=object)
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

        mineral = "Kar"
        #
        # Molar mass
        molar_mass_pure = 2*vanadium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.99, 13.98], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 6
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 202*10**9
        # Shear modulus
        G = 80*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = mineral
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_claudetite(self):  # As2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["O", "As"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["As", arsenic[1], 2, arsenic[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "All": []}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

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

        mineral = "Clau"
        #
        # Molar mass
        molar_mass_pure = 2*arsenic[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.339, 12.984, 4.5405], [94.26], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 59.37*10**9
        # Shear modulus
        G = 24.61*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_arsenolite(self):  # As2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        arsenic = PeriodicSystem(name="As").get_data()
        majors_name = ["O", "As"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["As", arsenic[1], 2, arsenic[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "All": []}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

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

        mineral = "Arsens"
        #
        # Molar mass
        molar_mass_pure = 2*arsenic[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.0457], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 16
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 51.65*10**9
        # Shear modulus
        G = 21.5*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results

    def create_valentinite(self):  # Sb2O3
        name = "Val"
        # Major elements
        majors_name = ["O", "Sb"]
        majors_data = np.array([
            ["O", self.oxygen[1], 3, self.oxygen[2]], ["Sb", self.antimony[1], 2, self.antimony[2]]], dtype=object)
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
                minors = ["As"]
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
        molar_mass_pure = 2*self.antimony[2] + 3*self.oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.92, 12.46, 5.42], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 19*10**9
        # Shear modulus
        G = 20*10**9
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
        U = pe * rho_e * 10 ** (-3)
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

        return results

    def create_senarmontite(self):  # Sb2O3
        name = "Sen"
        # Major elements
        majors_name = ["O", "Sb"]
        majors_data = np.array([
            ["O", self.oxygen[1], 3, self.oxygen[2]], ["Sb", self.antimony[1], 2, self.antimony[2]]], dtype=object)
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
                minors = ["As"]
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
        molar_mass_pure = 2*self.antimony[2] + 3*self.oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[11.14], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 16
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 64.03*10**9
        # Shear modulus
        G = 23.57*10**9
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
        U = pe * rho_e * 10 ** (-3)
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

        return results

    def create_bismite(self):  # Bi2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        bismuth = PeriodicSystem(name="Bi").get_data()
        majors_name = ["O", "Bi"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Bi", bismuth[1], 2, bismuth[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "All": []}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

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

        mineral = "Bism"
        #
        # Molar mass
        molar_mass_pure = 2*bismuth[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.83, 8.14, 7.48], [67.066], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 54*10**9
        # Shear modulus
        G = 30*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_sphaerobismite(self):  # Bi2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        bismuth = PeriodicSystem(name="Bi").get_data()
        majors_name = ["O", "Bi"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Bi", bismuth[1], 2, bismuth[2]]], dtype=object)
        # Minor elements
        element_traces = {
            "3+": ["As"],
            "All": ["As"]}
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
            var_state = "variable"
        else:
            self.impurity = "pure"
            var_state = "fixed"

        if self.impurity == "random":
            self.traces_list = []
            minors = ["As"]
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

        mineral = "Sphbis"
        #
        # Molar mass
        molar_mass_pure = 2*bismuth[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.08, 6.46], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 21*10**9
        # Shear modulus
        G = 24*10**9
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
        results["state"] = var_state
        results["M"] = molar_mass
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
        results["trace elements"] = element_traces
        if p != None:
            results["p"] = round(p, 4)
        else:
            results["p"] = p

        return results
    #
    def create_pyrolusite(self):  # MnO2
        #
        name = "Prl"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Mn"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = manganese[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.39, 2.86], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 151.32*10**9
        # Shear modulus
        G = 55.66*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_brucite(self):  # Mg(OH)2
        #
        name = "Bru"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        majors_name = ["H", "O", "Mg"]
        majors_data = np.array([["H", hydrogen[1], 2, hydrogen[2]], ["O", oxygen[1], 2, oxygen[2]],
                                ["Mg", magnesium[1], 1, magnesium[2]]], dtype=object)
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
                minors = ["Fe", "Mn", "Zn"]
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
        molar_mass_pure = magnesium[2] + 2*(oxygen[2] + hydrogen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.147, 4.769], [], "trigonal"])
        V = dataV.calculate_volume()
        Z = 1
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 46*10**9
        # Shear modulus
        G = 34*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_aluminium_spinel(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["O", "Mg", "Al", "Mn", "Fe", "Zn"]
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
                minors = ["Ti", "Ca"]
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
        name = "Spi"
        #
        # Molar mass
        w_a = round(rd.uniform(0, 1.0), 4)
        w_b = round(rd.uniform(0, (1-w_a)), 4)
        w_c = round(rd.uniform(0, (1-w_a-w_b)), 4)
        w_d = round(1-w_a-w_b-w_c, 4)
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], w_a, magnesium[2]],
                                ["Al", aluminium[1], 2, aluminium[2]], ["Mn", manganese[1], w_d, manganese[2]],
                                ["Fe", iron[1], w_b, iron[2]], ["Zn", zinc[1], w_c, zinc[2]]], dtype=object)
        #
        molar_mass_pure = w_a*magnesium[2] + w_b*iron[2] + w_c*zinc[2] + w_d*manganese[2] + 2*aluminium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_spi = CrystalPhysics([[8.08], [], "cubic"])
        V_spi = dataV_spi.calculate_volume()
        Z_spi = 8
        V_m_spi = MineralChemistry().calculate_molar_volume(volume_cell=V_spi, z=Z_spi)
        dataRho_spi = CrystalPhysics([molar_mass, Z_spi, V_spi])
        rho_spi = dataRho_spi.calculate_bulk_density()
        rho_e_spi = wg(amounts=amounts, elements=element, rho_b=rho_spi).calculate_electron_density()
        dataV_hc = CrystalPhysics([[8.136], [], "cubic"])
        V_hc = dataV_hc.calculate_volume()
        Z_hc = 8
        V_m_hc = MineralChemistry().calculate_molar_volume(volume_cell=V_hc, z=Z_hc)
        dataRho_hc = CrystalPhysics([molar_mass, Z_hc, V_hc])
        rho_hc = dataRho_hc.calculate_bulk_density()
        rho_e_hc = wg(amounts=amounts, elements=element, rho_b=rho_hc).calculate_electron_density()
        dataV_glx = CrystalPhysics([[8.29], [], "cubic"])
        V_glx = dataV_glx.calculate_volume()
        Z_glx = 8
        V_m_glx = MineralChemistry().calculate_molar_volume(volume_cell=V_glx, z=Z_glx)
        dataRho_glx = CrystalPhysics([molar_mass, Z_glx, V_glx])
        rho_glx = dataRho_glx.calculate_bulk_density()
        rho_e_glx = wg(amounts=amounts, elements=element, rho_b=rho_glx).calculate_electron_density()
        dataV_ghn = CrystalPhysics([[8.062], [], "cubic"])
        V_ghn = dataV_ghn.calculate_volume()
        Z_ghn = 8
        V_m_ghn = MineralChemistry().calculate_molar_volume(volume_cell=V_ghn, z=Z_ghn)
        dataRho_ghn = CrystalPhysics([molar_mass, Z_ghn, V_ghn])
        rho_ghn = dataRho_ghn.calculate_bulk_density()
        rho_e_ghn = wg(amounts=amounts, elements=element, rho_b=rho_ghn).calculate_electron_density()
        #
        V = w_a*V_spi + w_b*V_hc + w_c*V_ghn + w_d*V_glx
        V_m = w_a*V_m_spi + w_b*V_m_hc + w_c*V_m_ghn + w_d*V_m_glx
        rho = w_a*rho_spi + w_b*rho_hc + w_c*rho_ghn + w_d*rho_glx
        rho_e = w_a*rho_e_spi + w_b*rho_e_hc + w_c*rho_e_ghn + w_d*rho_e_glx
        #
        # Bulk modulus
        K_spi = 180*10**9
        K_hc = 177.30*10**9
        K_glx = 171.81*10**9
        K_ghn = 179.24*10**9
        K = w_a*K_spi + w_b*K_hc + w_c*K_ghn + w_d*K_glx
        # Shear modulus
        G_spi = 96*10**9
        G_hc = 93.57*10**9
        G_glx = 92.88*10**9
        G_ghn = 93.52*10**9
        G = w_a*G_spi + w_b*G_hc + w_c*G_ghn + w_d*G_glx
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
    def create_chromium_spinel(self):
        #
        name = "Spi"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["O", "Mg", "Cr", "Fe", "Zn"]
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
        # Molar mass
        w_a = round(rd.uniform(0, 1.0), 4)
        w_b = round(rd.uniform(0, (1-w_a)), 4)
        w_c = round(1-w_a-w_b, 4)
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], w_b, magnesium[2]],
                                ["Cr", chromium[1], 2, chromium[2]], ["Fe", iron[1], w_a, iron[2]],
                                ["Zn", zinc[1], w_c, zinc[2]]], dtype=object)
        #
        molar_mass_pure = w_a*iron[2] + w_b*magnesium[2] + w_c*zinc[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Chr = CrystalPhysics([[8.344], [], "cubic"])
        V_Chr = dataV_Chr.calculate_volume()
        Z_Chr = 8
        V_m_Chr = MineralChemistry().calculate_molar_volume(volume_cell=V_Chr, z=Z_Chr)
        dataRho_Chr = CrystalPhysics([molar_mass, Z_Chr, V_Chr])
        rho_Chr = dataRho_Chr.calculate_bulk_density()
        rho_e_Chr = wg(amounts=amounts, elements=element, rho_b=rho_Chr).calculate_electron_density()
        dataV_MgChr = CrystalPhysics([[8.277], [], "cubic"])
        V_MgChr = dataV_MgChr.calculate_volume()
        Z_MgChr = 8
        V_m_MgChr = MineralChemistry().calculate_molar_volume(volume_cell=V_MgChr, z=Z_MgChr)
        dataRho_MgChr = CrystalPhysics([molar_mass, Z_MgChr, V_MgChr])
        rho_MgChr = dataRho_MgChr.calculate_bulk_density()
        rho_e_MgChr = wg(amounts=amounts, elements=element, rho_b=rho_MgChr).calculate_electron_density()
        dataV_ZnChr = CrystalPhysics([[8.352], [], "cubic"])
        V_ZnChr = dataV_ZnChr.calculate_volume()
        Z_ZnChr = 8
        V_m_ZnChr = MineralChemistry().calculate_molar_volume(volume_cell=V_ZnChr, z=Z_ZnChr)
        dataRho_ZnChr = CrystalPhysics([molar_mass, Z_ZnChr, V_ZnChr])
        rho_ZnChr = dataRho_ZnChr.calculate_bulk_density()
        rho_e_ZnChr = wg(amounts=amounts, elements=element, rho_b=rho_ZnChr).calculate_electron_density()
        #
        V = w_a*V_Chr + w_b*V_MgChr + w_c*V_ZnChr
        V_m = w_a*V_m_Chr + w_b*V_m_MgChr + w_c*V_m_ZnChr
        rho = w_a*rho_Chr + w_b*rho_MgChr + w_c*rho_ZnChr
        rho_e = w_a*rho_e_Chr + w_b*rho_e_MgChr + w_c*rho_e_ZnChr
        #
        # Bulk modulus
        K_Chr = 147.21*10**9
        K_MgChr = 143.27*10**9
        K_ZnChr = 149.02*10**9
        K = w_a*K_Chr + w_b*K_MgChr + w_c*K_ZnChr
        # Shear modulus
        G_Chr = 55.77*10**9
        G_MgChr = 63.29*10**9
        G_ZnChr = 55.11*10**9
        G = w_a*G_Chr + w_b*G_MgChr + w_c*G_ZnChr
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
    def create_cassiterite(self):  # SnO2
        #
        name = "Cst"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        tin = PeriodicSystem(name="Sn").get_data()
        majors_name = ["O", "Sn"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Sn", tin[1], 1, tin[2]]], dtype=object)
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
                minors = ["Fe", "Ta", "Nb", "Zn", "W", "Mn", "Sc", "Ge", "In", "Ga"]
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
        molar_mass_pure = tin[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.738, 3.118], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 172*10**9
        # Shear modulus
        G = 87*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_magnesiochromite(self):  # MgCr2O4
        #
        name = "Mg-Chr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        majors_name = ["O", "Mg", "Cr"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]],
                                ["Cr", chromium[1], 2, chromium[2]]], dtype=object)
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
                minors = ["Fe", "Al"]
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
        molar_mass_pure = magnesium[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.277], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 143.27*10**9
        # Shear modulus
        G = 63.29*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_zincochromite(self):  # ZnCr2O4
        #
        name = "Zn-Chr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["O", "Cr", "Zn"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Cr", chromium[1], 2, chromium[2]],
                                ["Zn", zinc[1], 1, zinc[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = zinc[2] + 2*chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.352], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 149.02*10**9
        # Shear modulus
        G = 55.11*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_cuprospinel(self):  # CuFe2O4
        #
        name = "Cu-Spi"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        majors_name = ["O", "Fe", "Cu"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Fe", iron[1], 2, iron[2]],
                                ["Cu", copper[1], 1, copper[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = copper[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.369], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 159*10**9
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
        # Results
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_jacobsite(self):  # MnFe2O4
        #
        name = "Jcb"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mn", "Fe"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mn", manganese[1], 1, manganese[2]],
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
                minors = ["Zn"]
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
        molar_mass_pure = manganese[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.499], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 158*10**9
        # Shear modulus
        G = 59*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_magnesioferrite(self):  # MgFe2O4
        #
        name = "Mfr"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Mg", "Fe"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]],
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
        # Molar mass
        molar_mass_pure = magnesium[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.366], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 134.17*10**9
        # Shear modulus
        G = 54.36*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_trevorite(self):  # NiFe2O4
        #
        name = "Trv"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["O", "Fe", "Ni"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Fe", iron[1], 2, iron[2]],
                                ["Ni", nickel[1], 1, nickel[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = nickel[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.41], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 135.57*10**9
        # Shear modulus
        G = 47.81*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_franklinite(self):  # ZnFe2O4
        #
        name = "Frk"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["O", "Fe", "Zn"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Fe", iron[1], 2, iron[2]],
                                ["Zn", zinc[1], 1, zinc[2]]], dtype=object)
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
                minors = ["Mn", "Ti", "Al", "Ca"]
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
        molar_mass_pure = zinc[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.42], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 164*10**9
        # Shear modulus
        G = 72*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_ulvoespinel(self):  # TiFe2O4
        #
        name = "Usp"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Ti", "Fe"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Ti", titanium[1], 1, titanium[2]],
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
        # Molar mass
        molar_mass_pure = titanium[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.4596], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 169.07*10**9
        # Shear modulus
        G = 73.17*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_iron_spinel(self):
        #
        name = "Spi"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        titanium = PeriodicSystem(name="Ti").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        copper = PeriodicSystem(name="Cu").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["O", "Mg", "Ti", "Mn", "Fe", "Ni", "Cu", "Zn"]
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
        # Molar mass
        w_a = round(rd.uniform(0, 1.0), 4)
        w_b = round(rd.uniform(0, (1-w_a)), 4)
        w_c = round(rd.uniform(0, (1-w_a-w_b)), 4)
        w_d = round(rd.uniform(0, (1-w_a-w_b-w_c)), 4)
        w_e = round(rd.uniform(0, (1-w_a-w_b-w_c-w_d)), 4)
        w_f = round(rd.uniform(0, (1-w_a-w_b-w_c-w_d-w_e)), 4)
        w_g = round(1-w_a-w_b-w_c-w_d-w_e-w_f, 4)
        #
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Mg", magnesium[1], w_b, magnesium[2]],
                                ["Ti", titanium[1], w_c, titanium[2]], ["Mn", manganese[1], w_e, manganese[2]],
                                ["Fe", iron[1], 2+w_a, iron[2]], ["Ni", nickel[1], w_f, nickel[2]],
                                ["Cu", copper[1], w_g, copper[2]], ["Zn", zinc[1], w_d, zinc[2]]], dtype=object)
        #
        molar_mass_pure = w_a*iron[2] + w_b*magnesium[2] + w_c*titanium[2] + w_d*zinc[2] + w_e*manganese[2] \
                          + w_f*nickel[2] + w_g*copper[2] + 2*iron[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Mag = CrystalPhysics([[8.396], [], "cubic"])
        V_Mag = dataV_Mag.calculate_volume()
        Z_Mag = 8
        V_m_Mag = MineralChemistry().calculate_molar_volume(volume_cell=V_Mag, z=Z_Mag)
        dataRho_Mag = CrystalPhysics([molar_mass, Z_Mag, V_Mag])
        rho_Mag = dataRho_Mag.calculate_bulk_density()
        rho_e_Mag = wg(amounts=amounts, elements=element, rho_b=rho_Mag).calculate_electron_density()
        dataV_Mfr = CrystalPhysics([[8.366], [], "cubic"])
        V_Mfr = dataV_Mfr.calculate_volume()
        Z_Mfr = 8
        V_m_Mfr = MineralChemistry().calculate_molar_volume(volume_cell=V_Mfr, z=Z_Mfr)
        dataRho_Mfr = CrystalPhysics([molar_mass, Z_Mfr, V_Mfr])
        rho_Mfr = dataRho_Mfr.calculate_bulk_density()
        rho_e_Mfr = wg(amounts=amounts, elements=element, rho_b=rho_Mfr).calculate_electron_density()
        dataV_Usp = CrystalPhysics([[8.4596], [], "cubic"])
        V_Usp = dataV_Usp.calculate_volume()
        Z_Usp = 8
        V_m_Usp = MineralChemistry().calculate_molar_volume(volume_cell=V_Usp, z=Z_Usp)
        dataRho_Usp = CrystalPhysics([molar_mass, Z_Usp, V_Usp])
        rho_Usp = dataRho_Usp.calculate_bulk_density()
        rho_e_Usp = wg(amounts=amounts, elements=element, rho_b=rho_Usp).calculate_electron_density()
        dataV_Frk = CrystalPhysics([[8.42], [], "cubic"])
        V_Frk = dataV_Frk.calculate_volume()
        Z_Frk = 8
        V_m_Frk = MineralChemistry().calculate_molar_volume(volume_cell=V_Frk, z=Z_Frk)
        dataRho_Frk = CrystalPhysics([molar_mass, Z_Frk, V_Frk])
        rho_Frk = dataRho_Frk.calculate_bulk_density()
        rho_e_Frk = wg(amounts=amounts, elements=element, rho_b=rho_Frk).calculate_electron_density()
        dataV_Jcb = CrystalPhysics([[8.499], [], "cubic"])
        V_Jcb = dataV_Jcb.calculate_volume()
        Z_Jcb = 8
        V_m_Jcb = MineralChemistry().calculate_molar_volume(volume_cell=V_Jcb, z=Z_Jcb)
        dataRho_Jcb = CrystalPhysics([molar_mass, Z_Jcb, V_Jcb])
        rho_Jcb = dataRho_Jcb.calculate_bulk_density()
        rho_e_Jcb = wg(amounts=amounts, elements=element, rho_b=rho_Jcb).calculate_electron_density()
        dataV_Trv = CrystalPhysics([[8.41], [], "cubic"])
        V_Trv = dataV_Trv.calculate_volume()
        Z_Trv = 8
        V_m_Trv = MineralChemistry().calculate_molar_volume(volume_cell=V_Trv, z=Z_Trv)
        dataRho_Trv = CrystalPhysics([molar_mass, Z_Trv, V_Trv])
        rho_Trv = dataRho_Trv.calculate_bulk_density()
        rho_e_Trv = wg(amounts=amounts, elements=element, rho_b=rho_Trv).calculate_electron_density()
        dataV_CuSpi = CrystalPhysics([[8.369], [], "cubic"])
        V_CuSpi = dataV_CuSpi.calculate_volume()
        Z_CuSpi = 8
        V_m_CuSpi = MineralChemistry().calculate_molar_volume(volume_cell=V_CuSpi, z=Z_CuSpi)
        dataRho_CuSpi = CrystalPhysics([molar_mass, Z_CuSpi, V_CuSpi])
        rho_CuSpi = dataRho_CuSpi.calculate_bulk_density()
        rho_e_CuSpi = wg(amounts=amounts, elements=element, rho_b=rho_CuSpi).calculate_electron_density()
        #
        V = w_a*V_Mag + w_b*V_Mfr + w_c*V_Usp + w_d*V_Frk + w_e*V_Jcb + w_f*V_Trv + w_g*V_CuSpi
        V_m = w_a*V_m_Mag + w_b*V_m_Mfr + w_c*V_m_Usp + w_d*V_m_Frk + w_e*V_m_Jcb + w_f*V_m_Trv + w_g*V_m_CuSpi
        rho = w_a*rho_Mag + w_b*rho_Mfr + w_c*rho_Usp + w_d*rho_Frk + w_e*rho_Jcb + w_f*rho_Trv + w_g*rho_CuSpi
        rho_e = w_a*rho_e_Mag + w_b*rho_e_Mfr + w_c*rho_e_Usp + w_d*rho_e_Frk + w_e*rho_e_Jcb + w_f*rho_e_Trv + w_g*rho_e_CuSpi
        #
        K_Mag = 176*10**9
        K_Mfr = 134.17*10**9
        K_Usp = 169.07*10**9
        K_Frk = 136.02*10**9
        K_Jcb = 158*10**9
        K_Trv = 135.57*10**9
        K_CuSpi = 159*10**9
        K = w_a*K_Mag + w_b*K_Mfr + w_c*K_Usp + w_d*K_Frk + w_e*K_Jcb + w_f*K_Trv + w_g*K_CuSpi
        # Shear modulus
        G_Mag = 64*10**9
        G_Mfr = 54.36*10**9
        G_Usp = 73.17*10**9
        G_Frk = 47.81*10**9
        G_Jcb = 59*10**9
        G_Trv = 47.81*10**9
        G_CuSpi = 47*10**9
        G = w_a*G_Mag + w_b*G_Mfr + w_c*G_Usp + w_d*G_Frk + w_e*G_Jcb + w_f*G_Trv + w_g*G_CuSpi
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
    def create_litharge(self):  # PbO
        #
        name = "Lith"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Pb"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Pb", lead[1], 1, lead[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = lead[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.976, 5.023], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 11*10**9
        # Shear modulus
        G = 5*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_massicot(self):  # PbO
        #
        name = "Mas"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Pb"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Pb", lead[1], 1, lead[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = lead[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.459, 4.723, 5.859], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 54.65*10**9
        # Shear modulus
        G = 20.01*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_minium(self):  # Pb3O4
        #
        name = "Min"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Pb"]
        majors_data = np.array([["O", oxygen[1], 4, oxygen[2]], ["Pb", lead[1], 3, lead[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = 3*lead[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[8.815, 6.565], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 40*10**9
        # Shear modulus
        G = 11*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_plattnerite(self):  # PbO2
        #
        name = "Pltn"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Pb"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Pb", lead[1], 1, lead[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = lead[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.92, 3.367], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 1
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 116*10**9
        # Shear modulus
        G = 45*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_scrutinyite(self):  # PbO2
        #
        name = "Scrt"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Pb"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Pb", lead[1], 1, lead[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = lead[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.971, 5.956, 5.438], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 83.44*10**9
        # Shear modulus
        G = 29.06*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_zincite(self):  # ZnO
        #
        name = "Znc"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        zinc = PeriodicSystem(name="Zn").get_data()
        majors_name = ["O", "Zn"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Zn", zinc[1], 1, zinc[2]]], dtype=object)
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
                minors = ["Mn", "Fe"]
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
        molar_mass_pure = zinc[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.242, 5.176], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 130*10**9
        # Shear modulus
        G = 41*10**9
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_columbite(self):  # (Fe,Mg,Mn) Nb2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        niobium = PeriodicSystem(name="Nb").get_data()
        majors_name = ["O", "Mg", "Mn", "Fe", "Nb"]
        w_Fe = round(rd.uniform(0, 1.0), 4)
        w_Mg = round(rd.uniform(0, (1 - w_Fe)), 4)
        w_Mn = round(1 - w_Fe - w_Mg, 4)
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], w_Mg, magnesium[2]],
                                ["Mn", manganese[1], w_Mn, manganese[2]], ["Fe", iron[1], w_Fe, iron[2]],
                                ["Nb", niobium[1], 2, niobium[2]]], dtype=object)
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
                minors = ["Ta"]
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
        name = "Clb"
        #
        # Molar mass
        molar_mass_pure = (w_Fe*iron[2] + w_Mg*magnesium[2] + w_Mn*manganese[2]) + 2*niobium[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[5.746, 14.308, 5.075], [], "orthorhombic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 4
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Mg = CrystalPhysics([[5.02, 14.17, 5.65], [], "orthorhombic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 4
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[5.767, 14.434, 5.085], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 4
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()
        #
        V_m = w_Fe*V_m_Fe + w_Mg*V_m_Mg + w_Mn*V_m_Mn
        rho = w_Fe*rho_Fe + w_Mg*rho_Mg + w_Mn*rho_Mn
        rho_e = w_Fe*rho_e_Fe + w_Mg*rho_e_Mg + w_Mn*rho_e_Mn
        #
        # Bulk modulus
        K_Fe = 166.16*10**9
        K_Mg = 170.14*10**9
        K_Mn = 162.81*10**9
        K = w_Fe*K_Fe + w_Mg*K_Mg + w_Mn*K_Mn
        # Shear modulus
        G_Fe = 75.01*10**9
        G_Mg = 89.16*10**9
        G_Mn = 73.84*10**9
        G = w_Fe*G_Fe + w_Mg*G_Mg + w_Mn*G_Mn
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
    def create_tantalite(self):  # (Fe,Mg,Mn) Ta2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        tantalum = PeriodicSystem(name="Ta").get_data()
        majors_name = ["O", "Mg", "Mn", "Fe", "Ta"]
        w_Fe = round(rd.uniform(0, 1.0), 4)
        w_Mg = round(rd.uniform(0, (1 - w_Fe)), 4)
        w_Mn = round(1 - w_Fe - w_Mg, 4)
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], w_Mg, magnesium[2]],
                                ["Mn", manganese[1], w_Mn, manganese[2]], ["Fe", iron[1], w_Fe, iron[2]],
                                ["Ta", tantalum[1], 2, tantalum[2]]], dtype=object)
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
                minors = ["Nb"]
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
        name = "Tnt"
        #
        # Molar mass
        molar_mass_pure = (w_Fe*iron[2] + w_Mg*magnesium[2] + w_Mn*manganese[2]) + 2*tantalum[2] + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe = CrystalPhysics([[5.73, 14.24, 5.08], [], "orthorhombic"])
        V_Fe = dataV_Fe.calculate_volume()
        Z_Fe = 4
        V_m_Fe = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe, z=Z_Fe)
        dataRho_Fe = CrystalPhysics([molar_mass, Z_Fe, V_Fe])
        rho_Fe = dataRho_Fe.calculate_bulk_density()
        rho_e_Fe = wg(amounts=amounts, elements=element, rho_b=rho_Fe).calculate_electron_density()
        #
        dataV_Mg = CrystalPhysics([[14.355, 5.735, 5.058], [], "orthorhombic"])
        V_Mg = dataV_Mg.calculate_volume()
        Z_Mg = 4
        V_m_Mg = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg, z=Z_Mg)
        dataRho_Mg = CrystalPhysics([molar_mass, Z_Mg, V_Mg])
        rho_Mg = dataRho_Mg.calculate_bulk_density()
        rho_e_Mg = wg(amounts=amounts, elements=element, rho_b=rho_Mg).calculate_electron_density()
        #
        dataV_Mn = CrystalPhysics([[5.770, 14.465, 5.097], [], "orthorhombic"])
        V_Mn = dataV_Mn.calculate_volume()
        Z_Mn = 4
        V_m_Mn = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn, z=Z_Mn)
        dataRho_Mn = CrystalPhysics([molar_mass, Z_Mn, V_Mn])
        rho_Mn = dataRho_Mn.calculate_bulk_density()
        rho_e_Mn = wg(amounts=amounts, elements=element, rho_b=rho_Mn).calculate_electron_density()
        #
        V_m = w_Fe*V_m_Fe + w_Mg*V_m_Mg + w_Mn*V_m_Mn
        rho = w_Fe*rho_Fe + w_Mg*rho_Mg + w_Mn*rho_Mn
        rho_e = w_Fe*rho_e_Fe + w_Mg*rho_e_Mg + w_Mn*rho_e_Mn
        #
        # Bulk modulus
        K_Fe = 199.31*10**9
        K_Mg = 192.05*10**9
        K_Mn = 201.87*10**9
        K = w_Fe*K_Fe + w_Mg*K_Mg + w_Mn*K_Mn
        # Shear modulus
        G_Fe = 88.85*10**9
        G_Mg = 97.27*10**9
        G_Mn = 94.55*10**9
        G = w_Fe*G_Fe + w_Mg*G_Mg + w_Mn*G_Mn
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
    def create_coltan(self):  # (Fe,Mg,Mn) (Nb,Ta)2 O6
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        niobium = PeriodicSystem(name="Nb").get_data()
        tantalum = PeriodicSystem(name="Ta").get_data()
        majors_name = ["O", "Mg", "Mn", "Fe", "Nb", "Ta"]
        w_Fe = round(rd.uniform(0, 1.0), 4)
        w_Mg = round(rd.uniform(0, (1 - w_Fe)), 4)
        w_Mn = round(1 - w_Fe - w_Mg, 4)
        w_Nb = round(rd.uniform(0, 1.0), 4)
        w_Ta = round(1 - w_Nb, 4)
        majors_data = np.array([["O", oxygen[1], 6, oxygen[2]], ["Mg", magnesium[1], w_Mg, magnesium[2]],
                                ["Mn", manganese[1], w_Mn, manganese[2]], ["Fe", iron[1], w_Fe, iron[2]],
                                ["Nb", niobium[1], 2*w_Nb, niobium[2]], ["Ta", tantalum[1], 2*w_Ta, tantalum[2]]], dtype=object)
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
        name = "Clt"
        #
        # Molar mass
        molar_mass_pure = (w_Fe*iron[2] + w_Mg*magnesium[2] + w_Mn*manganese[2]) \
                          + 2*(w_Nb*niobium[2] + w_Ta*tantalum[2]) + 6*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV_Fe_Clb = CrystalPhysics([[5.746, 14.308, 5.075], [], "orthorhombic"])
        V_Fe_Clb = dataV_Fe_Clb.calculate_volume()
        Z_Fe_Clb = 4
        V_m_Fe_Clb = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe_Clb, z=Z_Fe_Clb)
        dataRho_Fe_Clb = CrystalPhysics([molar_mass, Z_Fe_Clb, V_Fe_Clb])
        rho_Fe_Clb = dataRho_Fe_Clb.calculate_bulk_density()
        rho_e_Fe_Clb = wg(amounts=amounts, elements=element, rho_b=rho_Fe_Clb).calculate_electron_density()
        #
        dataV_Mg_Clb = CrystalPhysics([[5.02, 14.17, 5.65], [], "orthorhombic"])
        V_Mg_Clb = dataV_Mg_Clb.calculate_volume()
        Z_Mg_Clb = 4
        V_m_Mg_Clb = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg_Clb, z=Z_Mg_Clb)
        dataRho_Mg_Clb = CrystalPhysics([molar_mass, Z_Mg_Clb, V_Mg_Clb])
        rho_Mg_Clb = dataRho_Mg_Clb.calculate_bulk_density()
        rho_e_Mg_Clb = wg(amounts=amounts, elements=element, rho_b=rho_Mg_Clb).calculate_electron_density()
        #
        dataV_Mn_Clb = CrystalPhysics([[5.767, 14.434, 5.085], [], "orthorhombic"])
        V_Mn_Clb = dataV_Mn_Clb.calculate_volume()
        Z_Mn_Clb = 4
        V_m_Mn_Clb = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn_Clb, z=Z_Mn_Clb)
        dataRho_Mn_Clb = CrystalPhysics([molar_mass, Z_Mn_Clb, V_Mn_Clb])
        rho_Mn_Clb = dataRho_Mn_Clb.calculate_bulk_density()
        rho_e_Mn_Clb = wg(amounts=amounts, elements=element, rho_b=rho_Mn_Clb).calculate_electron_density()
        #
        dataV_Fe_Tnt = CrystalPhysics([[5.746, 14.308, 5.075], [], "orthorhombic"])
        V_Fe_Tnt = dataV_Fe_Tnt.calculate_volume()
        Z_Fe_Tnt = 4
        V_m_Fe_Tnt = MineralChemistry().calculate_molar_volume(volume_cell=V_Fe_Tnt, z=Z_Fe_Tnt)
        dataRho_Fe_Tnt = CrystalPhysics([molar_mass, Z_Fe_Tnt, V_Fe_Tnt])
        rho_Fe_Tnt = dataRho_Fe_Tnt.calculate_bulk_density()
        rho_e_Fe_Tnt = wg(amounts=amounts, elements=element, rho_b=rho_Fe_Tnt).calculate_electron_density()
        #
        dataV_Mg_Tnt = CrystalPhysics([[5.02, 14.17, 5.65], [], "orthorhombic"])
        V_Mg_Tnt = dataV_Mg_Tnt.calculate_volume()
        Z_Mg_Tnt = 4
        V_m_Mg_Tnt = MineralChemistry().calculate_molar_volume(volume_cell=V_Mg_Tnt, z=Z_Mg_Tnt)
        dataRho_Mg_Tnt = CrystalPhysics([molar_mass, Z_Mg_Tnt, V_Mg_Tnt])
        rho_Mg_Tnt = dataRho_Mg_Tnt.calculate_bulk_density()
        rho_e_Mg_Tnt = wg(amounts=amounts, elements=element, rho_b=rho_Mg_Tnt).calculate_electron_density()
        #
        dataV_Mn_Tnt = CrystalPhysics([[5.767, 14.434, 5.085], [], "orthorhombic"])
        V_Mn_Tnt = dataV_Mn_Tnt.calculate_volume()
        Z_Mn_Tnt = 4
        V_m_Mn_Tnt = MineralChemistry().calculate_molar_volume(volume_cell=V_Mn_Tnt, z=Z_Mn_Tnt)
        dataRho_Mn_Tnt = CrystalPhysics([molar_mass, Z_Mn_Tnt, V_Mn_Tnt])
        rho_Mn_Tnt = dataRho_Mn_Tnt.calculate_bulk_density()
        rho_e_Mn_Tnt = wg(amounts=amounts, elements=element, rho_b=rho_Mn_Tnt).calculate_electron_density()
        #
        V_m = w_Fe*(w_Nb*V_m_Fe_Clb + w_Ta*V_m_Fe_Tnt) + w_Mg*(w_Nb*V_m_Mg_Clb + w_Ta*V_m_Mg_Tnt) \
              + w_Mn*(w_Nb*V_m_Mn_Clb + w_Ta*V_m_Mn_Tnt)
        rho = w_Fe*(w_Nb*rho_Fe_Clb + w_Ta*rho_Fe_Tnt) + w_Mg*(w_Nb*rho_Mg_Clb + w_Ta*rho_Mg_Tnt) \
              + w_Mn*(w_Nb*rho_Mn_Clb + w_Ta*rho_Mn_Tnt)
        rho_e = w_Fe*(w_Nb*rho_e_Fe_Clb + w_Ta*rho_e_Fe_Tnt) + w_Mg*(w_Nb*rho_e_Mg_Clb + w_Ta*rho_e_Mg_Tnt) \
                + w_Mn*(w_Nb*rho_e_Mn_Clb + w_Ta*rho_e_Mn_Tnt)
        #
        # Bulk modulus
        K_Fe_Clb = 166.16*10**9
        K_Mg_Clb = 170.14*10**9
        K_Mn_Clb = 162.81*10**9
        K_Fe_Tnt = 199.31*10**9
        K_Mg_Tnt = 192.05*10**9
        K_Mn_Tnt = 201.87*10**9
        K = w_Fe*(w_Nb*K_Fe_Clb + w_Ta*K_Fe_Tnt) + w_Mg*(w_Nb*K_Mg_Clb + w_Ta*K_Mg_Tnt) \
            + w_Mn*(w_Nb*K_Mn_Clb + w_Ta*K_Mn_Tnt)
        # Shear modulus
        G_Fe_Clb = 75.01*10**9
        G_Mg_Clb = 89.16*10**9
        G_Mn_Clb = 73.84*10**9
        G_Fe_Tnt = 88.85*10**9
        G_Mg_Tnt = 97.27*10**9
        G_Mn_Tnt = 94.55*10**9
        G = w_Fe*(w_Nb*G_Fe_Clb + w_Ta*G_Fe_Tnt) + w_Mg*(w_Nb*G_Mg_Clb + w_Ta*G_Mg_Tnt) \
            + w_Mn*(w_Nb*G_Mn_Clb + w_Ta*G_Mn_Tnt)
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
    def create_argutite(self):   # GeO2
        #
        name = "Argt"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        germanium = PeriodicSystem(name="Ge").get_data()
        majors_name = ["O", "Ge"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Ge", germanium[1], 1, germanium[2]]], dtype=object)
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
                minors = ["Zn", "Mn", "Fe"]
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
        molar_mass_pure = germanium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.3963, 2.8626], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 206*10**9
        # Shear modulus
        G = 123*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_paratellurite(self):   # TeO2
        #
        name = "Prtl"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        tellurium = PeriodicSystem(name="Te").get_data()
        majors_name = ["O", "Te"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Te", tellurium[1], 1, tellurium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = tellurium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.81, 7.613], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 30*10**9
        # Shear modulus
        G = 28*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_stishovite(self):   # SiO2
        #
        name = "Stv"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Si"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Si", silicon[1], 1, silicon[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = silicon[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.179, 2.6649], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 272*10**9
        # Shear modulus
        G = 203*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_baddeleyite(self):   # ZrO2
        #
        name = "Bdy"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        zirconium = PeriodicSystem(name="Zr").get_data()
        majors_name = ["O", "Zr"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Zr", zirconium[1], 1, zirconium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = zirconium[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.1477, 5.203, 5.3156], [99.38], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 183*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_bunsenite(self):   # NiO
        #
        name = "Bsn"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        nickel = PeriodicSystem(name="Ni").get_data()
        majors_name = ["O", "Ni"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Ni", nickel[1], 1, nickel[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = nickel[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.1768], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 131*10**9
        # Shear modulus
        G = 98*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_periclase(self):   # MgO
        #
        name = "Per"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        majors_name = ["O", "Mg"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Mg", magnesium[1], 1, magnesium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = magnesium[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.203], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 151*10**9
        # Shear modulus
        G = 119*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_manganosite(self):   # MnO
        #
        name = "Mns"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        manganese = PeriodicSystem(name="Mn").get_data()
        majors_name = ["O", "Mn"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Mn", manganese[1], 1, manganese[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = manganese[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.436], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 155.1*10**9
        # Shear modulus
        G = 79.9*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_monteponite(self):   # CdO
        #
        name = "Mntp"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        cadmium = PeriodicSystem(name="Cd").get_data()
        majors_name = ["O", "Mn"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Cd", cadmium[1], 1, cadmium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = cadmium[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.689], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 126*10**9
        # Shear modulus
        G = 45*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_lime(self):   # CaO
        #
        name = "Lm"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Ca"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Ca", calcium[1], 1, calcium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = calcium[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.797], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 105*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_wustite(self):   # FeO
        #
        name = "Wus"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Fe"]
        majors_data = np.array([["O", oxygen[1], 1, oxygen[2]], ["Fe", iron[1], 1, iron[2]]], dtype=object)
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
                minors = ["Mg", "Mn", "Ni"]
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
        molar_mass_pure = iron[2] + oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.296], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 180*10**9
        # Shear modulus
        G = 60*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_avicennite(self):   # Tl2O3
        #
        name = "Avcn"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        thallium = PeriodicSystem(name="Tl").get_data()
        majors_name = ["O", "Tl"]
        majors_data = np.array([["O", oxygen[1], 3, oxygen[2]], ["Tl", thallium[1], 2, thallium[2]]], dtype=object)
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
        # Molar mass
        molar_mass_pure = 2*thallium[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[10.5468], [], "cubic"])
        V = dataV.calculate_volume()
        Z = 10
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 109*10**9
        # Shear modulus
        G = 32*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_crocoite(self):   # PbCrO4
        #
        name = "Crc"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Cr", "Pb"]
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["Cr", chromium[1], 1, chromium[2]], ["Pb", lead[1], 1, lead[2]]],
            dtype=object)
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
                minors = ["Zn", "S"]
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
        molar_mass_pure = lead[2] + chromium[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[7.12, 7.44, 6.8], [77.55], "monoclinic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 33*10**9
        # Shear modulus
        G = 18*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_stolzite(self):   # PbWO4
        #
        name = "Sz"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        tungsten = PeriodicSystem(name="W").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "W", "Pb"]
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["W", tungsten[1], 1, tungsten[2]], ["Pb", lead[1], 1, lead[2]]],
            dtype=object)
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
                minors = ["Mo"]
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
        molar_mass_pure = lead[2] + tungsten[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.46, 12.05], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 81*10**9
        # Shear modulus
        G = 30*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_scheelite(self):  # CaWO4
        #
        name = "Sch"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        tungsten = PeriodicSystem(name="W").get_data()
        majors_name = ["O", "Ca", "W"]
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["Ca", calcium[1], 1, calcium[2]], ["W", tungsten[1], 1, tungsten[2]]],
            dtype=object)
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
                minors = ["Mo", "Nb", "Ta"]
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
        molar_mass_pure = calcium[2] + tungsten[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.242, 11.372], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 81*10**9
        # Shear modulus
        G = 33*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_wulfenite(self):   # PbMoO4
        #
        name = "Wul"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        molybdenum = PeriodicSystem(name="Mo").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        majors_name = ["O", "Mo", "Pb"]
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["Mo", molybdenum[1], 1, molybdenum[2]], ["Pb", lead[1], 1, lead[2]]],
            dtype=object)
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
                minors = ["W", "Ca", "V", "As", "Cr", "W", "Ti"]
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
        molar_mass_pure = lead[2] + molybdenum[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.435, 12.11], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 48.9958*10**9
        # Shear modulus
        G = 21.99845*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_powellite(self):  # CaMoO4
        #
        name = "Pwl"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        molybdenum = PeriodicSystem(name="Mo").get_data()
        majors_name = ["O", "Ca", "Mo"]
        majors_data = np.array(
            [["O", oxygen[1], 4, oxygen[2]], ["Ca", calcium[1], 1, calcium[2]],
             ["Mo", molybdenum[1], 1, molybdenum[2]]], dtype=object)
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
                minors = ["W"]
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
        molar_mass_pure = calcium[2] + molybdenum[2] + 4*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[5.23, 11.44], [], "tetragonal"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 75*10**9
        # Shear modulus
        G = 33*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_goethite(self):   # FeO(OH)
        #
        name = "Goe"
        #
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["H", "O", "Fe"]
        majors_data = np.array(
            [["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 2, oxygen[2]], ["Fe", iron[1], 1, iron[2]]],
            dtype=object)
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
                minors = ["Mn"]
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
        molar_mass_pure = iron[2] + oxygen[2] + (hydrogen[2] + oxygen[2])
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.596, 9.957, 3.021], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 4
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 93.2*10**9
        # Shear modulus
        G = 62.2*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
    def create_gold3oxide(self):   # Au2O3
        #
        name = "Au2O3"
        #
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        gold = PeriodicSystem(name="Au").get_data()
        majors_name = ["O", "Au"]
        majors_data = np.array(
            [["O", oxygen[1], 3, oxygen[2]], ["Au", gold[1], 2, gold[2]]],
            dtype=object)
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
                minors = []
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
        molar_mass_pure = 2*gold[2] + 3*oxygen[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[3.91, 10.54, 12.88], [], "orthorhombic"])
        V = dataV.calculate_volume()
        Z = 8
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 81*10**9
        # Shear modulus
        G = 38*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = var_state
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
class RutileGroup:
    #
    def __init__(self):
        self.data_argutite = Oxides(mineral="Argutite").get_data()
        self.data_cassiterite = Oxides(mineral="Cassiterite").get_data()
        self.data_paratellurite = Oxides(mineral="Paratellurite").get_data()
        self.data_plattnerite = Oxides(mineral="Plattnerite").get_data()
        self.data_pyrolusite = Oxides(mineral="Pyrolusite").get_data()
        self.data_rutile = Oxides(mineral="Rutile").get_data()
        self.data_stishovite = Oxides(mineral="Stishovite").get_data()
        self.data_minerals = {
            "Ge": self.data_argutite, "Sn": self.data_cassiterite, "Te": self.data_paratellurite,
            "Pb": self.data_plattnerite, "Mn": self.data_pyrolusite, "Ti": self.data_rutile, "Si": self.data_stishovite}
    #
    def create_rutile_group(self):   # XO2
        #
        name = "Rt-Group"
        #
        # Major elements
        list_X = ["Ge", "Sn", "Te", "Pb", "Mn", "Ti", "Si"]
        fractions_X = {}
        probabilities = list(np.around(np.random.dirichlet(np.ones(7), size=1)[0], 4))
        for index, value in enumerate(probabilities):
            fractions_X[list_X[index]] = value
        #
        oxygen = PeriodicSystem(name="O").get_data()
        majors_name = ["O", "Zr"]
        major_data = [["O", oxygen[1], 2, oxygen[2]]]
        for element in list_X:
            element_data = PeriodicSystem(name=element).get_data()
            major_data.append([element, element_data[1], fractions_X[element], element_data[2]])
        majors_data = np.array(major_data, dtype=object)
        #
        # Minor elements
        traces_data = []
        #
        # Molar mass
        molar_mass_pure = 2*oxygen[2]
        for element, value in fractions_X.items():
            molar_mass_pure += value*PeriodicSystem(name=element).get_data()[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = 0
        rho_e = 0
        V_m = 0
        for key, value in fractions_X.items():
            rho += value*self.data_minerals[key]["rho"]
            rho_e += value*self.data_minerals[key]["rho_e"]
            V_m += value*self.data_minerals[key]["V"]
        # Bulk modulus and Shear modulus
        K = 0
        G = 0
        for key, value in fractions_X.items():
            K += value*self.data_minerals[key]["K"]
            G += value*self.data_minerals[key]["G"]
        K = K*10**9
        G = G*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = "fixed"
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
class PericlaseGroup:
    #
    def __init__(self):
        self.data_periclase = Oxides(mineral="Periclase").get_data()
        self.data_wustite = Oxides(mineral="Wustite").get_data()
        self.data_manganosite = Oxides(mineral="Manganosite").get_data()
        self.data_bunsenite = Oxides(mineral="Bunsenite").get_data()
        self.data_monteponite = Oxides(mineral="Monteponite").get_data()
        self.data_lime = Oxides(mineral="Lime").get_data()
        self.data_minerals = {
            "Mg": self.data_periclase, "Fe": self.data_wustite, "Mn": self.data_manganosite,
            "Ni": self.data_bunsenite, "Cd": self.data_monteponite, "Ca": self.data_lime}
    #
    def create_periclase_group(self):   # XO
        #
        name = "Per-Group"
        #
        # Major elements
        list_X = ["Mg", "Fe", "Mn", "Ni", "Cd", "Ca"]
        fractions_X = {}
        probabilities = list(np.around(np.random.dirichlet(np.ones(len(list_X)), size=1)[0], 4))
        for index, value in enumerate(probabilities):
            fractions_X[list_X[index]] = value
        #
        oxygen = PeriodicSystem(name="O").get_data()
        major_data = [["O", oxygen[1], 1, oxygen[2]]]
        for element in list_X:
            element_data = PeriodicSystem(name=element).get_data()
            major_data.append([element, element_data[1], fractions_X[element], element_data[2]])
        majors_data = np.array(major_data, dtype=object)
        #
        # Minor elements
        traces_data = []
        #
        # Molar mass
        molar_mass_pure = oxygen[2]
        for element, value in fractions_X.items():
            molar_mass_pure += value*PeriodicSystem(name=element).get_data()[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = 0
        rho_e = 0
        V_m = 0
        for key, value in fractions_X.items():
            rho += value*self.data_minerals[key]["rho"]
            rho_e += value*self.data_minerals[key]["rho_e"]
            V_m += value*self.data_minerals[key]["V"]
        # Bulk modulus and Shear modulus
        K = 0
        G = 0
        for key, value in fractions_X.items():
            K += value*self.data_minerals[key]["K"]
            G += value*self.data_minerals[key]["G"]
        K = K*10**9
        G = G*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = "fixed"
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results
    #
class WulfeniteGroup:
    #
    def __init__(self):
        self.data_wulfenite = Oxides(mineral="Wulfenite").get_data()
        self.data_stolzite = Oxides(mineral="Stolzite").get_data()
        self.data_crocoite = Oxides(mineral="Crocoite").get_data()
        self.data_minerals = {
            "Cr": self.data_crocoite, "Mo": self.data_wulfenite, "W": self.data_stolzite}
    #
    def create_wulfenite_group(self):   # XYO4
        #
        name = "Wul-Group"
        #
        # Major elements
        list_X = ["Cr", "Mo", "W"]
        fractions_X = {}
        probabilities = list(np.around(np.random.dirichlet(np.ones(len(list_X)), size=1)[0], 4))
        for index, value in enumerate(probabilities):
            fractions_X[list_X[index]] = value
        #
        oxygen = PeriodicSystem(name="O").get_data()
        lead = PeriodicSystem(name="Pb").get_data()
        major_data = [["O", oxygen[1], 4, oxygen[2]], ["Pb", lead[1], 1, lead[2]]]
        for element in list_X:
            element_data = PeriodicSystem(name=element).get_data()
            major_data.append([element, element_data[1], fractions_X[element], element_data[2]])
        majors_data = np.array(major_data, dtype=object)
        #
        # Minor elements
        traces_data = []
        #
        # Molar mass
        molar_mass_pure = 4*oxygen[2] + lead[2]
        for element, value in fractions_X.items():
            molar_mass_pure += value*PeriodicSystem(name=element).get_data()[2]
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        rho = 0
        rho_e = 0
        V_m = 0
        for key, value in fractions_X.items():
            rho += value*self.data_minerals[key]["rho"]
            rho_e += value*self.data_minerals[key]["rho_e"]
            V_m += value*self.data_minerals[key]["V"]
        # Bulk modulus and Shear modulus
        K = 0
        G = 0
        for key, value in fractions_X.items():
            K += value*self.data_minerals[key]["K"]
            G += value*self.data_minerals[key]["G"]
        K = K*10**9
        G = G*10**9
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
        ## RESULTS
        results = {}
        results["mineral"] = name
        results["state"] = "fixed"
        results["M"] = molar_mass
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
        try:
            results["p"] = round(p, 4)
        except:
            results["p"] = p
        element_list = np.array(amounts)[:, 0]
        results["chemistry"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = amounts[index][2]
        #
        return results