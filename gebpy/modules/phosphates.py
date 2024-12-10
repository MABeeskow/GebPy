#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		phosphates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		16.11.2023

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from scipy import stats
from modules.chemistry import PeriodicSystem, DataProcessing, OxideCompounds
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry, TraceElements
#
# Phosphates
class Phosphates:
    """ Class that generates geophysical and geochemical data of phosphate minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure", data_type=False, mineral=None):
        self.traces_list = traces_list
        self.impurity = impurity
        self.data_type = data_type
        self.mineral = mineral
    #
    def generate_dataset(self, number):
        dataset = {}
        #
        for index in range(number):
            if self.mineral == "Apatite":
                data_mineral = self.create_aptite()
            elif self.mineral == "Fluoroapatite":
                data_mineral = self.create_aptite_f()
            elif self.mineral == "Chloroapatite":
                data_mineral = self.create_aptite_cl()
            elif self.mineral == "Hydroxyapatite":
                data_mineral = self.create_aptite_oh()
            #
            for key, value in data_mineral.items():
                if key in ["M", "rho", "rho_e", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U", "p"]:
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
    def create_aptite_f(self):   # Ca5 F (PO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        fluorine = PeriodicSystem(name="F").get_data()
        phosphorus = PeriodicSystem(name="P").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "F", "P", "Ca"]
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["F", fluorine[1], 1, fluorine[2]],
                                ["P", phosphorus[1], 3, phosphorus[2]], ["Ca", calcium[1], 5, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "pure":
            var_state = "fixed"
        else:
            var_state = "variable"
            self.traces_list = []
            minors_x = ["Ti", "Zr", "Hf", "Th"]         # mainly 4+
            minors_y = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Cr", "As"]  # mainly 3+
            minors_z = ["Cl", "H", "Rb"]                # mainly 1+
            minors_w = ["Mn", "Co", "Sr", "Ba", "Pb"]   # mainly 2+
            if self.impurity == "random":
                n_x = rd.randint(0, len(minors_x))
                n_y = rd.randint(0, len(minors_y))
                n_w = rd.randint(0, len(minors_w))
                if n_x > 0:
                    selection_x = rd.sample(minors_x, n_x)
                    self.traces_list.extend(selection_x)
                if n_y > 0 and n_w == 0:
                    n_z = rd.randint(1, n_y)
                    selection_y = rd.sample(minors_y, n_y)
                    selection_z = rd.sample(minors_z, n_z)
                    self.traces_list.extend(selection_y)
                    self.traces_list.extend(selection_z)
                if n_w > 0 and n_y == 0:
                    n_z = rd.randint(1, n_w)
                    selection_w = rd.sample(minors_w, n_w)
                    selection_z = rd.sample(minors_z, n_z)
                    self.traces_list.extend(selection_w)
                    self.traces_list.extend(selection_z)
                if n_y > 0 and n_w > 0:
                    if n_y + n_w <= len(minors_z):
                        n_z = rd.randint(1, (n_y + n_w))
                    else:
                        n_z = len(minors_z)
                    selection_y = rd.sample(minors_y, n_y)
                    selection_w = rd.sample(minors_w, n_w)
                    selection_z = rd.sample(minors_z, n_z)
                    self.traces_list.extend(selection_y)
                    self.traces_list.extend(selection_w)
                    self.traces_list.extend(selection_z)
            elif self.impurity != "random":
                self.traces_list = []
                for element in self.impurity:
                    if element in minors_x:
                        self.traces_list.append(element)
                    elif element in minors_y:
                        self.traces_list.append(element)
                    elif element in minors_z:
                        self.traces_list.append(element)
                    elif element in minors_w:
                        self.traces_list.append(element)
                # minors = ["Cl", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Mn", "H"]
                # n = rd.randint(1, len(minors))
                # while len(self.traces_list) < n:
                #     selection = rd.choice(minors)
                #     if selection not in self.traces_list and selection not in majors_name:
                #         self.traces_list.append(selection)
                #     else:
                #         continue
            traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
            x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
            for i in range(len(self.traces_list)):
                traces_data.append([str(self.traces_list[i]), int(traces[i][1]), float(x_traces[i])])
            if len(traces_data) > 0:
                traces_data = np.array(traces_data, dtype=object)
                traces_data = traces_data[traces_data[:, 1].argsort()]
        #
        data = []
        mineral = "Ap"
        #
        # Molar mass
        try:
            molar_mass_pure = 5*calcium[2] + fluorine[2] + 3*(phosphorus[2] + 4*oxygen[2])
            molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                                   majors=majors_data).calculate_molar_mass()
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        except:
            compositon_data = TraceElements(tracer=self.traces_list).calculate_composition_apatite_f()
            molar_mass = 0
            amounts = []
            for element in compositon_data:
                chem_data = PeriodicSystem(name=element).get_data()
                molar_mass += compositon_data[element]["x"]*chem_data[2]
                amounts.append([chem_data[0], chem_data[1], compositon_data[element]["w"]])
            element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.367, 6.884], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 83*10**9
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
        #
        if self.data_type == False:
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
            results["M"] = molar_mass
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
    def create_aptite_cl(self):   # Ca5 Cl (PO4)3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        phosphorus = PeriodicSystem(name="P").get_data()
        chlorine = PeriodicSystem(name="Cl").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "P", "Cl", "Ca"]
        majors_data = np.array([["O", oxygen[1], 12, oxygen[2]], ["P", phosphorus[1], 3, phosphorus[2]],
                                ["Cl", chlorine[1], 1, chlorine[2]], ["Ca", calcium[1], 5, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["F", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Mn", "H"]
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
        data = []
        mineral = "Ap"
        #
        # Molar mass
        molar_mass_pure = 5*calcium[2] + chlorine[2] + 3*(phosphorus[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.598, 6.776], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 98.70*10**9
        # Shear modulus
        G = 61.17*10**9
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
    def create_aptite_oh(self):   # Ca5 OH (PO4)3
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        phosphorus = PeriodicSystem(name="P").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "O", "P", "Ca"]
        majors_data = np.array([["H", hydrogen[1], 1, hydrogen[2]], ["O", oxygen[1], 13, oxygen[2]],
                                ["P", phosphorus[1], 3, phosphorus[2]], ["Ca", calcium[1], 5, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["F", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Mn", "Cl"]
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
        data = []
        mineral = "Ap"
        #
        # Molar mass
        molar_mass_pure = 5*calcium[2] + (hydrogen[2] + oxygen[2]) + 3*(phosphorus[2] + 4*oxygen[2])
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                               majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[9.418, 6.875], [], "hexagonal"])
        V = dataV.calculate_volume()
        Z = 2
        V_m = MineralChemistry().calculate_molar_volume(volume_cell=V, z=Z)
        dataRho = CrystalPhysics([molar_mass, Z, V*10**(6)])
        rho = dataRho.calculate_bulk_density()
        rho_e = wg(amounts=amounts, elements=element, rho_b=rho).calculate_electron_density()
        # Bulk modulus
        K = 105.01*10**9
        # Shear modulus
        G = 62.69*10**9
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
    def create_aptite(self):   # Ca5 (F,Cl,OH) (PO4)3
        mineral = "Ap"
        # Major elements
        hydrogen = PeriodicSystem(name="H").get_data()
        oxygen = PeriodicSystem(name="O").get_data()
        fluorine = PeriodicSystem(name="F").get_data()
        phosphorus = PeriodicSystem(name="P").get_data()
        chlorine = PeriodicSystem(name="Cl").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["H", "O", "F", "P", "Cl", "Ca"]
        a = round(rd.uniform(0, 1), 4)
        b = round(rd.uniform(0, (1-a)), 4)
        c = round(1-a-b, 4)
        majors_data = np.array([["H", hydrogen[1], c, hydrogen[2]], ["O", oxygen[1], 12+c, oxygen[2]],
                                ["F", fluorine[1], a, fluorine[2]], ["P", phosphorus[1], 3, phosphorus[2]],
                                ["Cl", chlorine[1], b, chlorine[2]], ["Ca", calcium[1], 5, calcium[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Y", "Er", "Mn"]
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
        molar_mass_pure = (5*calcium[2] + a*fluorine[2] + b*chlorine[2] + c*(hydrogen[2] + oxygen[2]) + 3*
                           (phosphorus[2] + 4*oxygen[2]))
        molar_mass, amounts = MineralChemistry(
            w_traces=traces_data, molar_mass_pure=molar_mass_pure, majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Oxide Composition
        list_oxides = ["H2O", "P2O5", "CaO", "Cl", "F"]
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
        dataV_F = CrystalPhysics([[9.367, 6.884], [], "hexagonal"])
        V_F = dataV_F.calculate_volume()
        Z_F = 2
        V_m_F = MineralChemistry().calculate_molar_volume(volume_cell=V_F, z=Z_F)*10**(6)
        dataRho_F = CrystalPhysics([molar_mass, Z_F, V_F*10**(6)])
        rho_F = dataRho_F.calculate_bulk_density()
        rho_e_F = wg(amounts=amounts, elements=element, rho_b=rho_F).calculate_electron_density()
        #
        dataV_Cl = CrystalPhysics([[9.598, 6.776], [], "hexagonal"])
        V_Cl = dataV_Cl.calculate_volume()
        Z_Cl = 2
        V_m_Cl = MineralChemistry().calculate_molar_volume(volume_cell=V_Cl, z=Z_Cl)*10**(6)
        dataRho_Cl = CrystalPhysics([molar_mass, Z_Cl, V_Cl*10**(6)])
        rho_Cl = dataRho_Cl.calculate_bulk_density()
        rho_e_Cl = wg(amounts=amounts, elements=element, rho_b=rho_Cl).calculate_electron_density()
        #
        dataV_OH = CrystalPhysics([[9.418, 6.875], [], "hexagonal"])
        V_OH = dataV_OH.calculate_volume()
        Z_OH = 2
        V_m_OH = MineralChemistry().calculate_molar_volume(volume_cell=V_OH, z=Z_OH)*10**(6)
        dataRho_OH = CrystalPhysics([molar_mass, Z_OH, V_OH*10**(6)])
        rho_OH = dataRho_OH.calculate_bulk_density()
        rho_e_OH = wg(amounts=amounts, elements=element, rho_b=rho_OH).calculate_electron_density()
        #
        V_m = a*V_m_F + b*V_m_Cl + c*V_m_OH
        rho = a*rho_F + b*rho_Cl + c*rho_OH
        rho_e = a*rho_e_F + b*rho_e_Cl + c*rho_e_OH
        # Bulk modulus
        K_F = 83*10**9
        K_Cl = 98.70*10**9
        K_OH = 105.01*10**9
        K = a*K_F + b*K_Cl + c*K_OH
        # Shear modulus
        G_F = 41*10**9
        G_Cl = 61.17*10**9
        G_OH = 62.69*10**9
        G = a*G_F + b*G_Cl + c*G_OH
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