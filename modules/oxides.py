#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		oxides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		04.06.2021

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry

# OXIDES
class Oxides():
    """ Class that generates geophysical and geochemical data of oxide minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure"):
        self.traces_list = traces_list
        self.impurity = impurity
    #
    def create_quartz(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        majors_name = ["O", "Si"]
        majors_data = np.array([["O", oxygen[1], 2, oxygen[2]], ["Si", silicon[1], 1, silicon[2]]], dtype=object)
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["H", "Al", "Li", "Fe", "Ti", "Na", "Mg", "Ge"]
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
        mineral = "Qz"
        #
        # Molar mass
        molar_mass_pure = silicon[2] + 2*oxygen[2]
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        dataV = CrystalPhysics([[4.9135, 5.4050], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([molar_mass, 3, V])
        rho = dataRho.calculate_bulk_density()
        # Bulk modulus
        K = 29*10**9
        # Shear modulus
        G = 44*10**9
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
        U = pe*rho*10**(-3)
        # Electrical resistivity
        p = 2*10**14
        #
        data.append(mineral)
        data.append(round(molar_mass, 3))
        data.append(round(rho, 2))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 4)])
        data.append([round(vP, 2), round(vS, 2), round(vPvS, 2)])
        data.append([round(gamma_ray, 2), round(pe, 2), round(U, 2), p])
        data.append(amounts)
        #
        return data
    #
    def create_uraninite(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        uranium = PeriodicSystem(name="U").get_data()
        majors_name = ["O", "U"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        #
        data = []
        #
        mineral = "Urn"
        #
        # Molar mass
        majors = [oxygen, uranium]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        for i in range(len(self.traces_list)):
            M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
            weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [2.0, (1-np.sum(x_traces))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[5.4682], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 28*10**(-8)
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data
    #
    def create_magnetite(self):
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Fe"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["Mg", "Zn", "Mn", "Ni", "Cr", "Ti", "V", "Al"]
            n = rd.randint(1, len(minors))
            while len(self.traces_list) < n:
                selection = rd.choice(minors)
                if selection not in self.traces_list and selection not in majors_name:
                    self.traces_list.append(selection)
                else:
                    continue
        #
        data = []
        #
        mineral = "Mag"
        #
        # Molar mass
        majors = [oxygen, iron]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        for i in range(len(self.traces_list)):
            M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
            weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [4.0, (3-np.sum(x_traces))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[8.396], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 176*10**9
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 2850
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data
    #
    def create_hematite(self):   # Fe2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Fe"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        #
        data = []
        #
        mineral = "Hem"
        #
        # Molar mass
        majors = [oxygen, iron]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        x_hydrogen = 0
        for i in range(len(self.traces_list)):
            if self.traces_list[i] == "H":
                x_hydrogen = x_traces[i]
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            else:
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [(3.0-x_hydrogen), (2.0-(np.sum(x_traces)-x_hydrogen))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[5.038, 13.772], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 5*10**6
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data
    #
    def create_corundum(self):   # Al2O3
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        majors_name = ["O", "Al"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        #
        data = []
        #
        mineral = "Crn"
        #
        # Molar mass
        majors = [oxygen, aluminium]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        x_hydrogen = 0
        for i in range(len(self.traces_list)):
            if self.traces_list[i] == "H":
                x_hydrogen = x_traces[i]
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            else:
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [(3.0-x_hydrogen), (2.0-(np.sum(x_traces)-x_hydrogen))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[4.75, 12.982], [], "trigonal"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 6, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = 5*10**6
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data
    #
    def create_wustite(self):   # FeO
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Fe"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        #
        data = []
        #
        mineral = "Wus"
        #
        # Molar mass
        majors = [oxygen, iron]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        x_hydrogen = 0
        for i in range(len(self.traces_list)):
            if self.traces_list[i] == "H":
                x_hydrogen = x_traces[i]
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            else:
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [(1.0-x_hydrogen), (1.0-(np.sum(x_traces)-x_hydrogen))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[4.296], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 4, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
        # Bulk modulus
        K = 147*10**9
        # Shear modulus
        G = 44*10**9
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data
    #
    def create_chromite(self):   # FeCr2O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        chromium = PeriodicSystem(name="Cr").get_data()
        iron = PeriodicSystem(name="Fe").get_data()
        majors_name = ["O", "Cr", "Fe"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        #
        data = []
        #
        mineral = "Chr"
        #
        # Molar mass
        majors = [oxygen, chromium, iron]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        x_hydrogen = 0
        for i in range(len(self.traces_list)):
            if self.traces_list[i] == "H":
                x_hydrogen = x_traces[i]
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            else:
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [(4.0-x_hydrogen), 2.0, (1.0-(np.sum(x_traces)-x_hydrogen))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[8.344], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data
    #
    def create_spinel(self):   # MgAl2O4
        # Major elements
        oxygen = PeriodicSystem(name="O").get_data()
        magnesium = PeriodicSystem(name="Mg").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        majors_name = ["O", "Mg", "Al"]
        # Minor elements
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        #
        data = []
        #
        mineral = "Spl"
        #
        # Molar mass
        majors = [oxygen, magnesium, aluminium]
        weights = []
        traces = [PeriodicSystem(name=i).get_data() for i in self.traces_list]
        x_traces = [round(rd.uniform(0., 0.001), 6) for i in range(len(self.traces_list))]
        M = 0
        x_hydrogen = 0
        for i in range(len(self.traces_list)):
            if self.traces_list[i] == "H":
                x_hydrogen = x_traces[i]
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
            else:
                M += x_traces[i]*PeriodicSystem(name=self.traces_list[i]).get_data()[2]
                weights.append([traces[i][0], int(traces[i][1]), float(x_traces[i])])
        x_majors = [(4.0-x_hydrogen), 2.0, (1.0-(np.sum(x_traces)-x_hydrogen))]
        for i in range(len(majors)):
            M += x_majors[i]*majors[i][2]
            weights.append([majors[i][0], majors[i][1], x_majors[i]])
        weights = np.array(weights, dtype=object)
        weights = weights[weights[:, 1].argsort()]
        element = [PeriodicSystem(name=weights[i][0]).get_data() for i in range(len(weights))]
        amounts = weights[:, 2]
        weights = weights.tolist()
        # Density
        dataV = CrystalPhysics([[8.0898], [], "cubic"])
        V = dataV.calculate_volume()
        dataRho = CrystalPhysics([M, 8, V])
        rho = dataRho.calculate_bulk_density()
        data_rho_e = CrystalPhysics([element, amounts, rho])
        rho_e = data_rho_e.calculate_electron_density()
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
        PE = bg.calculate_pe(self, x_list=amounts, elements_list=element)
        pe = wg(amounts=amounts, elements=element).calculate_pe()
        print("PE", round(pe, 2))
        U = PE*rho_e*10**(-3)
        # Electrical resistivity
        p = None
        #
        data.append(mineral)
        data.append(round(M, 2))
        data.append(round(rho, 1))
        data.append([round(K*10**(-9), 2), round(G*10**(-9), 2), round(E*10**(-9), 2), round(nu, 2), round(vPvS, 2)])
        data.append([round(vP, 1), round(vS, 1)])
        data.append([round(gamma_ray, 2), round(PE, 2), round(U, 2), p])
        data.append(weights)
        #
        return data