#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		silicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		31.10.2021

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry

# TECTOSILICATES
class Tectosilicates:
    """ Class that generates geophysical and geochemical data of tectosiliacte minerals"""
    #
    def __init__(self, traces_list=[], impurity="pure"):
        self.traces_list = traces_list
        self.impurity = impurity
    #
    def create_alkalifeldspar(self, enrichment=None):
        self.enrichment = enrichment
        #
        # Major Elements
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        potassium = PeriodicSystem(name="K").get_data()
        majors_name = ["O", "Na", "Al", "Si", "K"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
        if self.impurity == "random":
            self.traces_list = []
            minors = ["Fe", "Ca", "Na", "Li", "Cs", "Rb", "Pb"]
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
        #
        # Molar mass
        if self.enrichment == None:
            x = round(rd.uniform(0, 1), 2)
            mineral = "Kfs"
        elif self.enrichment == "Na":
            x = round(rd.uniform(0.75, 1), 2)
            mineral = "Kfs"
        elif self.enrichment == "K":
            x = round(rd.uniform(0, 0.25), 2)
            mineral = "Kfs"
        elif self.enrichment == "random":
            x = round(rd.uniform(0, 1), 2)
            if x >= 0.9:
                mineral = "Ab"
            elif x < 0.9 and x >= 0.63:
                mineral = "Ano"
            elif x < 0.63 and x > 0.1:
                mineral = "Sa"
            elif x <= 0.1:
                magicnumber = rd.randint(0, 1)
                if magicnumber == 0:
                    mineral = "Or"
                if magicnumber == 1:
                    mineral = "Mc"
        #
        majors_data = np.array([["O", oxygen[1], 8, oxygen[2]], ["Na", sodium[1], x, sodium[2]],
                                ["Al", aluminium[1], 1, aluminium[2]], ["Si", silicon[1], 3, silicon[2]],
                                ["K", potassium[1], (1-x), potassium[2]]], dtype=object)
        #
        molar_mass_pure = round(x*sodium[2] + (1-x)*potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
        dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
        V_Ab = dataV_Ab.calculate_volume()
        dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
        rho_Ab = dataRho_Ab.calculate_bulk_density()
        rho_e_Ab = wg(amounts=amounts, elements=element, rho_b=rho_Ab).calculate_electron_density()
        M_Or = round(potassium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
        dataV_Or = CrystalPhysics([[8.56, 12.96, 7.21], [116.1], "monoclinic"])
        V_Or = dataV_Or.calculate_volume()
        dataRho_Or = CrystalPhysics([M_Or, 4, V_Or])
        rho_Or = dataRho_Or.calculate_bulk_density()
        rho_e_Or = wg(amounts=amounts, elements=element, rho_b=rho_Or).calculate_electron_density()
        rho = x*rho_Ab + (1-x)*rho_Or
        rho_e = x*rho_e_Ab + (1-x)*rho_e_Or
        # Bulk modulus
        K_Ab = 103.45*10**9
        K_Or = 89.78*10**9
        K = (x*K_Ab + (1-x)*K_Or)
        # Shear modulus
        G_Ab = 69.0*10**9
        G_Or = 61.52*10**9
        G = (x*G_Ab + (1-x)*G_Or)
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
    def create_plagioclase(self, enrichment=None):
        self.enrichment = enrichment
        #
        # Major Elements
        oxygen = PeriodicSystem(name="O").get_data()
        sodium = PeriodicSystem(name="Na").get_data()
        aluminium = PeriodicSystem(name="Al").get_data()
        silicon = PeriodicSystem(name="Si").get_data()
        calcium = PeriodicSystem(name="Ca").get_data()
        majors_name = ["O", "Na", "Al", "Si", "Ca"]
        # Minor elements
        traces_data = []
        if len(self.traces_list) > 0:
            self.impurity = "impure"
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
        data = []
        #
        # Molar mass
        if self.enrichment == None:
            x = round(rd.uniform(0, 1), 2)
            mineral = "Pl"
        elif self.enrichment == "Na":
            x = round(rd.uniform(0.75, 1), 2)
            mineral = "Pl"
        elif self.enrichment == "Ca":
            x = round(rd.uniform(0, 0.25), 2)
            mineral = "Pl"
        elif self.enrichment == "random":
            x = round(rd.uniform(0, 1), 2)
            if x >= 0.9:
                mineral = "Ab"
            elif x < 0.9 and x >= 0.7:
                mineral = "Olg"
            elif x < 0.7 and x >= 0.5:
                mineral = "Andes"
            elif x < 0.5 and x >= 0.3:
                mineral = "Lab"
            elif x < 0.3 and x > 0.1:
                mineral = "Byt"
            elif x <= 0.1:
                mineral = "An"
        #
        majors_data = np.array([["O", oxygen[1], 8, oxygen[2]], ["Na", sodium[1], x, sodium[2]],
                                ["Al", aluminium[1], (2-x), aluminium[2]], ["Si", silicon[1], (2+x), silicon[2]],
                                ["Ca", calcium[1], (1-x), calcium[2]]], dtype=object)
        #
        molar_mass_pure = round(x*sodium[2] + (1-x)*calcium[2] + (2-x)*aluminium[2] + (2+x)*silicon[2] + 8*oxygen[2], 3)
        molar_mass, amounts = MineralChemistry(w_traces=traces_data, molar_mass_pure=molar_mass_pure,
                                      majors=majors_data).calculate_molar_mass()
        element = [PeriodicSystem(name=amounts[i][0]).get_data() for i in range(len(amounts))]
        # Density
        M_Ab = round(sodium[2] + aluminium[2] + 3*silicon[2] + 8*oxygen[2], 3)
        dataV_Ab = CrystalPhysics([[8.15, 12.835, 7.135], [93.85, 116.55, 88.95], "triclinic"])
        V_Ab = dataV_Ab.calculate_volume()
        dataRho_Ab = CrystalPhysics([M_Ab, 4, V_Ab])
        rho_Ab = dataRho_Ab.calculate_bulk_density()
        rho_e_Ab = wg(amounts=amounts, elements=element, rho_b=rho_Ab).calculate_electron_density()
        M_An = round(calcium[2] + 2*aluminium[2] + 2*silicon[2] + 8*oxygen[2], 3)
        dataV_An = CrystalPhysics([[8.18, 12.88, 14.17], [93.2, 115.8, 91.2], "triclinic"])
        V_An = dataV_An.calculate_volume()
        dataRho_An = CrystalPhysics([M_An, 8, V_An])
        rho_An = dataRho_An.calculate_bulk_density()
        rho_e_An = wg(amounts=amounts, elements=element, rho_b=rho_An).calculate_electron_density()
        rho = x*rho_Ab + (1-x)*rho_An
        rho_e = x*rho_e_Ab + (1-x)*rho_e_An
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