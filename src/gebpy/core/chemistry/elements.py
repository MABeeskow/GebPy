#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		elements.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		26.02.2020

# -----------------------------------------------

# MODULES

class elements:
    #
    def __init__(self):
        pass

    #
    def H(self):  # H
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemH = ["H", 1, 1.008, 0.01121, 0.0899, 0.0, 0.0, 0.0, 1270, 0.0]
        #
        return chemH

    #
    def Li(self):  # Li
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        data = ["Li", 3, 6.94, 1.297 * 10 ** (-5), 535, 11, 4.2, 4.9, 5570, 2802, 9.4 * 10 ** (-8)]
        #
        return data

    #
    def B(self):  # B
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        data = ["B", 5, 10.814, 0.439 * 10 ** (-5), 2460, 211, 200, 456, 13934, 9017]
        #
        return data

    #
    def C(self):  # C
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemC = ["C", 6, 12.011, 5.3146 * 10 ** (-6), 2260, 11.9, 8.9, 21.5, 3243, 1984, 1375]
        #
        return chemC
    #
    def N(self):  # N
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        data = ["N", 7, 14.007, 13.54*10**(-6), 1026, 2.44, 0.0, 0.0, 335.25, 0.0, 0.0*10**(-8)]
        #
        return data
    #
    def O(self):  # O
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemO = ["O", 8, 15.999, 1.736 * 10 ** (-5), 1495, 0.0, 0.0, 0.0, 317.5, 0.0, 0.0]
        #
        return chemO

    #
    def F(self):  # F
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMg = ["F", 9, 18.998, 1.120 * 10 ** (-5), 1700, 0, 0, 0, 0, 0.0]
        #
        return chemMg

    #
    def Na(self):  # Na
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemNa = ["Na", 11, 22.990, 2.378 * 10 ** (-5), 968, 6.3, 3.3, 10, 3325, 1847, 4.9 * 10 ** (-8)]
        #
        return chemNa

    #
    def Mg(self):  # Mg
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMg = ["Mg", 12, 24.305, 1.400 * 10 ** (-5), 1738, 45, 17, 45, 6240, 3128, 4.4 * 10 ** (-8)]
        #
        return chemMg

    #
    def Al(self):  # Al
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemAl = ["Al", 13, 26.982, 1.000 * 10 ** (-5), 2700, 76, 26, 70, 6402, 3103, 2.7 * 10 ** (-8)]
        #
        return chemAl

    #
    def Si(self):  # Si
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemSi = ["Si", 14, 28.085, 1.206 * 10 ** (-5), 2330, 100, 62, 160, 8854, 5158, 1]
        #
        return chemSi

    #
    def S(self):  # S
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMg = ["S", 16, 32.059, 1.553 * 10 ** (-5), 1960, 7.7, 3.9, 10.2, 2565, 1411, 1 * 10 ** 15 * 10 ** (-8)]
        #
        return chemMg

    #
    def Cl(self):  # Cl
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemSi = ["Cl", 17, 35.446, 1.739 * 10 ** (-5), 2030, 1.1, 0.0, 0.0, 206, 0.0]
        #
        return chemSi
    #
    def Ar(self):  # Ar
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        data = ["Ar", 18, 39.948, 2.256*10**(-5), 83.81, 8.53*10**(-3), 0.0, 0.0, 319, 0.0, 0.0*10**(-8)]
        #
        return data
    #
    def K(self):  # K
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemK = ["K", 19, 39.098, 4.594 * 10 ** (-5), 856, 3.1, 1.3, 0.0, 2376, 1232, 7.5 * 10 ** (-8)]
        #
        return chemK

    #
    def Ca(self):  # Ca
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCa = ["Ca", 20, 40.078, 2.620 * 10 ** (-5), 1550, 17, 7.4, 20, 4163, 2185, 3.4 * 10 ** (-8)]
        #
        return chemCa
    #
    def Ti(self):  # Ti
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        data = ["Ti", 22, 47.867, 10.64*10**(-6), 4507, 110, 44, 116, 6117, 3125]
        #
        return data
    #
    def Cr(self):  # Cr
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMn = ["Cr", 24, 51.996, 0.723 * 10 ** (-5), 7140, 160, 115, 279, 6625, 4013]
        #
        return chemMn

    #
    def Mn(self):  # Mn
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMn = ["Mn", 25, 54.938, 0.735 * 10 ** (-5), 7440, 120, 90, 198, 5150, 3477]
        #
        return chemMn

    #
    def Fe(self):  # Fe
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        chemFe = ["Fe", 26, 55.845, 7.09 * 10 ** (-6), 7874, 170, 82, 211, 5956, 3227, 10 * 10 ** (-8)]
        #
        return chemFe
    #
    def Co(self):  # Co
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        chem_Co = ["Co", 27, 58.933, 6.67*10**(-6), 8890, 180, 209, 75, 7183, 4849, 6*10**(-8)]
        #
        return chem_Co
    #
    def Cu(self):  # Cu
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCu = ["Cu", 29, 63.546, 0.711 * 10 ** (-5), 8920, 140, 48, 130, 4442, 2320]
        #
        return chemCu

    #
    def As(self):  # As
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemCu = ["As", 33, 74.922, 1.295 * 10 ** (-5), 5727, 22, 3, 8, 2131, 724]
        #
        return chemCu
    #
    def Nb(self):  # Nb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMo = ["Nb", 41, 92.906, 10.83*10**(-6), 8580, 170, 38, 105, 5071, 2104]
        #
        return chemMo

    #
    def Mo(self):  # Mo
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMo = ["Mo", 42, 95.95, 0.938 * 10 ** (-5), 10280, 230, 20, 329, 4882, 1395]
        #
        return chemMo

    #
    def Ag(self):  # Ag
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMo = ["Ag", 47, 107.87, 1.027 * 10 ** (-5), 10490, 100, 30, 83, 3653, 1691]
        #
        return chemMo
    #
    def Sb(self):  # Sb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemMo = ["Sb", 51, 121.76, 1.819 * 10 ** (-5), 6690, 42, 20, 55, 3202, 1728]
        #
        return chemMo
    #
    def Ba(self):  # Ba
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        data = ["Ba", 56, 137.33, 38.16*10**(-6), 3650, 9.6, 4.9, 13, 2102, 1159, 34*10**(-8)]
        #
        return data
    #
    def Ta(self):  # Ta
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        data = ["Ta", 73, 180.95, 10.85*10**(-6), 16680, 200, 69, 186, 4184, 2034]
        #
        return data
    #
    def W(self):  # W
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        data = ["W", 74, 183.84, 9.47 * 10 ** (-6), 19260, 304, 147, 380, 5095, 2763]
        #
        return data
    #
    def Au(self):  # Au
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        data = ["Au", 79, 196.97, 10.21 * 10 ** (-6), 19320, 137, 15, 43, 2851, 881, 2.2 * 10 ** (-8)]
        #
        return data
    #
    def Pb(self):  # Pb
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemPb = ["Pb", 82, 207.2, 1.826 * 10 ** (-5), 11340, 46, 5.75, 16, 2106, 712]
        #
        return chemPb

    #
    def Bi(self):  # Bi
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS]
        chemPb = ["Bi", 83, 208.98, 2.131 * 10 ** (-5), 9780, 31, 12, 32, 2192, 1108]
        #
        return chemPb
    #
    def U(self):  # U
        # [symbol, atomic number, atomic mass, molar volume, density, bulk modulus, shear modulus, young's modulus, vP, vS, resistivity]
        chem_U = ["U", 92, 238.03, 12.49*10**(-6), 18970, 100, 111, 207, 3616, 2419, 28*10**(-8)]
        #
        return chem_U

class PeriodicSystem():
    #
    def __init__(self, name):
        self.name = name
    #
    def get_data(self):
        if self.name in ["H", "Hydrogen", "hydrogen"]:
            data = ["H", 1, 1.008, 0.01121, 0.0899, 0.0, 0.0, 0.0, 1270, 0.0]
        elif self.name in ["Li", "Lithium", "lithium"]:
            data = ["Li", 3, 6.94, 1.297 * 10 ** (-5), 535, 11, 4.2, 4.9, 5570, 2802, 9.4 * 10 ** (-8)]
        elif self.name in ["B", "Boron", "boron"]:
            data = ["B", 5, 10.814, 0.439 * 10 ** (-5), 2460, 211, 200, 456, 13934, 9017]
        elif self.name in ["C", "Carbon", "carbon"]:
            data = ["C", 6, 12.011, 5.3146 * 10 ** (-6), 2260, 11.9, 8.9, 21.5, 3243, 1984, 1375]
        elif self.name in ["N", "Nitrogen", "nitrogen"]:
            data = ["N", 7, 14.007, 13.54*10**(-6), 1026, 2.44, 0.0, 0.0, 335.25, 0.0, 0.0*10**(-8)]
        elif self.name in ["O", "Oxygen", "oxygen"]:
            data = ["O", 8, 15.999, 1.736 * 10 ** (-5), 1495, 0.0, 0.0, 0.0, 317.5, 0.0, 0.0]
        elif self.name in ["F", "Flourine", "flourine"]:
            data = ["F", 9, 18.998, 1.120 * 10 ** (-5), 1700, 0, 0, 0, 0, 0.0]
        elif self.name in ["Na", "Sodium", "sodium"]:
            data = ["Na", 11, 22.990, 2.378 * 10 ** (-5), 968, 6.3, 3.3, 10, 3325, 1847, 4.9*10**(-8)]
        elif self.name in ["Mg", "Magnesium", "magnesium"]:
            data = ["Mg", 12, 24.305, 1.400 * 10 ** (-5), 1738, 45, 17, 45, 6240, 3128, 4.4*10**(-8)]
        elif self.name in ["Al", "Aluminium", "aluminium"]:
            data = ["Al", 13, 26.982, 1.000 * 10 ** (-5), 2700, 76, 26, 70, 6402, 3103, 2.7 * 10 ** (-8)]
        elif self.name in ["Si", "Silicon", "silicon"]:
            data = ["Si", 14, 28.085, 1.206 * 10 ** (-5), 2330, 100, 62, 160, 8854, 5158, 1]
        elif self.name in ["S", "Sulfur", "sulfur"]:
            data = ["S", 16, 32.059, 1.553 * 10 ** (-5), 1960, 7.7, 3.9, 10.2, 2565, 1411, 1*10**(-8)]
        elif self.name in ["Cl", "Chlorine", "chlorine"]:
            data = ["Cl", 17, 35.446, 1.739 * 10 ** (-5), 2030, 1.1, 0.0, 0.0, 206, 0.0]
        elif self.name in ["Ar", "Argon", "argon"]:
            data = ["Ar", 18, 39.948, 2.256*10**(-5), 83.81, 8.53*10**(-3), 0.0, 0.0, 319, 0.0, 0.0*10**(-8)]
        elif self.name in ["K", "Potassium", "potassium"]:
            data = ["K", 19, 39.098, 4.594 * 10 ** (-5), 856, 3.1, 1.3, 0.0, 2376, 1232, 7.5*10**(-8)]
        elif self.name in ["Ca", "Calcium", "calcium"]:
            data = ["Ca", 20, 40.078, 2.620 * 10 ** (-5), 1550, 17, 7.4, 20, 4163, 2185, 3.4*10**(-8)]
        elif self.name in ["Ti", "Titanium", "titanium"]:
            data = ["Ti", 22, 47.867, 10.64*10**(-6), 4507, 110, 44, 116, 6117, 3125]
        elif self.name in ["Cr", "Chromium", "chromium"]:
            data = ["Cr", 24, 51.996, 0.723 * 10 ** (-5), 7140, 160, 115, 279, 6625, 4013]
        elif self.name in ["Mn", "Manganese", "manganese"]:
            data = ["Mn", 25, 54.938, 0.735 * 10 ** (-5), 7440, 120, 90, 198, 5150, 3477]
        elif self.name in ["Fe", "Iron", "iron"]:
            data = ["Fe", 26, 55.845, 7.09 * 10 ** (-6), 7874, 170, 82, 211, 5956, 3227, 10*10**(-8)]
        elif self.name in ["Co", "Cobalt", "cobalt"]:
            data = ["Co", 27, 58.933, 6.67*10**(-6), 8890, 180, 209, 75, 7183, 4849, 6*10**(-8)]
        elif self.name in ["Cu", "Copper", "copper"]:
            data = ["Cu", 29, 63.546, 0.711 * 10 ** (-5), 8920, 140, 48, 130, 4442, 2320]
        elif self.name in ["As", "Arsenic", "arsenic"]:
            data = ["As", 33, 74.922, 1.295 * 10 ** (-5), 5727, 22, 3, 8, 2131, 724]
        elif self.name in ["Nb", "Niobium", "niobium"]:
            data = ["Nb", 41, 92.906, 10.83*10**(-6), 8580, 170, 38, 105, 5071, 2104]
        elif self.name in ["Mo", "Molybdenum", "molybdenum"]:
            data = ["Mo", 42, 95.95, 0.938 * 10 ** (-5), 10280, 230, 20, 329, 4882, 1395]
        elif self.name in ["Ag", "Silver", "silver"]:
            data = ["Ag", 47, 107.87, 1.027 * 10 ** (-5), 10490, 100, 30, 83, 3653, 1691]
        elif self.name in ["Sb", "Antimony", "antimony"]:
            data = ["Sb", 51, 121.76, 1.819 * 10 ** (-5), 6690, 42, 20, 55, 3202, 1728]
        elif self.name in ["Ba", "Barium", "barium"]:
            data = ["Ba", 56, 137.33, 38.16*10**(-6), 3650, 9.6, 4.9, 13, 2102, 1159, 34*10**(-8)]
        elif self.name in ["Ta", "Tantalum", "tantalum"]:
            data = ["Ta", 73, 180.95, 10.85*10**(-6), 16680, 200, 69, 186, 4184, 2034]
        elif self.name in ["W", "Tungsten", "tungsten"]:
            data = ["W", 74, 183.84, 9.47 * 10 ** (-6), 19260, 304, 147, 380, 5095, 2763]
        elif self.name in ["Au", "Gold", "gold"]:
            data = ["Au", 79, 196.97, 10.21 * 10 ** (-6), 19320, 137, 15, 43, 2851, 881, 2.2*10**(-8)]
        elif self.name in ["Pb", "Lead", "lead"]:
            data = ["Pb", 82, 207.2, 1.826 * 10 ** (-5), 11340, 46, 5.75, 16, 2106, 712]
        elif self.name in ["Bi", "Bismuth", "bismuth"]:
            data = ["Bi", 83, 208.98, 2.131 * 10 ** (-5), 9780, 31, 12, 32, 2192, 1108]
        elif self.name in ["U", "Uranium", "uranium"]:
            data = ["U", 92, 238.03, 12.49*10**(-6), 18970, 100, 111, 207, 3616, 2419, 28*10**(-8)]
        #
        return data