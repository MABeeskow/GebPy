#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		chemistry.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.05.2021

# -----------------------------------------------

# MODULES
import numpy as np

# CLASSES
class PeriodicSystem():
    """ Class that generates the chemical and physical data of the elements"""
    #
    def __init__(self, name=None, atomicnumber=None):
        """
        :param name: name of the chemical element
        :param atomicnumber: atomic number of the chemical element
        """
        self.name = name
        self.atomicnumber = atomicnumber
    #
    def get_data(self):
        data = []
        # [symbol, atomic number, atomic mass, density, bulk modulus, shear modulus, young's modulus, vP,
        # vS, resistivity]
        if self.name in ["H", "Hydrogen", "hydrogen"] or self.atomicnumber == 1:
            mass_molar = 1.008
            density = 0.0899
            bulk_mod = round(1.4*10**(-4), 5)
            shear_mod = 0.0
            young_mod = (9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod)
            v_p = round(np.sqrt((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(density)), 2)
            v_s = np.sqrt((shear_mod*10**9)/(density))
            resistivity = 0.0
            thermal_cond = 0.1805
            data = ["H", self.atomicnumber, mass_molar, density, bulk_mod, shear_mod, young_mod, v_p, v_s, resistivity,
                    thermal_cond]
        elif self.name in ["He", "Helium", "helium"] or self.atomicnumber == 2:
            mass_molar = 4.0026
            density = 0.1785
            bulk_mod = round(1.01*10**(-5), 7)
            shear_mod = 0.0
            young_mod = (9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod)
            v_p = round(np.sqrt((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(density)), 2)
            v_s = np.sqrt((shear_mod*10**9)/(density))
            resistivity = 0.0
            thermal_cond = 0.1513
            data = ["He", self.atomicnumber, mass_molar, density, bulk_mod, shear_mod, young_mod, v_p, v_s, resistivity,
                    thermal_cond]
        elif self.name in ["Li", "Lithium", "lithium"] or self.atomicnumber == 3:
            mass_molar = 6.938
            density = 530
            bulk_mod = 11
            shear_mod = 4.2
            young_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 1)
            v_p = round(np.sqrt((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(density)), 2)
            v_s = round(np.sqrt((shear_mod*10**9)/(density)), 2)
            resistivity = round(9.4*10**(-8), 9)
            thermal_cond = 85
            data = ["Li", self.atomicnumber, mass_molar, density, bulk_mod, shear_mod, young_mod, v_p, v_s, resistivity,
                    thermal_cond]
        elif self.name in ["Be", "Beryllium", "beryllium"] or self.atomicnumber == 4:
            mass_molar = 9.0122
            density = 1850
            bulk_mod = 130
            shear_mod = 132
            young_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 1)
            v_p = round(np.sqrt((bulk_mod*10**9 + 4/3*shear_mod*10**9)/(density)), 2)
            v_s = round(np.sqrt((shear_mod*10**9)/(density)), 2)
            resistivity = round(4.0*10**(-8), 9)
            thermal_cond = 190
            data = ["Be", self.atomicnumber, mass_molar, density, bulk_mod, shear_mod, young_mod, v_p, v_s, resistivity,
                    thermal_cond]
        elif self.name in ["B", "Boron", "boron"] or self.atomicnumber == 5:
            mass_molar = 10.806
            density = 2460
            bulk_mod = 320
            v_p = 16200
            shear_mod = round((3*(density*v_p**2 - bulk_mod*10**9))/(4)*10**(-9), 1)
            young_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 1)
            v_s = round(np.sqrt((shear_mod*10**9)/(density)), 2)
            resistivity = 10000
            thermal_cond = 27
            data = ["B", self.atomicnumber, mass_molar, density, bulk_mod, shear_mod, young_mod, v_p, v_s, resistivity,
                    thermal_cond]
        elif self.name in ["C", "Carbon", "carbon"] or self.atomicnumber == 6:
            mass_molar = 12.009
            density = 3510
            bulk_mod = 33
            v_p = 18350
            shear_mod = round((3*(density*v_p**2 - bulk_mod*10**9))/(4)*10**(-9), 1)
            young_mod = round((9*bulk_mod*shear_mod)/(3*bulk_mod + shear_mod), 1)
            v_s = round(np.sqrt((shear_mod*10**9)/(density)), 2)
            resistivity = 1*10**(-5)
            thermal_cond = 140
            data = ["C", self.atomicnumber, mass_molar, density, bulk_mod, shear_mod, young_mod, v_p, v_s, resistivity,
                    thermal_cond]
        elif self.name in ["N", "Nitrogen", "nitrogen"] or self.atomicnumber == 7:
            data = ["N", 7, 14.007, 13.54*10**(-6), 1026, 2.44, 0.0, 0.0, 335.25, 0.0, 0.0*10**(-8)]
        elif self.name in ["O", "Oxygen", "oxygen"] or self.atomicnumber == 8:
            data = ["O", 8, 15.999, 1.736 * 10 ** (-5), 1495, 0.0, 0.0, 0.0, 317.5, 0.0, 0.0]
        elif self.name in ["F", "Flourine", "flourine"] or self.atomicnumber == 9:
            data = ["F", 9, 18.998, 1.120 * 10 ** (-5), 1700, 0, 0, 0, 0, 0.0]
        elif self.name in ["Na", "Sodium", "sodium"] or self.atomicnumber == 11:
            data = ["Na", 11, 22.990, 2.378 * 10 ** (-5), 968, 6.3, 3.3, 10, 3325, 1847, 4.9*10**(-8)]
        elif self.name in ["Mg", "Magnesium", "magnesium"] or self.atomicnumber == 12:
            data = ["Mg", 12, 24.305, 1.400 * 10 ** (-5), 1738, 45, 17, 45, 6240, 3128, 4.4*10**(-8)]
        elif self.name in ["Al", "Aluminium", "aluminium"] or self.atomicnumber == 13:
            data = ["Al", 13, 26.982, 1.000 * 10 ** (-5), 2700, 76, 26, 70, 6402, 3103, 2.7 * 10 ** (-8)]
        elif self.name in ["Si", "Silicon", "silicon"] or self.atomicnumber == 14:
            data = ["Si", 14, 28.085, 1.206 * 10 ** (-5), 2330, 100, 62, 160, 8854, 5158, 1]
        elif self.name in ["S", "Sulfur", "sulfur"] or self.atomicnumber == 16:
            data = ["S", 16, 32.059, 1.553 * 10 ** (-5), 1960, 7.7, 3.9, 10.2, 2565, 1411, 1*10**(-8)]
        elif self.name in ["Cl", "Chlorine", "chlorine"] or self.atomicnumber == 17:
            data = ["Cl", 17, 35.446, 1.739 * 10 ** (-5), 2030, 1.1, 0.0, 0.0, 206, 0.0]
        elif self.name in ["Ar", "Argon", "argon"] or self.atomicnumber == 18:
            data = ["Ar", 18, 39.948, 2.256*10**(-5), 83.81, 8.53*10**(-3), 0.0, 0.0, 319, 0.0, 0.0*10**(-8)]
        elif self.name in ["K", "Potassium", "potassium"] or self.atomicnumber == 19:
            data = ["K", 19, 39.098, 4.594 * 10 ** (-5), 856, 3.1, 1.3, 0.0, 2376, 1232, 7.5*10**(-8)]
        elif self.name in ["Ca", "Calcium", "calcium"] or self.atomicnumber == 20:
            data = ["Ca", 20, 40.078, 2.620 * 10 ** (-5), 1550, 17, 7.4, 20, 4163, 2185, 3.4*10**(-8)]
        elif self.name in ["Ti", "Titanium", "titanium"] or self.atomicnumber == 22:
            data = ["Ti", 22, 47.867, 10.64*10**(-6), 4507, 110, 44, 116, 6117, 3125]
        elif self.name in ["Cr", "Chromium", "chromium"] or self.atomicnumber == 24:
            data = ["Cr", 24, 51.996, 0.723 * 10 ** (-5), 7140, 160, 115, 279, 6625, 4013]
        elif self.name in ["Mn", "Manganese", "manganese"] or self.atomicnumber == 25:
            data = ["Mn", 25, 54.938, 0.735 * 10 ** (-5), 7440, 120, 90, 198, 5150, 3477]
        elif self.name in ["Fe", "Iron", "iron"] or self.atomicnumber == 26:
            data = ["Fe", 26, 55.845, 7.09 * 10 ** (-6), 7874, 170, 82, 211, 5956, 3227, 10*10**(-8)]
        elif self.name in ["Co", "Cobalt", "cobalt"] or self.atomicnumber == 27:
            data = ["Co", 27, 58.933, 6.67*10**(-6), 8890, 180, 209, 75, 7183, 4849, 6*10**(-8)]
        elif self.name in ["Cu", "Copper", "copper"] or self.atomicnumber == 29:
            data = ["Cu", 29, 63.546, 0.711 * 10 ** (-5), 8920, 140, 48, 130, 4442, 2320]
        elif self.name in ["As", "Arsenic", "arsenic"] or self.atomicnumber == 33:
            data = ["As", 33, 74.922, 1.295 * 10 ** (-5), 5727, 22, 3, 8, 2131, 724]
        elif self.name in ["Nb", "Niobium", "niobium"] or self.atomicnumber == 41:
            data = ["Nb", 41, 92.906, 10.83*10**(-6), 8580, 170, 38, 105, 5071, 2104]
        elif self.name in ["Mo", "Molybdenum", "molybdenum"] or self.atomicnumber == 42:
            data = ["Mo", 42, 95.95, 0.938 * 10 ** (-5), 10280, 230, 20, 329, 4882, 1395]
        elif self.name in ["Ag", "Silver", "silver"] or self.atomicnumber == 47:
            data = ["Ag", 47, 107.87, 1.027 * 10 ** (-5), 10490, 100, 30, 83, 3653, 1691]
        elif self.name in ["Sb", "Antimony", "antimony"] or self.atomicnumber == 51:
            data = ["Sb", 51, 121.76, 1.819 * 10 ** (-5), 6690, 42, 20, 55, 3202, 1728]
        elif self.name in ["Ba", "Barium", "barium"] or self.atomicnumber == 56:
            data = ["Ba", 56, 137.33, 38.16*10**(-6), 3650, 9.6, 4.9, 13, 2102, 1159, 34*10**(-8)]
        elif self.name in ["Ta", "Tantalum", "tantalum"] or self.atomicnumber == 73:
            data = ["Ta", 73, 180.95, 10.85*10**(-6), 16680, 200, 69, 186, 4184, 2034]
        elif self.name in ["W", "Tungsten", "tungsten"] or self.atomicnumber == 74:
            data = ["W", 74, 183.84, 9.47 * 10 ** (-6), 19260, 304, 147, 380, 5095, 2763]
        elif self.name in ["Au", "Gold", "gold"] or self.atomicnumber == 79:
            data = ["Au", 79, 196.97, 10.21 * 10 ** (-6), 19320, 137, 15, 43, 2851, 881, 2.2*10**(-8)]
        elif self.name in ["Pb", "Lead", "lead"] or self.atomicnumber == 82:
            data = ["Pb", 82, 207.2, 1.826 * 10 ** (-5), 11340, 46, 5.75, 16, 2106, 712]
        elif self.name in ["Bi", "Bismuth", "bismuth"] or self.atomicnumber == 83:
            data = ["Bi", 83, 208.98, 2.131 * 10 ** (-5), 9780, 31, 12, 32, 2192, 1108]
        elif self.name in ["U", "Uranium", "uranium"] or self.atomicnumber == 92:
            data = ["U", 92, 238.03, 12.49*10**(-6), 18970, 100, 111, 207, 3616, 2419, 28*10**(-8)]
        #
        return data