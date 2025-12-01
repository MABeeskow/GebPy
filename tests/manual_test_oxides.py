#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_oxides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		02.12.2025

#-----------------------------------------------

"""
Module: manual_test_oxides.py
Manual test file related to module oxides.py
"""

# PACKAGES
import time

from numpy.core.defchararray import lower

# MODULES
from src.gebpy.core.minerals.oxides import Oxides
from gebpy.modules.oxides import Oxides as Oxides_old

# CODE
n_datasets = 10
print("\n--- Manual test for: oxides.py ---")
print(f"\nDEFAULT_DATA (QUARTZ):")
start = time.time()
DEFAULT_DATA = Oxides(name="Quartz", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nOLD_DATA (QUARTZ):")
start = time.time()
OLD_DATA = Oxides_old(mineral="Quartz", data_type=True, traces_list=[]).generate_dataset(number=n_datasets)
end = time.time()
delta_old = end - start
print(f"Runtime: {delta_old:.5f} seconds")

if n_datasets < 20:
    print("Results:", OLD_DATA)

speed_ratio = delta_old/delta_new
try:
    speed_boost = (delta_old - delta_new)/delta_old*100
except:
    speed_boost = (delta_old - delta_new)/0.000001*100
print("\nSpeed ratio:", round(speed_ratio, 4))
print("Speed boost:", round(speed_boost, 2), "%")

n_datasets = 5

print(f"\nDATA (ANATASE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Anatase", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MANGANITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Manganite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (GROUTITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Groutite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (PYROPHANITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Pyrophanite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (GEIKIELITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Geikielite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ESKOLAITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Eskolaite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (KARELINAITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Karelianite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CLAUDETITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Claudetite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ARSENOLITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Arsenolite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (VALENTINITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Valentinite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SENARMONTITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Senarmontite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (BISMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Bismite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (PYROLUSITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Pyrolusite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (BRUCITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Brucite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (GALAXITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Galaxite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (GAHNITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Gahnite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (HERCYNITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Hercynite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SPINEL):")
start = time.time()
DEFAULT_DATA = Oxides(name="Spinel", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (AL-SPINEL):")
start = time.time()
DEFAULT_DATA = Oxides(name="Al-Spinel", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CHROMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Chromite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ZINCOCHROMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Zincochromite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MAGNESIOCHROMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Magnesiochromite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CR-SPINEL):")
start = time.time()
DEFAULT_DATA = Oxides(name="Cr-Spinel", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CASSITERITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Cassiterite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CUPROSPINEL):")
start = time.time()
DEFAULT_DATA = Oxides(name="Cuprospinel", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (JACOBSITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Jacobsite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MAGNESIOFERRITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Magnesioferrite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (TREVORITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Trevorite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (FRANKLINITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Franklinite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MAGNETITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Magnetite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (FE-SPINEL):")
start = time.time()
DEFAULT_DATA = Oxides(name="Fe-Spinel", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (HUEBNERITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Huebnerite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (FERBERITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Ferberite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (WOLFRAMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Wolframite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (URANINITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Uraninite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (HEMATITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Hematite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CORUNDUM):")
start = time.time()
DEFAULT_DATA = Oxides(name="Corundum", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (TISTARITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Tistarite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CORUNDUM-GROUP):")
start = time.time()
DEFAULT_DATA = Oxides(name="Corundum-Group", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (NICHROMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Nichromite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MANGANOCHROMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Manganochromite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (COCHROMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Cochromite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CHROMITE-GROUP):")
start = time.time()
DEFAULT_DATA = Oxides(name="Chromite-Group", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (BOEHMITE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Boehmite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (DIASPORE):")
start = time.time()
DEFAULT_DATA = Oxides(name="Diaspore", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)