#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_phyllosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: manual_test_phyllosilicates.py
Manual test file related to module phyllosilicates.py
"""

# PACKAGES
import time

from numpy.core.defchararray import lower

# MODULES
from src.gebpy.core.minerals.phyllosilicates import Phyllosilicates
from gebpy_legacy.modules.silicates import Phyllosilicates as Phyllosilicates_old

# CODE
n_datasets = 10
print("\n--- Manual test for: phyllosilicates.py ---")
print(f"\nDEFAULT_DATA (ANNITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Annite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nOLD_DATA (ANNITE):")
start = time.time()
OLD_DATA = Phyllosilicates_old(mineral="Annite", data_type=True, traces_list=[]).generate_dataset(number=n_datasets)
end = time.time()
delta_old = end - start
print(f"Runtime: {delta_old:.5f} seconds")

if n_datasets < 20:
    print("Results:", OLD_DATA)

speed_ratio = delta_old/delta_new
speed_boost = (delta_old - delta_new)/delta_old*100
print("\nSpeed ratio:", round(speed_ratio, 4))
print("Speed boost:", round(speed_boost, 2), "%")

n_datasets = 5

print(f"\nDATA (PHLOGOPITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Phlogopite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SIDEROPHYLLITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Siderophyllite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (EASTONITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Eastonite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ILLITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Illite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (KAOLINITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Kaolinite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CHAMOSITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Chamosite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CLINOCHLORE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Clinochlore", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (PENNANTITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Pennantite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (NIMITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Nimite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MUSCOVITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Muscovite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (TALC):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Talc", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CHRYSOTILE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Chrysotile", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ANTIGORITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Antigorite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (PYROPHYLLITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Pyrophyllite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (BIOTITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Biotite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MONTMORILLONITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Montmorillonite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (NONTRONITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Nontronite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SAPONITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Saponite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (GLAUCONITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Glauconite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (VERMICULITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Vermiculite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CHLORITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Chlorite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (FE-CHLORITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="FeChlorite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MG-CHLORITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="MgChlorite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MN-CHLORITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="MnChlorite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (NI-CHLORITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="NiChlorite", random_seed=42, variability=True).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)