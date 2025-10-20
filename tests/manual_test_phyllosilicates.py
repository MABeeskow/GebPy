#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_phyllosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		21.10.2025

#-----------------------------------------------

"""
Module: manual_test_phyllosilicates.py
Manual test file related to module phyllosilicates.py
"""

# PACKAGES
import time

# MODULES
from src.gebpy.core.minerals.phyllosilicates import Phyllosilicates
from gebpy.modules.silicates import Phyllosilicates as Phyllosilicates_old

# CODE
n_datasets = 10
print("\n--- Manual test for: phyllosilicates.py ---")
print(f"\nDEFAULT_DATA (ANNITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Annite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.3f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nOLD_DATA (ANNITE):")
start = time.time()
OLD_DATA = Phyllosilicates_old(mineral="Annite", data_type=True, traces_list=[]).generate_dataset(number=n_datasets)
end = time.time()
delta_old = end - start
print(f"Runtime: {delta_old:.3f} seconds")

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
print(f"Runtime: {delta_new:.3f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SIDEROPHYLLITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Siderophyllite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.3f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (EASTONITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Eastonite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.3f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ILLITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Illite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.3f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (KAOLINITE):")
start = time.time()
DEFAULT_DATA = Phyllosilicates(name="Kaolinite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.3f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)