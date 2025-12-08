#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_carbonates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		08.12.2025

#-----------------------------------------------

"""
Module: manual_test_carbonates.py
Manual test file related to module carbonates.py
"""

# PACKAGES
import time

from numpy.core.defchararray import lower

# MODULES
from src.gebpy.core.minerals.carbonates import Carbonates
from gebpy.modules.carbonates import Carbonates as Carbonates_old

# CODE
n_datasets = 10
print("\n--- Manual test for: carbonates.py ---")
print(f"\nDEFAULT_DATA (CALCITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Calcite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nOLD_DATA (CALCITE):")
start = time.time()
OLD_DATA = Carbonates_old(mineral="Calcite", data_type=True, traces_list=[]).generate_dataset(number=n_datasets)
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

print(f"\nDATA (DOLOMITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Dolomite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MAGNESITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Magnesite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SIDERITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Siderite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (RHODOCHROSITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Rhodochrosite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ARAGONITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Aragonite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (CERUSSITE):")
start = time.time()
DEFAULT_DATA = Carbonates(name="Cerussite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)