#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_sulfides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		28.01.2026

#-----------------------------------------------

"""
Module: manual_test_sulfides.py
Manual test file related to module sulfides.py
"""

# PACKAGES
import time

# MODULES
from src.gebpy.core.minerals.sulfides import Sulfides
from gebpy_legacy.modules.sulfides import Sulfides as Sulfides_old

# CODE
def benchmark_sulfide(name, n_datasets, seed=42, verbose=True):
    print(f"\nDATA ({name.upper()}):")
    start = time.perf_counter()
    data = Sulfides(name=name, random_seed=seed).generate_dataset(number=n_datasets)
    delta = time.perf_counter() - start
    print(f"Runtime: {delta:.5f} seconds")

    if verbose and n_datasets < 20:
        print("Results:", data)

    return delta, data

n_datasets = 10
print("\n--- Manual test for: sulfides.py ---")
print(f"\nDEFAULT_DATA (PYRITE):")
start = time.time()
DEFAULT_DATA = Sulfides(name="Pyrite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nOLD_DATA (PYRITE):")
start = time.time()
OLD_DATA = Sulfides_old(mineral="Pyrite", data_type=True, traces_list=[]).generate_dataset(number=n_datasets)
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

minerals = [
    "Pyrite", "Acanthite", "Bornite", "Cattierite", "Chalcocite", "Galena", "Sphalerite", "Cinnabar", "Cobaltite",
    "Covellite", "Tetrahedrite", "Tennantite", "Fahlore", "Gallite", "Laforetite", "Lenaite", "Roquesite",
    "Chalcopyrite-Group", "Marcasite"]

results = {}
for mineral in sorted(minerals):
    dt, data = benchmark_sulfide(mineral, n_datasets)
    results[mineral] = dt