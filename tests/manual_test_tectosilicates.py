#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_tectosilicates.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: manual_test_tectosilicates.py
Manual test file related to module phyllosilicates.py
"""

# PACKAGES
import time

# MODULES
from src.gebpy.core.minerals.tectosilicates import Tectosilicates
from gebpy_legacy.modules.silicates import Tectosilicates as Tectosilicates_old

# CODE
n_datasets = 10
print("\n--- Manual test for: tectosilicates.py ---")
print(f"\nDEFAULT_DATA (ORTHOCLASE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Orthoclase", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nOLD_DATA (ORTHOCLASE):")
start = time.time()
OLD_DATA = Tectosilicates_old(mineral="Orthoclase", data_type=True, traces_list=[]).generate_dataset(number=n_datasets)
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

print(f"\nDATA (ALBITE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Albite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MICROCLINE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Microcline", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ANORTHITE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Anorthite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (ALKALI FELDSPAR):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Alkali feldspar", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (PLAGIOCLASE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Plagioclase", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MARIALITE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Marialite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (MEIONITE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Meionite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (SCAPOLITE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Scapolite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (KALSILITE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Kalsilite", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (NA-NEPHELINE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="NaNepheline", random_seed=42, variability=True).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)

print(f"\nDATA (NEPHELINE):")
start = time.time()
DEFAULT_DATA = Tectosilicates(name="Nepheline", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")

if n_datasets < 20:
    print("Results:", DEFAULT_DATA)