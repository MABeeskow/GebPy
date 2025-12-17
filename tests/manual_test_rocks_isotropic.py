#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_rocks_isotropic.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: manual_test_rocks_isotropic.py
Manual test file related to module isotropic_rocks.py
"""

# PACKAGES
import time
import pandas as pd

# MODULES
from src.gebpy.core.rocks.isotropic_rocks import IsotropicRocks
from gebpy_legacy.modules.siliciclastics import SiliciclasticRocks

# CODE
n_datasets = 75
print("\n--- Manual test for: sedimentary.py ---")
print(f"\nDEFAULT_DATA (SANDSTONE, WATER):")
start = time.time()
DEFAULT_DATA = IsotropicRocks(name="Sandstone", random_seed=42).generate_dataset(number=n_datasets)
end = time.time()
delta_new = end - start
print(f"Runtime: {delta_new:.5f} seconds")
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)
if n_datasets < 20:
    print("Results:", DEFAULT_DATA.describe())

print(f"\nOLD_DATA (SANDSTONE, WATER):")
start = time.time()
OLD_DATA = data = SiliciclasticRocks(fluid="water", actualThickness=0).create_sandstone(
    number=n_datasets, porosity=[0.1, 0.35])
end = time.time()
delta_old = end - start
print(f"Runtime: {delta_old:.5f} seconds")

if n_datasets < 20:
    print("Results:", OLD_DATA)

speed_ratio = delta_old/delta_new
speed_boost = (delta_old - delta_new)/delta_old*100
print("\nSpeed ratio:", round(speed_ratio, 4))
print("Speed boost:", round(speed_boost, 2), "%")

n_datasets = 15

print(f"\nDATA(SANDSTONE, WATER):")
data_rock = IsotropicRocks(name="Sandstone", random_seed=42).generate_dataset(number=n_datasets, fluid="water")
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(SANDSTONE, OIL):")
data_rock = IsotropicRocks(name="Sandstone", random_seed=42).generate_dataset(number=n_datasets, fluid="oil")
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(SANDSTONE, NATURAL GAS):")
data_rock = IsotropicRocks(name="Sandstone", random_seed=42).generate_dataset(number=n_datasets, fluid="natural gas")
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(SANDSTONE, CUSTOM FLUID DENSITY):")
data_rock = IsotropicRocks(name="Sandstone", random_seed=42).generate_dataset(number=n_datasets, density_fluid=725)
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(LIMESTONE, WATER):")
data_rock = IsotropicRocks(name="Limestone", random_seed=42).generate_dataset(number=n_datasets)
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(LIMESTONE, OIL):")
data_rock = IsotropicRocks(name="Limestone", random_seed=42).generate_dataset(number=n_datasets, fluid="oil")
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(LIMESTONE, NATURAL GAS):")
data_rock = IsotropicRocks(name="Limestone", random_seed=42).generate_dataset(number=n_datasets, fluid="natural gas")
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(LIMESTONE, CUSTOM FLUID DENSITY):")
data_rock = IsotropicRocks(name="Limestone", random_seed=42).generate_dataset(number=n_datasets, density_fluid=725)
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(DOLOSTONE, WATER):")
data_rock = IsotropicRocks(name="Dolostone", random_seed=42).generate_dataset(number=n_datasets)
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(MARL, WATER):")
data_rock = IsotropicRocks(name="Marl", random_seed=42, variability=True).generate_dataset(number=n_datasets)
if n_datasets < 20:
    print("Results:", data_rock.describe())