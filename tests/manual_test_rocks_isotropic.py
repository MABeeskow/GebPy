#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_rocks_isotropic.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.01.2026

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

print(f"\nDATA(MARL, WATER, CONSTRAINED ELEMENT AMOUNT):")
data_rock = IsotropicRocks(name="Marl", random_seed=42, variability=True).generate_dataset(
    number=n_datasets, element_constraints={"C": (0.10, 0.11)})
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(MARL, WATER, CONSTRAINED MINERAL COMPOSITION):")
data_rock = IsotropicRocks(name="Marl", random_seed=42, variability=True).generate_dataset(
    number=n_datasets, mineral_comp={
        "Calcite": 0.8, "Dolomite": 0.05, "Illite": 0.05, "Chlorite": 0.0, "Kaolinite": 0.0, "Montmorillonite": 0.0,
        "Quartz": 0.05, "Plagioclase": 0.05})
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(LIMESTONE, WATER, ADDITIONAL MINERAL ASSEMBLAGE):")
data_rock = IsotropicRocks(name="Limestone", random_seed=42).generate_dataset(
    number=n_datasets, additional_assemblage={"volume_fraction": 0.15, "mineralogy": {
        "Galena": 0.4, "Sphalerite": 0.5, "Chalcopyrite": 0.1}, "rescaling_host": True})
if n_datasets < 20:
    print("Results:", data_rock.describe())

print(f"\nDATA(CUSTOM SANDSTONE, WATER):")
data_rock = IsotropicRocks(name="Sandstone XY", random_seed=42).generate_dataset(
    number=n_datasets, comp_constrained=True, mineral_constrained={
        "Quartz": [0.8, 1.0], "Alkali feldspar": [0.0, 0.2], "Calcite": [0.0, 0.1], "Pyrite": [0.0, 0.05],
        "Sphalerite": [0.0, 0.1], "Galena": [0.0, 0.1]}, porosity_constrained={"min": 0.1, "max": 0.25})
if n_datasets < 20:
    print("Results:", data_rock.describe())

# # TEST
# print(f"\nTEST(YAML ARCHITECTURE):")
# IsotropicRocks._yaml_cache.clear()
# IsotropicRocks._mineralogy_cache.clear()
# IsotropicRocks._mineral_groups_cache.clear()
# print("=== First call ===")
# rock1 = IsotropicRocks(name="Sandstone", random_seed=32)
# data1 = rock1.generate_dataset(number=5)
# print("YAML cache after first call:",
#       IsotropicRocks._yaml_cache.keys())
#
# print("=== Second call (same rock) ===")
# rock2 = IsotropicRocks(name="Sandstone", random_seed=64)
# data2 = rock2.generate_dataset(number=5)
# print("YAML cache after second call:",
#       IsotropicRocks._yaml_cache.keys())

# print(f"\nTEST(MINERAL SAMPLING):")
# print("=== TEST 1: Basic sampling ===")
# rock = IsotropicRocks(name="Sandstone", random_seed=42)
# data = rock.generate_dataset(number=10)
# print(data.head())
# print("Shape:", data.shape)
#
# print("\n=== TEST 2: Mineral fractions sum to 1 ===")
# phi_cols = [c for c in data.columns if c.startswith("phi.")]
# row_sums = data[phi_cols].sum(axis=1)
# print("Row sums:")
# print(row_sums)
# print("Min:", row_sums.min(), "Max:", row_sums.max())
# print("\n=== TEST 3: Reproducibility ===")
# r1 = IsotropicRocks("Sandstone", random_seed=123)
# r2 = IsotropicRocks("Sandstone", random_seed=123)
# d1 = r1.generate_dataset(number=5)
# d2 = r2.generate_dataset(number=5)
# phi_cols = [c for c in d1.columns if c.startswith("phi.")]
# print("Difference in mineral fractions:")
# print((d1[phi_cols] - d2[phi_cols]).abs().max())
# print("\n=== TEST 4: Mineral bounds ===")
# for col in phi_cols:
#     print(
#         col,
#         "min =", data[col].min(),
#         "max =", data[col].max()
#     )
# print("\n=== TEST 5: Element constraints ===")
# constraints = {"Si": (0.2, 0.4)}
# rock_c = IsotropicRocks("Sandstone", random_seed=7)
# data_c = rock_c.generate_dataset(
#     number=5,
#     element_constraints=constraints
# )
# print(data_c[["w.Si"]])
# print("Min Si:", data_c["w.Si"].min())
# print("Max Si:", data_c["w.Si"].max())


