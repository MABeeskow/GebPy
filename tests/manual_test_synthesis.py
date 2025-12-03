#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_synthesis.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		03.12.2025

#-----------------------------------------------

"""
Module: manual_test_synthesis.py
Manual test file related to synthesis.py
"""

# PACKAGES
import time
import psutil, os
import pandas as pd
pd.set_option("display.max_columns", 10)
pd.set_option("display.max_rows", 10)

# MODULES
from src.gebpy.core.minerals.synthesis import MineralDataGeneration, DEFAULT_DATA

print("\n--- Manual test for: synthesis.py ---")
print(f"\nDEFAULT_DATA:")
print(
    f"Columns in dataset: {list(DEFAULT_DATA.columns)}\n"
    f"Number of datapoints: {len(DEFAULT_DATA)}\n"
    f"Complete dataset: {DEFAULT_DATA.head(5)}\n"
)

data_init = MineralDataGeneration("Biotite", 15)
print(f"Initialized: {data_init}")
print(
    f"Mineral name: {data_init.name}\n"
    f"Datapoints: {data_init.n_datapoints}\n"
    f"Random seed: {data_init.random_seed}\n"
)


try:
    start = time.time()
    data_exp = data_init.generate_data()
    end = time.time()
    print(f"Runtime: {end - start:.3f} seconds")
    print(f"\nGenerated dataset for {data_init.name}:")
    print(f"Columns in dataset: {list(data_exp.columns)}")
    print(f"Number of datapoints: {len(data_exp)}")
    print(f"\nDataset preview for {data_init.name}:")
    print(data_exp.head(5))
    print(f"\nDataset statistics:\n{data_exp.describe()}\n")
except Exception as e:
    print(f"\n❌ Error while generating data for {data_init.name}: {e}")

data_init = MineralDataGeneration("Plagioclase", 15)
print(f"Initialized: {data_init}")
print(
    f"Mineral name: {data_init.name}\n"
    f"Datapoints: {data_init.n_datapoints}\n"
    f"Random seed: {data_init.random_seed}\n"
)


try:
    start = time.time()
    data_exp = data_init.generate_data()
    end = time.time()
    print(f"Runtime: {end - start:.3f} seconds")
    print(f"\nGenerated dataset for {data_init.name}:")
    print(f"Columns in dataset: {list(data_exp.columns)}")
    print(f"Number of datapoints: {len(data_exp)}")
    print(f"\nDataset preview for {data_init.name}:")
    print(data_exp.head(5))
    print(f"\nDataset statistics:\n{data_exp.describe()}\n")
except Exception as e:
    print(f"\n❌ Error while generating data for {data_init.name}: {e}")

data_init = MineralDataGeneration("Coltan", 15)
print(f"Initialized: {data_init}")
print(
    f"Mineral name: {data_init.name}\n"
    f"Datapoints: {data_init.n_datapoints}\n"
    f"Random seed: {data_init.random_seed}\n"
)


try:
    start = time.time()
    data_exp = data_init.generate_data()
    end = time.time()
    print(f"Runtime: {end - start:.3f} seconds")
    print(f"\nGenerated dataset for {data_init.name}:")
    print(f"Columns in dataset: {list(data_exp.columns)}")
    print(f"Number of datapoints: {len(data_exp)}")
    print(f"\nDataset preview for {data_init.name}:")
    print(data_exp.head(5))
    print(f"\nDataset statistics:\n{data_exp.describe()}\n")
except Exception as e:
    print(f"\n❌ Error while generating data for {data_init.name}: {e}")

data_init = MineralDataGeneration("Beeskowite", 15)
print(f"Initialized: {data_init}")
print(
    f"Mineral name: {data_init.name}\n"
    f"Datapoints: {data_init.n_datapoints}\n"
    f"Random seed: {data_init.random_seed}\n"
)

try:
    start = time.time()
    data_exp = data_init.generate_data()
    end = time.time()
    print(f"Runtime: {end - start:.3f} seconds")
    print(f"\nGenerated dataset for {data_init.name}:")
    print(f"Columns in dataset: {list(data_exp.columns)}")
    print(f"Number of datapoints: {len(data_exp)}")
    print(f"\nDataset preview for {data_init.name}:")
    print(data_exp.head(5))
    print(f"\nDataset statistics:\n{data_exp.describe()}\n")
except Exception as e:
    print(f"❌ Error while generating data for {data_init.name}: {e}")

process = psutil.Process(os.getpid())
mem = process.memory_info().rss / (1024 ** 3)
print(f"\nMemory used: {mem:.3f} GB")

print("\n--- Test completed successfully ---")