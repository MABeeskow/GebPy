#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_sulfides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		17.12.2025

#-----------------------------------------------

"""
Module: manual_test_sulfides.py
Manual test file related to module sulfides.py
"""

# PACKAGES
import time

from numpy.core.defchararray import lower

# MODULES
from src.gebpy.core.minerals.sulfides import Sulfides
from gebpy_legacy.modules.sulfides import Sulfides as Sulfides_old

# CODE
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