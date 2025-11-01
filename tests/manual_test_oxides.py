#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_oxides.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		01.11.2025

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
speed_boost = (delta_old - delta_new)/delta_old*100
print("\nSpeed ratio:", round(speed_ratio, 4))
print("Speed boost:", round(speed_boost, 2), "%")

n_datasets = 5