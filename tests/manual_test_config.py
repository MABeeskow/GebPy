#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_config.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		18.10.2025

#-----------------------------------------------

"""
Module: manual_test_config.py
Manual test file related to config.py
"""

# MODULES
from src.gebpy.core.minerals.config import MineralConfiguration, DEFAULT_CONFIG

print("DEFAULT_CONFIG:", DEFAULT_CONFIG)

cfg = MineralConfiguration("Olivine", 500)
print("Created:", cfg)

cfg.set_random_seed(None)
print("After removing seed:", cfg)

try:
    cfg.set_number_of_datapoints(-5)
except ValueError as e:
    print("Caught expected error:", e)