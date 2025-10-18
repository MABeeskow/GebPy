#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_synthesis.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		18.10.2025

#-----------------------------------------------

"""
Module: manual_test_synthesis.py
Manual test file related to synthesis.py
"""

# MODULES
from src.gebpy.core.minerals.synthesis import MineralDataGeneration, DEFAULT_DATA

print("DEFAULT_DATA:", DEFAULT_DATA)

data_init = MineralDataGeneration("Biotite", 15)
data_exp = data_init.generate_data()
print("Created:", data_exp)