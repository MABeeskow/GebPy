#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		sedimentary.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		03.12.2025

#-----------------------------------------------

"""
Module: sedimentary.py
This module controls the generation of the synthetic data of the oxides minerals.
"""

# PACKAGES
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

from soupsieve.util import lower
from asteval import Interpreter

# MODULES

# Code
class SedimentaryRocks:

    def __init__(self, name, random_seed, rounding: int = 3) -> None:
        self.name = name
        self.random_seed = random_seed
        self.rounding = rounding

    def _generate_data(self, number: int = 1, as_dataframe=False) -> None:
        siliciclastics = {"Sandstone"}
