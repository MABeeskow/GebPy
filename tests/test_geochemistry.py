#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_geochemistry.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		28.05.2021

# -----------------------------------------------

## MODULES
from modules import geochemistry
from modules.chemistry import PeriodicSystem
from modules.oxides import Quartz

class TestingGeochemistry():

    def __init__(self, name, traces):
        self.name = name
        self.traces = traces
        TestingGeochemistry.test_mass_spectrometry(self)

    def test_mass_spectrometry(self):
        if self.name in ["Qz", "Quartz", "quartz"]:
            mineral = Quartz(traces_list=self.traces).create_quartz()
        else:
            print("No mineral found!")
        geochemistry.MassSpectrometry(data=mineral)

# TESTING
TestingGeochemistry(name="Quartz", traces=["Al", "Ti"])