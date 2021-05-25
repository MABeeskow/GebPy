#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_chemistry.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.05.2021

# -----------------------------------------------

## MODULES
from modules import chemistry

class TestingChemistry():
    #
    def __init__(self, numbers_list=[1]):
        """
        :param atomicnumber: list of atomic numbers
        """
        self.numbers_list = numbers_list
        #
        for i in self.numbers_list:
            print(chemistry.PeriodicSystem(atomicnumber=i).get_data())

# RUN
numbers_list = range(1, 7)
TestingChemistry(numbers_list=numbers_list)