#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_chemistry.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		26.05.2021

# -----------------------------------------------

## MODULES
from modules import chemistry

class TestingChemistry():
    #
    def __init__(self, numbers_list=None, name=None):
        """
        :param atomicnumber: list of atomic numbers
        """
        self.numbers_list = numbers_list
        self.name = name
        #
        if self.numbers_list != None:
            for i in self.numbers_list:
                print(chemistry.PeriodicSystem(atomicnumber=i).get_data())
        elif self.name != None:
            print(chemistry.PeriodicSystem(name=self.name).get_data())

# RUN
numbers_list = range(1, 99)
TestingChemistry(numbers_list=numbers_list)
TestingChemistry(name="O")