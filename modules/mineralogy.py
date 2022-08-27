#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		mineralogy.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		27.08.2020

#-----------------------------------------------

## MODULES
from modules.oxides import Oxides

class Mineralogy:
    #
    def __init__(self, keyword):
        self.keyword = keyword
    #
    def compare_minerals(self, number=10):
        if self.keyword == "Spinel Minerals":
            data_minerals = {"Al-Spi": {}, "Fe-Spi": {}, "Cr-Spi": {}}
            categories = ["M", "rho", "V", "vP", "vS", "vP/vS", "K", "G", "E", "nu", "GR", "PE", "U", "chemistry"]
            for key, value in data_minerals.items():
                for category in categories:
                    if category != "chemistry":
                        value[category] = []
                    else:
                        value[category] = {}
            #
            n = 0
            list_elements = []
            while n < number:
                data_al_sp = Oxides(impurity="pure", data_type=True).create_aluminium_spinel()
                data_fe_sp = Oxides(impurity="pure", data_type=True).create_iron_spinel()
                data_cr_sp = Oxides(impurity="pure", data_type=True).create_chromium_spinel()
                #
                if n == 0:
                    list_elements_al = list(data_al_sp["chemistry"].keys())
                    list_elements_fe = list(data_fe_sp["chemistry"].keys())
                    list_elements_cr = list(data_cr_sp["chemistry"].keys())
                    list_elements = list(set(list_elements_al).intersection(list_elements_fe))
                    list_elements = list(set(list_elements).intersection(list_elements_cr))
                    list_elements.remove("O")
                    for element in list_elements:
                        data_minerals["Al-Spi"]["chemistry"][element] = []
                        data_minerals["Fe-Spi"]["chemistry"][element] = []
                        data_minerals["Cr-Spi"]["chemistry"][element] = []
                #
                for category in categories:
                    if category != "chemistry":
                        data_minerals["Al-Spi"][category].append(round(data_al_sp[category], 4))
                        data_minerals["Fe-Spi"][category].append(round(data_fe_sp[category], 4))
                        data_minerals["Cr-Spi"][category].append(round(data_cr_sp[category], 4))
                    else:
                        for element in list_elements:
                            data_minerals["Al-Spi"][category][element].append(round(data_al_sp[category][element], 4))
                            data_minerals["Fe-Spi"][category][element].append(round(data_al_sp[category][element], 4))
                            data_minerals["Cr-Spi"][category][element].append(round(data_al_sp[category][element], 4))
                #
                n += 1
            #
            data_minerals["Al-Spi"]["color"] = "tab:blue"
            data_minerals["Fe-Spi"]["color"] = "tab:red"
            data_minerals["Cr-Spi"]["color"] = "tab:green"
        #
        return data_minerals