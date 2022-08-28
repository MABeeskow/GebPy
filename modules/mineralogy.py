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
        if self.keyword == "Spinel Group":
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
                    list_elements.extend(list_elements_al)
                    list_elements.extend(list_elements_fe)
                    list_elements.extend(list_elements_cr)
                    list_elements = list(dict.fromkeys(list_elements))
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
                            if element in data_al_sp[category]:
                                data_minerals["Al-Spi"][category][element].append(round(data_al_sp[category][element], 4))
                            else:
                                data_minerals["Al-Spi"][category][element].append(round(0.0, 4))
                            if element in data_fe_sp[category]:
                                data_minerals["Fe-Spi"][category][element].append(round(data_fe_sp[category][element], 4))
                            else:
                                data_minerals["Fe-Spi"][category][element].append(round(0.0, 4))
                            if element in data_cr_sp[category]:
                                data_minerals["Cr-Spi"][category][element].append(round(data_cr_sp[category][element], 4))
                            else:
                                data_minerals["Cr-Spi"][category][element].append(round(0.0, 4))
                #
                n += 1
            #
            data_minerals["Al-Spi"]["color"] = "tab:blue"
            data_minerals["Fe-Spi"]["color"] = "tab:red"
            data_minerals["Cr-Spi"]["color"] = "tab:green"
        #
        elif self.keyword == "Hematite Group":
            data_minerals = {"Hem": {}, "Crn": {}, "Ilm": {}}
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
                data_hem = Oxides(impurity="pure", data_type=True).create_hematite()
                data_crn = Oxides(impurity="pure", data_type=True).create_corundum()
                data_ilm = Oxides(impurity="pure", data_type=True).create_ilmenite()
                #
                if n == 0:
                    list_elements_hem = list(data_hem["chemistry"].keys())
                    list_elements_crn = list(data_crn["chemistry"].keys())
                    list_elements_ilm = list(data_ilm["chemistry"].keys())
                    list_elements.extend(list_elements_hem)
                    list_elements.extend(list_elements_crn)
                    list_elements.extend(list_elements_ilm)
                    list_elements = list(dict.fromkeys(list_elements))
                    for element in list_elements:
                        data_minerals["Hem"]["chemistry"][element] = []
                        data_minerals["Crn"]["chemistry"][element] = []
                        data_minerals["Ilm"]["chemistry"][element] = []
                #
                for category in categories:
                    if category != "chemistry":
                        data_minerals["Hem"][category].append(round(data_hem[category], 4))
                        data_minerals["Crn"][category].append(round(data_crn[category], 4))
                        data_minerals["Ilm"][category].append(round(data_ilm[category], 4))
                    else:
                        for element in list_elements:
                            if element in data_hem[category]:
                                data_minerals["Hem"][category][element].append(round(data_hem[category][element], 4))
                            else:
                                data_minerals["Hem"][category][element].append(round(0.0, 4))
                            if element in data_crn[category]:
                                data_minerals["Crn"][category][element].append(round(data_crn[category][element], 4))
                            else:
                                data_minerals["Crn"][category][element].append(round(0.0, 4))
                            if element in data_ilm[category]:
                                data_minerals["Ilm"][category][element].append(round(data_ilm[category][element], 4))
                            else:
                                data_minerals["Ilm"][category][element].append(round(0.0, 4))
                #
                n += 1
            #
            data_minerals["Hem"]["color"] = "tab:blue"
            data_minerals["Crn"]["color"] = "tab:red"
            data_minerals["Ilm"]["color"] = "tab:green"
        #
        return data_minerals