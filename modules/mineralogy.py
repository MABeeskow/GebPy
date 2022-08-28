#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		mineralogy.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		28.08.2020

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
            data_minerals = {"Hem": {}, "Crn": {}, "Esk": {}, "Tta": {}, "Kar": {}, "Hem-Group": {}}
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
            #
            data_hem = Oxides(impurity="pure", data_type=True).create_hematite()
            data_crn = Oxides(impurity="pure", data_type=True).create_corundum()
            data_esk = Oxides(impurity="pure", data_type=True).create_eskolaite()
            data_tta = Oxides(impurity="pure", data_type=True).create_tistarite()
            data_kar = Oxides(impurity="pure", data_type=True).create_karelianite()
            #
            while n < number:
                data_group = Oxides(impurity="pure", data_type=True).create_hematite_group()
                #
                if n == 0:
                    list_elements_hem = list(data_hem["chemistry"].keys())
                    list_elements_crn = list(data_crn["chemistry"].keys())
                    list_elements_esk = list(data_esk["chemistry"].keys())
                    list_elements_tta = list(data_tta["chemistry"].keys())
                    list_elements_kar = list(data_kar["chemistry"].keys())
                    list_elements_group = list(data_group["chemistry"].keys())
                    list_elements.extend(list_elements_hem)
                    list_elements.extend(list_elements_crn)
                    list_elements.extend(list_elements_esk)
                    list_elements.extend(list_elements_tta)
                    list_elements.extend(list_elements_kar)
                    list_elements.extend(list_elements_group)
                    list_elements = list(dict.fromkeys(list_elements))
                    for element in list_elements:
                        data_minerals["Hem"]["chemistry"][element] = []
                        data_minerals["Crn"]["chemistry"][element] = []
                        data_minerals["Esk"]["chemistry"][element] = []
                        data_minerals["Tta"]["chemistry"][element] = []
                        data_minerals["Kar"]["chemistry"][element] = []
                        data_minerals["Hem-Group"]["chemistry"][element] = []
                #
                for category in categories:
                    if category != "chemistry":
                        data_minerals["Hem"][category].append(round(data_hem[category], 4))
                        data_minerals["Crn"][category].append(round(data_crn[category], 4))
                        data_minerals["Esk"][category].append(round(data_esk[category], 4))
                        data_minerals["Tta"][category].append(round(data_tta[category], 4))
                        data_minerals["Kar"][category].append(round(data_kar[category], 4))
                        data_minerals["Hem-Group"][category].append(round(data_group[category], 4))
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
                            if element in data_esk[category]:
                                data_minerals["Esk"][category][element].append(round(data_esk[category][element], 4))
                            else:
                                data_minerals["Esk"][category][element].append(round(0.0, 4))
                            if element in data_tta[category]:
                                data_minerals["Tta"][category][element].append(round(data_tta[category][element], 4))
                            else:
                                data_minerals["Tta"][category][element].append(round(0.0, 4))
                            if element in data_kar[category]:
                                data_minerals["Kar"][category][element].append(round(data_kar[category][element], 4))
                            else:
                                data_minerals["Kar"][category][element].append(round(0.0, 4))
                            if element in data_group[category]:
                                data_minerals["Hem-Group"][category][element].append(round(data_group[category][element], 4))
                            else:
                                data_minerals["Hem-Group"][category][element].append(round(0.0, 4))
                #
                n += 1
            #
            data_minerals["Hem"]["color"] = "tab:blue"
            data_minerals["Crn"]["color"] = "tab:red"
            data_minerals["Esk"]["color"] = "tab:green"
            data_minerals["Tta"]["color"] = "tab:orange"
            data_minerals["Kar"]["color"] = "tab:purple"
            data_minerals["Hem-Group"]["color"] = "tab:gray"
        #
        return data_minerals