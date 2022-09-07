#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		mineralogy.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		08.09.2022

#-----------------------------------------------

## MODULES
from modules.oxides import Oxides, RutileGroup, PericlaseGroup, WulfeniteGroup

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
        elif self.keyword == "Rutile Group":
            data_minerals = {"Argt": {}, "Cst": {}, "Prtl": {}, "Pltn": {}, "Prl": {}, "Rt": {}, "Stv": {},
                             "Rt-Group": {}}
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
            data_argt = Oxides(impurity="pure", data_type=True).create_argutite()
            data_cst = Oxides(impurity="pure", data_type=True).create_cassiterite()
            data_prtl = Oxides(impurity="pure", data_type=True).create_paratellurite()
            data_pltn = Oxides(impurity="pure", data_type=True).create_plattnerite()
            data_prl = Oxides(impurity="pure", data_type=True).create_pyrolusite()
            data_rt = Oxides(impurity="pure", data_type=True).create_rutile()
            data_stv = Oxides(impurity="pure", data_type=True).create_stishovite()
            #
            while n < number:
                data_group = RutileGroup().create_rutile_group()
                #
                if n == 0:
                    list_elements_argt = list(data_argt["chemistry"].keys())
                    list_elements_cst = list(data_cst["chemistry"].keys())
                    list_elements_prtl = list(data_prtl["chemistry"].keys())
                    list_elements_pltn = list(data_pltn["chemistry"].keys())
                    list_elements_prl = list(data_prl["chemistry"].keys())
                    list_elements_rt = list(data_rt["chemistry"].keys())
                    list_elements_stv = list(data_stv["chemistry"].keys())
                    list_elements_group = list(data_group["chemistry"].keys())
                    list_elements.extend(list_elements_argt)
                    list_elements.extend(list_elements_cst)
                    list_elements.extend(list_elements_prtl)
                    list_elements.extend(list_elements_pltn)
                    list_elements.extend(list_elements_prl)
                    list_elements.extend(list_elements_rt)
                    list_elements.extend(list_elements_stv)
                    list_elements.extend(list_elements_group)
                    list_elements = list(dict.fromkeys(list_elements))
                    for element in list_elements:
                        data_minerals["Argt"]["chemistry"][element] = []
                        data_minerals["Cst"]["chemistry"][element] = []
                        data_minerals["Prtl"]["chemistry"][element] = []
                        data_minerals["Pltn"]["chemistry"][element] = []
                        data_minerals["Prl"]["chemistry"][element] = []
                        data_minerals["Rt"]["chemistry"][element] = []
                        data_minerals["Stv"]["chemistry"][element] = []
                        data_minerals["Rt-Group"]["chemistry"][element] = []
                #
                for category in categories:
                    if category != "chemistry":
                        data_minerals["Argt"][category].append(round(data_argt[category], 4))
                        data_minerals["Cst"][category].append(round(data_cst[category], 4))
                        data_minerals["Prtl"][category].append(round(data_prtl[category], 4))
                        data_minerals["Pltn"][category].append(round(data_pltn[category], 4))
                        data_minerals["Prl"][category].append(round(data_prl[category], 4))
                        data_minerals["Rt"][category].append(round(data_rt[category], 4))
                        data_minerals["Stv"][category].append(round(data_stv[category], 4))
                        data_minerals["Rt-Group"][category].append(round(data_group[category], 4))
                    else:
                        for element in list_elements:
                            if element in data_argt[category]:
                                data_minerals["Argt"][category][element].append(round(data_argt[category][element], 4))
                            else:
                                data_minerals["Argt"][category][element].append(round(0.0, 4))
                            if element in data_cst[category]:
                                data_minerals["Cst"][category][element].append(round(data_cst[category][element], 4))
                            else:
                                data_minerals["Cst"][category][element].append(round(0.0, 4))
                            if element in data_prtl[category]:
                                data_minerals["Prtl"][category][element].append(round(data_prtl[category][element], 4))
                            else:
                                data_minerals["Prtl"][category][element].append(round(0.0, 4))
                            if element in data_pltn[category]:
                                data_minerals["Pltn"][category][element].append(round(data_pltn[category][element], 4))
                            else:
                                data_minerals["Pltn"][category][element].append(round(0.0, 4))
                            if element in data_prl[category]:
                                data_minerals["Prl"][category][element].append(round(data_prl[category][element], 4))
                            else:
                                data_minerals["Prl"][category][element].append(round(0.0, 4))
                            if element in data_rt[category]:
                                data_minerals["Rt"][category][element].append(round(data_rt[category][element], 4))
                            else:
                                data_minerals["Rt"][category][element].append(round(0.0, 4))
                            if element in data_stv[category]:
                                data_minerals["Stv"][category][element].append(round(data_stv[category][element], 4))
                            else:
                                data_minerals["Stv"][category][element].append(round(0.0, 4))
                            if element in data_group[category]:
                                data_minerals["Rt-Group"][category][element].append(round(data_group[category][element], 4))
                            else:
                                data_minerals["Rt-Group"][category][element].append(round(0.0, 4))
                #
                n += 1
            #
            data_minerals["Argt"]["color"] = "tab:blue"
            data_minerals["Cst"]["color"] = "tab:red"
            data_minerals["Prtl"]["color"] = "tab:green"
            data_minerals["Pltn"]["color"] = "tab:orange"
            data_minerals["Prl"]["color"] = "tab:purple"
            data_minerals["Rt"]["color"] = "tab:brown"
            data_minerals["Stv"]["color"] = "tab:olive"
            data_minerals["Rt-Group"]["color"] = "tab:gray"
        #
        elif self.keyword == "Periclase Group":
            data_minerals = {"Per": {}, "Wus": {}, "Mns": {}, "Bsn": {}, "Mntp": {}, "Lm": {},
                             "Per-Group": {}}
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
            data_per = Oxides(impurity="pure", data_type=True).create_periclase()
            data_wus = Oxides(impurity="pure", data_type=True).create_wustite()
            data_mns = Oxides(impurity="pure", data_type=True).create_manganosite()
            data_bsn = Oxides(impurity="pure", data_type=True).create_bunsenite()
            data_mntp = Oxides(impurity="pure", data_type=True).create_monteponite()
            data_lm = Oxides(impurity="pure", data_type=True).create_lime()
            #
            while n < number:
                data_group = PericlaseGroup().create_periclase_group()
                #
                if n == 0:
                    list_elements_per = list(data_per["chemistry"].keys())
                    list_elements_wus = list(data_wus["chemistry"].keys())
                    list_elements_mns = list(data_mns["chemistry"].keys())
                    list_elements_bsn = list(data_bsn["chemistry"].keys())
                    list_elements_mntp = list(data_mntp["chemistry"].keys())
                    list_elements_lm = list(data_lm["chemistry"].keys())
                    list_elements_group = list(data_group["chemistry"].keys())
                    list_elements.extend(list_elements_per)
                    list_elements.extend(list_elements_wus)
                    list_elements.extend(list_elements_mns)
                    list_elements.extend(list_elements_bsn)
                    list_elements.extend(list_elements_mntp)
                    list_elements.extend(list_elements_lm)
                    list_elements.extend(list_elements_group)
                    list_elements = list(dict.fromkeys(list_elements))
                    for element in list_elements:
                        data_minerals["Per"]["chemistry"][element] = []
                        data_minerals["Wus"]["chemistry"][element] = []
                        data_minerals["Mns"]["chemistry"][element] = []
                        data_minerals["Bsn"]["chemistry"][element] = []
                        data_minerals["Mntp"]["chemistry"][element] = []
                        data_minerals["Lm"]["chemistry"][element] = []
                        data_minerals["Per-Group"]["chemistry"][element] = []
                #
                for category in categories:
                    if category != "chemistry":
                        data_minerals["Per"][category].append(round(data_per[category], 4))
                        data_minerals["Wus"][category].append(round(data_wus[category], 4))
                        data_minerals["Mns"][category].append(round(data_mns[category], 4))
                        data_minerals["Bsn"][category].append(round(data_bsn[category], 4))
                        data_minerals["Mntp"][category].append(round(data_mntp[category], 4))
                        data_minerals["Lm"][category].append(round(data_lm[category], 4))
                        data_minerals["Per-Group"][category].append(round(data_group[category], 4))
                    else:
                        for element in list_elements:
                            if element in data_per[category]:
                                data_minerals["Per"][category][element].append(round(data_per[category][element], 4))
                            else:
                                data_minerals["Per"][category][element].append(round(0.0, 4))
                            if element in data_wus[category]:
                                data_minerals["Wus"][category][element].append(round(data_wus[category][element], 4))
                            else:
                                data_minerals["Wus"][category][element].append(round(0.0, 4))
                            if element in data_mns[category]:
                                data_minerals["Mns"][category][element].append(round(data_mns[category][element], 4))
                            else:
                                data_minerals["Mns"][category][element].append(round(0.0, 4))
                            if element in data_bsn[category]:
                                data_minerals["Bsn"][category][element].append(round(data_bsn[category][element], 4))
                            else:
                                data_minerals["Bsn"][category][element].append(round(0.0, 4))
                            if element in data_mntp[category]:
                                data_minerals["Mntp"][category][element].append(round(data_mntp[category][element], 4))
                            else:
                                data_minerals["Mntp"][category][element].append(round(0.0, 4))
                            if element in data_lm[category]:
                                data_minerals["Lm"][category][element].append(round(data_lm[category][element], 4))
                            else:
                                data_minerals["Lm"][category][element].append(round(0.0, 4))
                            if element in data_group[category]:
                                data_minerals["Per-Group"][category][element].append(round(data_group[category][element], 4))
                            else:
                                data_minerals["Per-Group"][category][element].append(round(0.0, 4))
                #
                n += 1
            #
            data_minerals["Per"]["color"] = "tab:blue"
            data_minerals["Wus"]["color"] = "tab:red"
            data_minerals["Mns"]["color"] = "tab:green"
            data_minerals["Bsn"]["color"] = "tab:orange"
            data_minerals["Mntp"]["color"] = "tab:purple"
            data_minerals["Lm"]["color"] = "tab:brown"
            data_minerals["Per-Group"]["color"] = "tab:gray"
        #
        elif self.keyword == "Wulfenite Group":
            data_minerals = {"Wul": {}, "Crc": {}, "Sz": {}, "Wul-Group": {}}
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
            data_wul = Oxides(impurity="pure", data_type=True).create_wulfenite()
            data_crc = Oxides(impurity="pure", data_type=True).create_crocoite()
            data_sz = Oxides(impurity="pure", data_type=True).create_stolzite()
            #
            while n < number:
                data_group = WulfeniteGroup().create_wulfenite_group()
                #
                if n == 0:
                    list_elements_wul = list(data_wul["chemistry"].keys())
                    list_elements_crc = list(data_crc["chemistry"].keys())
                    list_elements_sz = list(data_sz["chemistry"].keys())
                    list_elements_group = list(data_group["chemistry"].keys())
                    list_elements.extend(list_elements_wul)
                    list_elements.extend(list_elements_crc)
                    list_elements.extend(list_elements_sz)
                    list_elements.extend(list_elements_group)
                    list_elements = list(dict.fromkeys(list_elements))
                    for element in list_elements:
                        data_minerals["Wul"]["chemistry"][element] = []
                        data_minerals["Crc"]["chemistry"][element] = []
                        data_minerals["Sz"]["chemistry"][element] = []
                        data_minerals["Wul-Group"]["chemistry"][element] = []
                #
                for category in categories:
                    if category != "chemistry":
                        data_minerals["Wul"][category].append(round(data_wul[category], 4))
                        data_minerals["Crc"][category].append(round(data_crc[category], 4))
                        data_minerals["Sz"][category].append(round(data_sz[category], 4))
                        data_minerals["Wul-Group"][category].append(round(data_group[category], 4))
                    else:
                        for element in list_elements:
                            if element in data_wul[category]:
                                data_minerals["Wul"][category][element].append(round(data_wul[category][element], 4))
                            else:
                                data_minerals["Wul"][category][element].append(round(0.0, 4))
                            if element in data_crc[category]:
                                data_minerals["Crc"][category][element].append(round(data_crc[category][element], 4))
                            else:
                                data_minerals["Crc"][category][element].append(round(0.0, 4))
                            if element in data_sz[category]:
                                data_minerals["Sz"][category][element].append(round(data_sz[category][element], 4))
                            else:
                                data_minerals["Sz"][category][element].append(round(0.0, 4))
                            if element in data_group[category]:
                                data_minerals["Wul-Group"][category][element].append(round(data_group[category][element], 4))
                            else:
                                data_minerals["Wul-Group"][category][element].append(round(0.0, 4))
                #
                n += 1
            #
            data_minerals["Wul"]["color"] = "tab:blue"
            data_minerals["Crc"]["color"] = "tab:red"
            data_minerals["Sz"]["color"] = "tab:green"
            data_minerals["Wul-Group"]["color"] = "tab:gray"
        #
        return data_minerals