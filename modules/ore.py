#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		ore.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		25.05.2022

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.silicates import Phyllosilicates
from modules.oxides import Oxides
from modules.carbonates import Carbonates
from modules.sulfides import Sulfides
from modules import fluids
from modules.fluids import Water
from modules.geophysics import Elasticity as elast
from modules.chemistry import PeriodicSystem, DataProcessing
from modules.minerals import CrystalPhysics
from modules.geophysics import BoreholeGeophysics as bg
from modules.geophysics import WellLog as wg
from modules.geochemistry import MineralChemistry
from modules.pyllosilicates import Pyllosilicates

class Ores:
    #
    def __init__(self, fluid, actualThickness, porosity=None, data_type=False):
        self.fluid = fluid
        self.actualThickness = actualThickness
        self.porosity = porosity
        self.data_type = data_type
        self.data_water = Water.water("")
    #
    def create_kupferschiefer(self, w_Cu=None, amounts=None):
        #
        results = {}
        results["rock"] = "Kupferschiefer"
        #
        self.w_Cu = w_Cu
        self.amounts = amounts
        #
        # Mineralogy
        quartz = Oxides(data_type=True).create_quartz()
        calcite = Carbonates(data_type=True).create_calcite()
        pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()
        chalcopyrite = Sulfides(impurity="pure", data_type=True).create_chalcopyrite()
        galena = Sulfides(impurity="pure", data_type=True).create_galena()
        bornite = Sulfides(impurity="pure", data_type=True).create_bornite()
        sphalerite = Sulfides(impurity="pure", data_type=True).create_sphalerite()
        chalcocite = Sulfides(impurity="pure", data_type=True).create_chalcocite()
        covellite = Sulfides(impurity="pure", data_type=True).create_covellite()
        digenite = Sulfides(impurity="pure", data_type=True).create_digenite()
        illite = Pyllosilicates(impurity="pure", dict=True).create_illite()
        kaolinite = Pyllosilicates(impurity="pure", dict=True).create_kaolinite()
        montmorillonite = Pyllosilicates(impurity="pure", dict=True).create_montmorillonite()
        #
        mineralogy = [illite, kaolinite, montmorillonite, quartz, calcite, chalcopyrite, bornite, chalcocite, covellite,
                      digenite, galena, sphalerite, pyrite]
        #
        water = fluids.Water.water("")
        #
        data = []
        #
        cond = False
        composition = []
        while cond == False:
            if self.w_Cu == None and self.amounts == None:
                w_clay = round(rd.uniform(0.4, 0.6), 4)
                w_ore = round(rd.uniform(0.4, (1-w_clay)), 4)
                w_misc = round(rd.uniform(0.0, (1-w_clay-w_ore)), 4)
                #
                w_ilt = round(w_clay*rd.uniform(0.5, 1), 4)
                w_kln = round(w_clay*rd.uniform(0, (1-w_ilt)), 4)
                w_mnt = round(w_clay-w_ilt-w_kln, 4)
                #
                magicnumber = rd.randint(1, 7)
                if magicnumber == 1:
                    w_cv = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_bn = round(w_ore*rd.uniform(0, (1-w_cv)), 4)
                    w_cc = round(w_ore*rd.uniform(0, (1-w_cv-w_bn)), 4)
                    w_ccp = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc)), 4)
                    w_dg = round(w_ore-w_cv-w_bn-w_cc-w_ccp, 4)
                    #
                    w_gn = 0.0
                    w_py = 0.0
                    w_sp = 0.0
                elif magicnumber == 2:
                    w_bn = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_ccp = round(w_ore*rd.uniform(0, (1-w_bn)), 4)
                    w_sp = round(w_ore*rd.uniform(0, (1-w_bn-w_ccp)), 4)
                    w_dg = round(w_ore-w_bn-w_ccp-w_sp, 4)
                    #
                    w_cv = 0.0
                    w_cc = 0.0
                    w_gn = 0.0
                    w_py = 0.0
                    w_dg = 0.0
                elif magicnumber == 3:
                    w_cc = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_gn = round(w_ore*rd.uniform(0, (1-w_cc)), 4)
                    w_sp = round(w_ore-w_cc-w_gn, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_ccp = 0.0
                    w_py = 0.0
                    w_dg = 0.0
                elif magicnumber == 4:
                    w_gn = round(w_ore*rd.uniform(0.33, 0.67), 4)
                    w_ccp = round(w_ore*rd.uniform(0.05, (1-w_gn)), 4)
                    w_py = round(w_ore*rd.uniform(0, (1-w_gn-w_ccp)), 4)
                    w_sp = round(w_ore-w_gn-w_ccp-w_py, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_cc = 0.0
                    w_dg = 0.0
                elif magicnumber == 5:
                    w_py = round(w_ore*rd.uniform(0.33, 0.67), 4)
                    w_gn = round(w_ore*rd.uniform(0, (1-w_py)), 4)
                    w_sp = round(w_ore*rd.uniform(0, (1-w_py-w_gn)), 4)
                    w_ccp = round(w_ore-w_py-w_gn-w_sp, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_cc = 0.0
                    w_dg = 0.0
                elif magicnumber == 6:
                    w_sp = round(w_ore*rd.uniform(0.33, 0.67), 4)
                    w_ccp = round(w_ore*rd.uniform(0.05, (1-w_sp)), 4)
                    w_gn = round(w_ore*rd.uniform(0, (1-w_sp-w_ccp)), 4)
                    w_py = round(w_ore-w_sp-w_ccp-w_gn, 4)
                    #
                    w_cv = 0.0
                    w_bn = 0.0
                    w_cc = 0.0
                    w_dg = 0.0
                elif magicnumber == 7:
                    w_dg = round(w_ore*rd.uniform(0.5, 1), 4)
                    w_cv = round(w_ore*rd.uniform(0, (1-w_dg)), 4)
                    w_bn = round(w_ore*rd.uniform(0, (1-w_dg-w_cv)), 4)
                    w_cc = round(w_ore*rd.uniform(0, (1-w_dg-w_cv-w_bn)), 4)
                    w_ccp = round(w_ore-w_dg-w_cv-w_bn-w_cc, 4)
                    #
                    w_gn = 0.0
                    w_py = 0.0
                    w_sp = 0.0
                #
                w_qz = round(w_misc*rd.uniform(0, 1), 4)
                w_cal = round(1 - w_clay - w_ore - w_qz, 4)
            elif self.w_Cu != None:
                w_clay = round(rd.uniform(0.33, 0.67), 4)
                w_ore = round(rd.uniform(0.33, (1-w_clay)), 4)
                w_misc = round(rd.uniform(0.0, (1-w_clay-w_ore)), 4)
                #
                w_ilt = round(w_clay*rd.uniform(0, 1), 4)
                w_kln = round(w_clay*rd.uniform(0, (1-w_ilt)), 4)
                w_mnt = round(w_clay-w_ilt-w_kln, 4)
                #
                w_cv = round(w_ore*rd.uniform(0, 1), 4)
                w_bn = round(w_ore*rd.uniform(0, (1-w_cv)), 4)
                w_cc = round(w_ore*rd.uniform(0, (1-w_cv-w_bn)), 4)
                w_ccp = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc)), 4)
                w_gn = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc-w_ccp)), 4)
                w_py = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc-w_ccp-w_gn)), 4)
                w_dg = round(w_ore*rd.uniform(0, (1-w_cv-w_bn-w_cc-w_ccp-w_gn-w_py)), 4)
                w_sp = round(w_ore-w_cv-w_bn-w_cc-w_ccp-w_gn-w_py-w_dg, 4)
                #
                w_qz = round(w_misc*rd.uniform(0, 1), 4)
                w_cal = 1 - w_clay - w_ore - w_qz
            elif type(self.amounts) is list:
                w_hl = round(abs(np.random.normal(self.amounts[0], 0.025)), 4)
                w_anh = round(abs(np.random.normal(self.amounts[1], 0.025)), 4)
                w_gp = round(abs(np.random.normal(self.amounts[2], 0.025)), 4)
                w_syl = round(1-w_hl-w_anh-w_gp, 4)
            #
            if w_ilt >= 0.0 and w_kln >= 0.0 and w_mnt >= 0.0 and w_cv >= 0.0 and w_bn >= 0.0 and w_cc >= 0.0 \
                    and w_ccp >= 0.0 and w_gn >= 0.0 and w_py >= 0.0 and w_dg >= 0.0 and w_sp >= 0.0 and w_qz >= 0.0 \
                    and w_cal >= 0.0:
                sumMin = round(w_ilt + w_kln + w_mnt + w_cv + w_bn + w_cc + w_ccp + w_gn + w_py + w_dg + w_sp + w_qz + w_cal, 4)
            else:
                sumMin = 0
            #
            w_H = round(w_ilt*illite["chemistry"]["H"] + w_kln*kaolinite["chemistry"]["H"] + w_mnt*montmorillonite["chemistry"]["H"], 4)
            w_C = round(w_cal*calcite["chemistry"]["C"], 4)
            w_Na = round(w_mnt*montmorillonite["chemistry"]["Na"], 4)
            w_Mg = round(w_ilt*illite["chemistry"]["Mg"] + w_mnt*montmorillonite["chemistry"]["Mg"], 4)
            w_Al = round(w_ilt*illite["chemistry"]["Al"] + w_kln*kaolinite["chemistry"]["Al"] + w_mnt*montmorillonite["chemistry"]["Al"], 4)
            w_Si = round(w_ilt*illite["chemistry"]["Si"] + w_kln*kaolinite["chemistry"]["Si"] + w_mnt*montmorillonite["chemistry"]["Si"] + w_qz*quartz["chemistry"]["Si"], 4)
            w_S = round(w_cv*covellite["chemistry"]["S"] + w_bn*bornite["chemistry"]["S"] + w_cc*chalcocite["chemistry"]["S"] + w_ccp*chalcopyrite["chemistry"]["S"] + w_gn*galena["chemistry"]["S"] + w_py*pyrite["chemistry"]["S"] + w_sp*sphalerite["chemistry"]["S"] + w_dg*digenite["chemistry"]["S"], 4)
            w_K = round(w_ilt*illite["chemistry"]["K"], 4)
            w_Ca = round(w_mnt*montmorillonite["chemistry"]["Ca"] + w_cal*calcite["chemistry"]["Ca"], 4)
            w_Fe = round(w_ilt*illite["chemistry"]["Fe"] + w_bn*bornite["chemistry"]["Fe"] + w_ccp*chalcopyrite["chemistry"]["Fe"] + w_py*pyrite["chemistry"]["Fe"], 4)
            w_Cu = round(w_cv*covellite["chemistry"]["Cu"] + w_bn*bornite["chemistry"]["Cu"] + w_cc*chalcocite["chemistry"]["Cu"] + w_ccp*chalcopyrite["chemistry"]["Cu"] + w_dg*digenite["chemistry"]["Cu"], 4)
            w_Zn = round(w_sp*sphalerite["chemistry"]["Zn"], 4)
            w_Pb = round(w_gn*galena["chemistry"]["Pb"], 4)
            w_O = round(1 - w_H - w_C - w_Na - w_Mg - w_Al - w_Si - w_S - w_K - w_Ca - w_Fe - w_Cu - w_Zn - w_Pb, 4)
            sumConc = round(w_H + w_C + w_O + w_Na + w_Mg + w_Al + w_Si + w_S + w_K + w_Ca + w_Fe + w_Cu + w_Zn + w_Pb, 4)
            #print("Amount:", sumMin, "C:", sumConc)
            #
            if sumMin == 1 and sumConc == 1:
                cond = True
                composition.extend((["Ilt", "Kln", "Mnt", "Qz", "Cal", "Ccp", "Bn", "Cc", "Cv", "Dg", "Gn", "Sp", "Py"]))
                concentrations = [w_H, w_C, w_O, w_Na, w_Mg, w_Al, w_Si, w_S, w_K, w_Ca, w_Fe, w_Cu, w_Zn, w_Pb]
                amounts = [w_ilt, w_kln, w_mnt, w_qz, w_cal, w_ccp, w_bn, w_cc, w_cv, w_dg, w_gn, w_sp, w_py]
            else:
                cond = False
        #
        element_list = ["H", "C", "O", "Na", "Mg", "Al", "Si", "S", "K", "Ca", "Fe", "Cu", "Zn", "Pb"]
        mineral_list = ["Ilt", "Kln", "Mnt", "Qz", "Cal", "Ccp", "Bn", "Cc", "Cv", "Dg", "Gn", "Sp", "Py"]
        data.append(composition)
        results["chemistry"] = {}
        results["mineralogy"] = {}
        for index, element in enumerate(element_list, start=0):
            results["chemistry"][element] = concentrations[index]
        for index, mineral in enumerate(mineral_list, start=0):
            results["mineralogy"][mineral] = amounts[index]
        #
        rhoSolid = (w_ilt*illite["rho"] + w_kln*kaolinite["rho"] + w_mnt*montmorillonite["rho"] + w_qz*quartz["rho"]
                    + w_cal*calcite["rho"] + w_ccp*chalcopyrite["rho"] + w_bn*bornite["rho"] + w_cc*chalcocite["rho"]
                    + w_cv*covellite["rho"] + w_dg*digenite["rho"] + w_gn*galena["rho"] + w_sp*sphalerite["rho"]
                    + w_py*pyrite["rho"]) / 1000
        X = [w_ilt, w_kln, w_mnt, w_qz, w_cal, w_ccp, w_bn, w_cc, w_cv, w_dg, w_gn, w_sp, w_py]
        K_list = [mineralogy[i]["K"] for i in range(len(mineralogy))]
        G_list = [mineralogy[i]["G"] for i in range(len(mineralogy))]
        K_geo = elast.calc_geometric_mean(self, X, K_list)
        G_geo = elast.calc_geometric_mean(self, X, G_list)
        K_solid = K_geo
        G_solid = G_geo
        vP_solid = np.sqrt((K_solid*10**9+4/3*G_solid*10**9)/(rhoSolid*10**3))
        vS_solid = np.sqrt((G_solid*10**9)/(rhoSolid*10**3))
        #
        if self.porosity == None:
            if self.actualThickness <= 1000:
                phi = rd.uniform(0.0, 0.1)
            elif self.actualThickness > 1000 and self.actualThickness <= 2000:
                phi = rd.uniform(0.0, 0.075)
            elif self.actualThickness > 2000 and self.actualThickness <= 3000:
                phi = rd.uniform(0.0, 0.05)
            elif self.actualThickness > 3000 and self.actualThickness <= 4000:
                phi = rd.uniform(0.0, 0.025)
            elif self.actualThickness > 4000:
                phi = rd.uniform(0.0, 0.0125)
        else:
            phi = self.porosity
        #
        results["phi"] = phi
        results["fluid"] = self.fluid
        #
        rho = (1 - phi) * rhoSolid + phi * water[2] / 1000
        vP = (1-phi)*vP_solid + phi*water[4][0]
        vS = (1 - phi) * vS_solid
        G_bulk = vS**2 * rho
        K_bulk = vP**2 * rho - 4/3*G_bulk
        E_bulk = (9*K_bulk*G_bulk)/(3*K_bulk+G_bulk)
        phiD = (rhoSolid - rho) / (rhoSolid - water[2] / 1000)
        phiN = (2 * phi ** 2 - phiD ** 2) ** (0.5)
        GR = w_ilt*illite["GR"] + w_kln*kaolinite["GR"] + w_mnt*montmorillonite["GR"] + w_qz*quartz["GR"] \
             + w_cal*calcite["GR"] + w_ccp*chalcopyrite["GR"] + w_bn*bornite["GR"] + w_cc*chalcocite["GR"] \
             + w_cv*covellite["GR"] + w_dg*digenite["GR"] + w_gn*galena["GR"] + w_sp*sphalerite["GR"] + w_py*pyrite["GR"]
        PE = w_ilt*illite["PE"] + w_kln*kaolinite["PE"] + w_mnt*montmorillonite["PE"] + w_qz*quartz["PE"] \
             + w_cal*calcite["PE"] + w_ccp*chalcopyrite["PE"] + w_bn*bornite["PE"] + w_cc*chalcocite["PE"] \
             + w_cv*covellite["PE"] + w_dg*digenite["PE"] + w_gn*galena["PE"] + w_sp*sphalerite["PE"] + w_py*pyrite["PE"]
        poisson_seismic = 0.5*(vP**2 - 2*vS**2)/(vP**2 - vS**2)
        poisson_elastic = (3*K_bulk - 2*G_bulk)/(6*K_bulk + 2*G_bulk)
        poisson_mineralogical = w_ilt*illite["nu"] + w_kln*kaolinite["nu"] + w_mnt*montmorillonite["nu"] + w_qz*quartz["nu"] \
                                + w_cal*calcite["nu"] + w_ccp*chalcopyrite["nu"] + w_bn*bornite["nu"] + w_cc*chalcocite["nu"] \
                                + w_cv*covellite["nu"] + w_dg*digenite["nu"] + w_gn*galena["nu"] + w_sp*sphalerite["nu"] + w_py*pyrite["nu"]
        #
        if self.data_type == False:
            #
            data.append([round(rho, 3), round(rhoSolid, 3), round(water[2] / 1000, 6)])
            data.append([round(K_bulk*10**(-6), 2), round(G_bulk*10**(-6), 2), round(E_bulk*10**(-6), 2), round(poisson_mineralogical, 3)])
            data.append([round(vP, 2), round(vS, 2), round(vP_solid, 2), round(water[4][0], 2)])
            data.append([round(phi, 3), round(phiD, 3), round(phiN, 3)])
            data.append("water")
            data.append([round(GR, 3), round(PE, 3)])
            data.append(concentrations)
            data.append(amounts)
            #
            return data
        else:
            #
            results["rho"] = round(rho*1000, 4)
            results["vP"] = round(vP, 4)
            results["vS"] = round(vS, 4)
            results["vP/vS"] = round(vP/vS, 4)
            results["G"] = round(G_bulk*10**(-6), 4)
            results["K"] = round(K_bulk*10**(-6), 4)
            results["E"] = round(E_bulk*10**(-6), 4)
            results["nu"] = round(poisson_mineralogical, 4)
            results["GR"] = round(GR, 4)
            results["PE"] = round(PE, 4)
            #
            return results
    #
    def create_compact_hematite_ore(self, number, porosity=None):
        #
        data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()           # fixed
        data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()       # fixed
        data_magnetite = Oxides(impurity="pure", data_type=True).create_magnetite()     # fixed
        data_pyrite = Sulfides(impurity="pure", data_type=True).create_pyrite()         # fixed
        data_rutile = Oxides(impurity="pure", data_type=True).create_rutile()           # fixed
        #
        assemblage_minerals = [data_quartz, data_hematite, data_magnetite, data_pyrite, data_rutile]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage_minerals:
            amounts_mineralogy[mineral["mineral"]] = []
            mineral_list.append(mineral["mineral"])
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        amounts_helper = []
        while n < number:
            w_total = 0
            n_minerals = 0
            for mineral in mineral_list:
                if mineral == "Qz":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.10, 0.35), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.35:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.45, 0.80), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.45 <= value <= 0.80:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Mag":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.10, 0.20), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.20:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Py":
                    if n_minerals < len(mineral_list)-1:
                        value = round(rd.uniform(0.0, 0.10), 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.10:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Rt":
                    if n_minerals < len(mineral_list)-1:
                        value = round(1-w_total, 4)
                    else:
                        value = round(1-w_total, 4)
                    if value >= 0.0 and 0.0 <= value <= 0.05:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
            #
            if np.sum(amounts_helper) == 1.0 and n_minerals == len(mineral_list):
                for index, mineral in enumerate(mineral_list):
                    amounts_mineralogy[mineral].append(amounts_helper[index])
                n += 1
                amounts_helper.clear()
            else:
                n += 0
                amounts_helper.clear()
        #
        n = 0
        amounts_helper = {}
        while n < number:
            amounts_helper.clear()
            w_total = 0
            n_elements = 0
            rho_s_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.05), 4)
            else:
                phi_helper = porosity
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in data_quartz["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Qz"][n]*data_quartz["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_hematite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Hem"][n]*data_hematite["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_magnetite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Mag"][n]*data_magnetite["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_pyrite["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Py"][n]*data_pyrite["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_rutile["chemistry"]:
                    if n_elements < len(elements)-1:
                        value = amounts_mineralogy["Rt"][n]*data_rutile["chemistry"][element]
                    else:
                        value = 1-w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                #
                n_elements += 1
            #
            value_total = 0
            for element, value in amounts_helper.items():
                amounts_helper[element] = round(value, 4)
                value_total += int(round(amounts_helper[element]*10000, 4))
            value_total = int(value_total)
            if value_total == 10000:
                condition_sum = True
            else:
                condition_sum = False
            if sum(amounts_helper.values()) == 1.0 or condition_sum == True or w_total == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_quartz["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_quartz["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_hematite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_hematite["PE"], 3)
                    elif mineral == "Mag":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_magnetite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_magnetite["PE"], 3)
                    elif mineral == "Py":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_pyrite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_pyrite["PE"], 3)
                    elif mineral == "Rt":
                        rho_s_helper += round(amounts_mineralogy[mineral][n]*data_rutile["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n]*data_rutile["K"], 3)
                        shearmod_helper += round(0.67*amounts_mineralogy[mineral][n]*data_rutile["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n]*data_rutile["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n]*data_rutile["PE"], 3)
                #
                rho_helper = round((1-phi_helper)*rho_s_helper + phi_helper*self.data_water[2]/1000, 3)
                youngsmod_helper = round((9*bulkmod_helper*shearmod_helper)/(3*bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round((3*bulkmod_helper - 2*shearmod_helper)/(6*bulkmod_helper + 2*shearmod_helper), 3)
                vP_helper = round(((bulkmod_helper*10**9 + 4/3*shearmod_helper*10**9)/(rho_helper))**0.5, 3)
                vS_helper = round(((shearmod_helper*10**9)/(rho_helper))**0.5, 3)
                vPvS_helper_helper = round(vP_helper/vS_helper, 3)
                #
                bulk_properties["rho_s"].append(round(rho_s_helper, 3))
                bulk_properties["rho"].append(rho_helper)
                bulk_properties["K"].append(round(bulkmod_helper, 3))
                bulk_properties["G"].append(round(shearmod_helper, 3))
                bulk_properties["E"].append(youngsmod_helper)
                bulk_properties["nu"].append(poisson_helper)
                bulk_properties["vP"].append(vP_helper)
                bulk_properties["vS"].append(vS_helper)
                bulk_properties["vPvS"].append(vPvS_helper_helper)
                bulk_properties["GR"].append(round(gr_helper, 3))
                bulk_properties["PE"].append(round(pe_helper, 3))
                bulk_properties["phi"].append(round(phi_helper, 3))
                n += 1
            #else:
            #    break
        #
        results = {}
        results["rock"] = "Compact Hematite Ore"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = bulk_properties["phi"][0]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"][0]
            results["rho"] = bulk_properties["rho"][0]
            results["vP"] = bulk_properties["vP"][0]
            results["vS"] = bulk_properties["vS"][0]
            results["vP/vS"] = bulk_properties["vPvS"][0]
            results["K"] = bulk_properties["K"][0]
            results["G"] = bulk_properties["G"][0]
            results["E"] = bulk_properties["E"][0]
            results["nu"] = bulk_properties["nu"][0]
            results["GR"] = bulk_properties["GR"][0]
            results["PE"] = bulk_properties["PE"][0]
        #
        return results

    #
    def create_banded_iron_formation(self, number, porosity=None):
        #
        data_quartz = Oxides(impurity="pure", data_type=True).create_quartz()  # fixed
        data_hematite = Oxides(impurity="pure", data_type=True).create_hematite()  # fixed
        data_magnetite = Oxides(impurity="pure", data_type=True).create_magnetite()  # fixed
        data_kaolinite = Phyllosilicates(impurity="pure", data_type=True).create_kaolinite()  # fixed
        data_goethite = Oxides(impurity="pure", data_type=True).create_goethite()  # fixed
        #
        assemblage_minerals = [data_quartz, data_hematite, data_magnetite, data_kaolinite, data_goethite]
        #
        amounts_mineralogy = {}
        amounts_chemistry = {}
        bulk_properties = {}
        properties = ["rho_s", "rho", "K", "G", "E", "nu", "vP", "vS", "vPvS", "GR", "PE", "phi"]
        for property in properties:
            bulk_properties[property] = []
        mineral_list = []
        elements = []
        for mineral in assemblage_minerals:
            amounts_mineralogy[mineral["mineral"]] = []
            mineral_list.append(mineral["mineral"])
            elements_mineral = list(mineral["chemistry"].keys())
            for element in elements_mineral:
                if element not in elements:
                    elements.append(element)
                    amounts_chemistry[element] = []
        mineral_list.sort()
        elements.sort()
        #
        n = 0
        amounts_helper = []
        while n < number:
            w_total = 0
            n_minerals = 0
            for mineral in mineral_list:
                if mineral == "Qz":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(1 - w_total, 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.35:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Hem":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.25, 0.70), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.25 <= value <= 0.70:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Mag":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.10, 0.25), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.15:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Kln":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.10, 0.25), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.10 <= value <= 0.25:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
                elif mineral == "Goe":
                    if n_minerals < len(mineral_list) - 1:
                        value = round(rd.uniform(0.05, 0.25), 4)
                    else:
                        value = round(1 - w_total, 4)
                    if value >= 0.0 and 0.05 <= value <= 0.25:
                        amounts_helper.append(value)
                        w_total += value
                        n_minerals += 1
            #
            if np.sum(amounts_helper) == 1.0 and n_minerals == len(mineral_list):
                for index, mineral in enumerate(mineral_list):
                    amounts_mineralogy[mineral].append(amounts_helper[index])
                n += 1
                amounts_helper.clear()
            else:
                n += 0
                amounts_helper.clear()
        #
        n = 0
        amounts_helper = {}
        while n < number:
            amounts_helper.clear()
            w_total = 0
            n_elements = 0
            rho_s_helper = 0
            bulkmod_helper = 0
            shearmod_helper = 0
            gr_helper = 0
            pe_helper = 0
            if porosity == None:
                phi_helper = round(rd.uniform(0.0, 0.05), 4)
            else:
                phi_helper = porosity
            #
            for element in elements:
                amounts_helper[element] = 0
                if element in data_quartz["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Qz"][n] * data_quartz["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_hematite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Hem"][n] * data_hematite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_magnetite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Mag"][n] * data_magnetite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_kaolinite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Kln"][n] * data_kaolinite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                if element in data_goethite["chemistry"]:
                    if n_elements < len(elements) - 1:
                        value = amounts_mineralogy["Goe"][n] * data_goethite["chemistry"][element]
                    else:
                        value = 1 - w_total
                    amounts_helper[element] += round(value, 4)
                    w_total += round(value, 4)
                #
                n_elements += 1
            #
            value_total = 0
            for element, value in amounts_helper.items():
                amounts_helper[element] = round(value, 4)
                value_total += int(round(amounts_helper[element] * 10000, 4))
            value_total = int(value_total)
            if value_total == 10000:
                condition_sum = True
            else:
                condition_sum = False
            if sum(amounts_helper.values()) == 1.0 or condition_sum == True or w_total == 1.0:
                for key, value in amounts_helper.items():
                    amounts_chemistry[key].append(round(value, 4))
                for mineral in mineral_list:
                    if mineral == "Qz":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_quartz["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_quartz["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_quartz["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_quartz["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_quartz["PE"], 3)
                    elif mineral == "Hem":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_hematite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_hematite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_hematite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_hematite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_hematite["PE"], 3)
                    elif mineral == "Mag":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_magnetite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_magnetite["PE"], 3)
                    elif mineral == "Kln":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_kaolinite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_kaolinite["PE"], 3)
                    elif mineral == "Goe":
                        rho_s_helper += round(amounts_mineralogy[mineral][n] * data_goethite["rho"], 3)
                        bulkmod_helper += round(amounts_mineralogy[mineral][n] * data_goethite["K"], 3)
                        shearmod_helper += round(0.67 * amounts_mineralogy[mineral][n] * data_goethite["G"], 3)
                        gr_helper += round(amounts_mineralogy[mineral][n] * data_goethite["GR"], 3)
                        pe_helper += round(amounts_mineralogy[mineral][n] * data_goethite["PE"], 3)
                #
                rho_helper = round((1 - phi_helper) * rho_s_helper + phi_helper * self.data_water[2] / 1000, 3)
                youngsmod_helper = round(
                    (9 * bulkmod_helper * shearmod_helper) / (3 * bulkmod_helper + shearmod_helper), 3)
                poisson_helper = round(
                    (3 * bulkmod_helper - 2 * shearmod_helper) / (6 * bulkmod_helper + 2 * shearmod_helper), 3)
                vP_helper = round(
                    ((bulkmod_helper * 10 ** 9 + 4 / 3 * shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vS_helper = round(((shearmod_helper * 10 ** 9) / (rho_helper)) ** 0.5, 3)
                vPvS_helper_helper = round(vP_helper / vS_helper, 3)
                #
                bulk_properties["rho_s"].append(round(rho_s_helper, 3))
                bulk_properties["rho"].append(rho_helper)
                bulk_properties["K"].append(round(bulkmod_helper, 3))
                bulk_properties["G"].append(round(shearmod_helper, 3))
                bulk_properties["E"].append(youngsmod_helper)
                bulk_properties["nu"].append(poisson_helper)
                bulk_properties["vP"].append(vP_helper)
                bulk_properties["vS"].append(vS_helper)
                bulk_properties["vPvS"].append(vPvS_helper_helper)
                bulk_properties["GR"].append(round(gr_helper, 3))
                bulk_properties["PE"].append(round(pe_helper, 3))
                bulk_properties["phi"].append(round(phi_helper, 3))
                n += 1
            #else:
            #    break
        #
        results = {}
        results["rock"] = "Banded Iron Formation"
        if number > 1:
            results["mineralogy"] = amounts_mineralogy
            results["chemistry"] = amounts_chemistry
            results["phi"] = bulk_properties["phi"]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"]
            results["rho"] = bulk_properties["rho"]
            results["vP"] = bulk_properties["vP"]
            results["vS"] = bulk_properties["vS"]
            results["vP/vS"] = bulk_properties["vPvS"]
            results["K"] = bulk_properties["K"]
            results["G"] = bulk_properties["G"]
            results["E"] = bulk_properties["E"]
            results["nu"] = bulk_properties["nu"]
            results["GR"] = bulk_properties["GR"]
            results["PE"] = bulk_properties["PE"]
        else:
            single_amounts_mineralogy = {}
            single_amounts_chemistry = {}
            for mineral, value in amounts_mineralogy.items():
                single_amounts_mineralogy[mineral] = value[0]
            for element, value in amounts_chemistry.items():
                single_amounts_chemistry[element] = value[0]
            results["mineralogy"] = single_amounts_mineralogy
            results["chemistry"] = single_amounts_chemistry
            results["phi"] = bulk_properties["phi"][0]
            results["fluid"] = "water"
            results["rho_s"] = bulk_properties["rho_s"][0]
            results["rho"] = bulk_properties["rho"][0]
            results["vP"] = bulk_properties["vP"][0]
            results["vS"] = bulk_properties["vS"][0]
            results["vP/vS"] = bulk_properties["vPvS"][0]
            results["K"] = bulk_properties["K"][0]
            results["G"] = bulk_properties["G"][0]
            results["E"] = bulk_properties["E"][0]
            results["nu"] = bulk_properties["nu"][0]
            results["GR"] = bulk_properties["GR"][0]
            results["PE"] = bulk_properties["PE"][0]
        #
        return results
    #
## TEST
# print("Test: Compact Hematite Ore")
# print(Ores(fluid="water", actualThickness=0).create_compact_hematite_ore(number=1))
# print(Ores(fluid="water", actualThickness=0).create_compact_hematite_ore(number=100))