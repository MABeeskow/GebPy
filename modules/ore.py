#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		ore.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		20.11.2021

# -----------------------------------------------

# MODULES
import numpy as np
import random as rd
from modules.oxides import Oxides
from modules.carbonates import Carbonates
from modules.sulfides import Sulfides
from modules import fluids
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
        calcite = Carbonates(dict=True).create_calcite()
        pyrite = Sulfides(impurity="pure", dict=True).create_pyrite()
        chalcopyrite = Sulfides(impurity="pure", dict=True).create_chalcopyrite()
        galena = Sulfides(impurity="pure", dict=True).create_galena()
        bornite = Sulfides(impurity="pure", dict=True).create_bornite()
        sphalerite = Sulfides(impurity="pure", dict=True).create_sphalerite()
        chalcocite = Sulfides(impurity="pure", dict=True).create_chalcocite()
        covellite = Sulfides(impurity="pure", dict=True).create_covellite()
        digenite = Sulfides(impurity="pure", dict=True).create_digenite()
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