#!/usr/bin/env python
# -*-coding: utf-8 -*-

# -----------------------------------------------

# Name:		test_rocks.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		21.11.2022

# -----------------------------------------------

## MODULES
import sys, time
import numpy as np
import random as rd
import  matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from modules import sequences, geophysics
from modules import igneous
from modules.siliciclastics import Sandstone
from modules.evaporites import Evaporites
from modules.carbonates import CarbonateRocks

# ## TESTING
# # Test soil generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_soil = data.create_soil(grainsize=True, thickness=10)
# for i in range(len(data_soil)):
#     print(data_soil[i])
#
# print("")
# # Test sand generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_sand = data.create_sand(grainsize=True, thickness=10)
# for i in range(len(data_sand)):
#     print(data_sand[i])
#
# print("")
# # Test sandstone generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_sandstone = data.create_sandstone(thickness=10)
# for i in range(len(data_sandstone)):
#     print(data_sandstone[i])
#
# print("")
# data = sequences.SedimentaryBasin()
# data_sandstone = data.create_sandstone(fluid="oil", thickness=10)
# for i in range(len(data_sandstone)):
#     print(data_sandstone[i])
#
# print("")
# data = sequences.SedimentaryBasin()
# data_sandstone = data.create_sandstone(fluid="gas", thickness=10)
# for i in range(len(data_sandstone)):
#     print(data_sandstone[i])
#
# print("")
# data = sequences.SedimentaryBasin()
# data_sandstone = data.create_sandstone(fluid="air", thickness=10)
# for i in range(len(data_sandstone)):
#     print(data_sandstone[i])
#
# print("")
# data = sequences.SedimentaryBasin()
# data_sandstone = data.create_sandstone(fluid="air", thickness=10, phi=0.15, Qz_rich=True)
# for i in range(len(data_sandstone)):
#     print(data_sandstone[i])
#
# print("")
# print("Feldspathic Sandstone Generation:")
# # Test feldspathic sandstone generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_sandstone = data.create_sandstone(thickness=10, keyword="feldspathic")
# for i in range(len(data_sandstone)):
#     print(data_sandstone[i])
#
# print("")
# # Test limestone generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_limestone = data.create_limestone(thickness=10)
# for i in range(len(data_limestone)):
#     print(data_limestone[i])
#
# print("")
# data = sequences.SedimentaryBasin()
# data_limestone = data.create_limestone(fluid="oil", thickness=10)
# for i in range(len(data_limestone)):
#     print(data_limestone[i])
#
# print("")
# data = sequences.SedimentaryBasin()
# data_limestone = data.create_limestone(fluid="gas", thickness=10)
# for i in range(len(data_limestone)):
#     print(data_limestone[i])
#
# print("")
# # Test shale generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_shale = data.create_shale(thickness=10)
# for i in range(len(data_shale)):
#     print(data_shale[i])
#
# #print("")
# # Test shale generation within SedimentaryBasin class
# #data = sequences.SedimentaryBasin()
# #data_shale = data.create_shale_complex(thickness=10)
# #for i in range(len(data_shale)):
# #    print(data_shale[i])
#
# print("")
# # Test rock salt generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_rocksalt = data.create_rocksalt(thickness=10)
# for i in range(len(data_rocksalt)):
#     print(data_rocksalt[i])
#
# print("")
# print(":: IGNEOUS ROCKS ::")
# # Test granite generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_granite = data.create_granite(thickness=10)
# for i in range(len(data_granite)):
#     print(data_granite[i])
#
# print("")
# # Test basalt generation within SedimentaryBasin class
# data = sequences.SedimentaryBasin()
# data_basalt = data.create_basalt(thickness=10)
# for i in range(len(data_basalt)):
#     print(data_basalt[i])
#
# print("")
# print(":: IGNEOUS ROCKS (PLUTONITE) ::")
# # Test granite generation within Plutonite class
# data = sequences.Plutonite()
# data_rock = data.create_granite(thickness=10)
# for i in range(len(data_rock)):
#     print(data_rock[i])
# print("")
# # Test gabbro generation within Plutonite class
# data = sequences.Plutonite()
# data_rock = data.create_gabbro(thickness=10)
# for i in range(len(data_rock)):
#     print(data_rock[i])
# print("")
# # Test felsic rock generation within Plutonite class
# data = sequences.Plutonite()
# data_rock = data.create_felsic(thickness=10)
# for i in range(len(data_rock)):
#     print(data_rock[i])
# print("")
# # Test intermediate rock generation within Plutonite class
# data = sequences.Plutonite()
# data_rock = data.create_intermediate(thickness=10)
# for i in range(len(data_rock)):
#     print(data_rock[i])

## SEDIMENTARY ROCKS ##
## SILICICLASTICS
# for i in range(10):
#     data_conglomerate = Sandstone().create_conglomerate(number=10)
#     print(i, data_conglomerate)
#     data_conglomerate = Sandstone().create_conglomerate_alt(number=10)
#     print(i, data_conglomerate)
    # data_mudstone = Sandstone().create_mudstone(number=1)
    # print(i, data_mudstone)
## CARBONATES
start = time.time()
data_container = []
for i in range(100):
    data_limestone = CarbonateRocks().create_limestone_alternative(number=1)
    data_container.append(data_limestone)
end1 = time.time()
data_limestone = CarbonateRocks().create_limestone_alternative(number=100)
end2 = time.time()
print(end1 - start)
print(end2 - end1)
#     data_dolomite = CarbonateRocks().create_dolomite(number=1)
#     print(i, data_dolomite)
#     data_limestone = CarbonateRocks().create_limestone(number=1)
#     print(i, data_limestone)

## EVAPORITE ROCKS ##
# for i in range(10):
#     data_anhydrite = Evaporites(fluid="water", actualThickness=0).create_anhydrite_rock(number=1)
#     print(i, data_anhydrite)
#     data_rocksalt = Evaporites(fluid="water", actualThickness=0).create_rocksalt(number=1)
#     print(i, data_rocksalt)
#     data_potash = Evaporites(fluid="water", actualThickness=0).create_potash(number=1)
#     print(i, data_potash)

## IGNEOUS ROCKS ##
## PLUTONICS
# for i in range(10):
    # data_granite = igneous.Plutonic(fluid="water", actualThickness=100).create_granite_streckeisen()
    # print(i, data_granite)
    # data_granodiorite = igneous.Plutonic(fluid="water", actualThickness=100).create_granodiorite_streckeisen()
    # print(i, data_granodiorite)
    # data_tonalite = igneous.Plutonic(fluid="water", actualThickness=100).create_tonalite_streckeisen()
    # print(i, data_tonalite)
    # data_gabbro = igneous.Plutonic(fluid="water", actualThickness=100).create_gabbro_streckeisen()
    # print(i, data_gabbro)
    # data_diorite = igneous.Plutonic(fluid="water", actualThickness=100).create_diorite_streckeisen()
    # print(i, data_diorite)
    # data_monzonite = igneous.Plutonic(fluid="water", actualThickness=100).create_monzonite_streckeisen()
    # print(i, data_monzonite)
    # data_syenite = igneous.Plutonic(fluid="water", actualThickness=100).create_syenite_streckeisen()
    # print(i, data_syenite)
    # data_granitoid = igneous.Plutonic(fluid="water", actualThickness=100).create_granitoid_streckeisen()
    # print(i, data_granitoid)
    # data_quarzolite = igneous.Plutonic(fluid="water", actualThickness=100).create_quarzolite_streckeisen()
    # print(i, data_quarzolite)
#
## VOLCANICS
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_volcanic_rock()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_rhyolite()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_alkaline_rhyolite()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_dacite()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_basalt()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_andesite()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_trachyte()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_alkaline_trachyte()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_quartzitic_trachyte()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_quartzitic_alkaline_trachyte()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_latite()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_quartzitic_latite()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_rhyolite_generalized()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_trachyte_generalized()
#     print(i, data)
# for i in range(10):
#     data = igneous.Volcanic(fluid="water", actualThickness=100).create_latite_generalized()
#     print(i, data)
## PYROCLASTIC ROCKS ##
# for i in range(10):
#     data = igneous.Pyroclastic(fluid="water", actualThickness=100).create_pyroclastic_rock()
#     print(i, data)