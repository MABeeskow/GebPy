#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		manual_test_elastic_calibration.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		31.01.2026

#-----------------------------------------------

"""
Module: manual_test_elastic_calibration.py
Manual test file related to module common.py
"""

# PACKAGES
import time
import numpy as np

# MODULES
from src.gebpy.core.rocks.common import ElasticCalibration as EC

# CODE
print("\n--- Manual test for: ElasticCalibration in rocks/common.py ---\n")
n_dataset = 5
sample_rho_bulk = np.random.normal(loc=3540, scale=20, size=n_dataset)
sample_vp_bulk = np.random.normal(loc=5830, scale=30, size=n_dataset)
sample_vs_bulk = np.random.normal(loc=2980, scale=30, size=n_dataset)

print("--- Input data:")
print("Results: bulk density rho")
print(min(sample_rho_bulk), max(sample_rho_bulk))
print("Results: p-wave velocity vp")
print(min(sample_vp_bulk), max(sample_vp_bulk))
print("Results: s-wave velocity vs")
print(min(sample_vs_bulk), max(sample_vs_bulk), "\n")

sample_g_bulk = EC.determine_shear_modulus(rho_bulk=sample_rho_bulk, vs_bulk=sample_vs_bulk)
sample_k_bulk = EC.determine_bulk_modulus(rho_bulk=sample_rho_bulk, vp_bulk=sample_vp_bulk, vs_bulk=sample_vs_bulk)
sample_k_bulk *= 1e-9
sample_g_bulk *= 1e-9

print("--- Results: single modulus determination")
print("Results: shear modulus G")
print(min(sample_g_bulk), np.mean(sample_g_bulk), max(sample_g_bulk))
print("Results: bulk modulus K")
print(min(sample_k_bulk), np.mean(sample_k_bulk), max(sample_k_bulk), "\n")

sample_k_bulk, sample_g_bulk = EC.determine_elastic_moduli(
    rho_bulk=sample_rho_bulk, vp_bulk=sample_vp_bulk, vs_bulk=sample_vs_bulk)
sample_k_bulk *= 1e-9
sample_g_bulk *= 1e-9

print("--- Results: combined modulus determination")
print("Results: shear modulus G")
print(min(sample_g_bulk), np.mean(sample_g_bulk), max(sample_g_bulk))
print("Results: bulk modulus K")
print(min(sample_k_bulk), np.mean(sample_k_bulk), max(sample_k_bulk), "\n")

scl_rho_bulk = EC.determine_scaling_factors(reference_data=3540, sample_data=sample_rho_bulk)
scl_rho_vp = EC.determine_scaling_factors(reference_data=5830, sample_data=sample_vp_bulk)
scl_rho_vs = EC.determine_scaling_factors(reference_data=2980, sample_data=sample_vs_bulk)

print("--- Results: scaling factor determination")
print("Results: bulk density rho_b")
print(min(scl_rho_bulk), np.mean(scl_rho_bulk), max(scl_rho_bulk))
print("Results: p-wave velocity vp")
print(min(scl_rho_vp), np.mean(scl_rho_vp), max(scl_rho_vp))
print("Results: p-wave velocity vs")
print(min(scl_rho_vs), np.mean(scl_rho_vs), max(scl_rho_vs), "\n")

ref_rho_bulk = np.random.normal(loc=3540, scale=30, size=n_dataset)
ref_k_bulk = np.random.normal(loc=105, scale=1.5, size=n_dataset)
ref_g_bulk = np.random.normal(loc=25, scale=1.5, size=n_dataset)
ref_vp_bulk = ((ref_k_bulk*1e9 + 4/3*ref_g_bulk*1e9)/(ref_rho_bulk))**0.5
ref_vs_bulk = ((ref_g_bulk*1e9)/(ref_rho_bulk))**0.5

print("--- Reference model data:")
print("Results: bulk modulus K")
print(min(ref_k_bulk), np.mean(ref_k_bulk), max(ref_k_bulk))
print("Results: shear modulus G")
print(min(ref_g_bulk), np.mean(ref_g_bulk), max(ref_g_bulk))
print("Results: p-wave velocity vp")
print(min(ref_vp_bulk), np.mean(ref_vp_bulk), max(ref_vp_bulk))
print("Results: s-wave velocity vs")
print(min(ref_vs_bulk), np.mean(ref_vs_bulk), max(ref_vs_bulk), "\n")

opt_elastic_params = EC.adjust_elastic_model_parameters(
    rho_smpl=sample_rho_bulk, vp_smpl=sample_vp_bulk, vs_smpl=sample_vs_bulk, k_ref=ref_k_bulk, g_ref=ref_g_bulk)
k_opt = opt_elastic_params["K_opt"]
g_opt = opt_elastic_params["G_opt"]
opt_vp_bulk = ((k_opt*1e9 + 4/3*g_opt*1e9)/(ref_rho_bulk))**0.5
opt_vs_bulk = ((g_opt*1e9)/(ref_rho_bulk))**0.5
difference = EC.determine_difference_from_ideality(
    w_k=opt_elastic_params["K_weight"], w_g=opt_elastic_params["G_weight"])

print("--- Results: elastic moduli optimization:")
print("Results: bulk modulus K")
print(min(k_opt), np.mean(k_opt), max(k_opt))
print("Results: shear modulus G")
print(min(g_opt), np.mean(g_opt), max(g_opt))
print("Results: p-wave velocity vp")
print(min(opt_vp_bulk), np.mean(opt_vp_bulk), max(opt_vp_bulk))
print("Results: s-wave velocity vs")
print(min(opt_vs_bulk), np.mean(opt_vs_bulk), max(opt_vs_bulk), "\n")
print("Results: difference from isotropic ideality")
print(difference, "\n")