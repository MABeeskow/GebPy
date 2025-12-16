#!/usr/bin/env python
# -*-coding: utf-8 -*-

#-----------------------------------------------

# Name:		test_imports.py
# Author:	Maximilian A. Beeskow
# Version:	1.0
# Date:		15.12.2025

#-----------------------------------------------

"""
Module: test_imports.py
Test file related to general import processes.
"""

def test_basic_import():
    import gebpy


def test_version_exists():
    import gebpy
    assert hasattr(gebpy, "__version__")
    assert isinstance(gebpy.__version__, str)


def test_public_api_imports():
    from gebpy import (
        Carbonates,
        Phyllosilicates,
        Tectosilicates,
        Oxides,
        Sulfides,
        SedimentaryRocks,
    )


def test_public_api_all():
    import gebpy
    for name in gebpy.__all__:
        assert hasattr(gebpy, name), f"{name} missing from gebpy namespace"


def test_deep_imports():
    from gebpy.core.minerals.carbonates import Carbonates
    from gebpy.core.rocks.isotropic_rocks import SedimentaryRocks


def test_carbonates_smoke():
    from gebpy.core.minerals.carbonates import Carbonates

    mineral = Carbonates(name="Calcite", random_seed=0)
    data = mineral.generate_dataset(number=1)

    assert isinstance(data, dict)
    assert "rho" in data