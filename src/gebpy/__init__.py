"""
GebPy – Synthetic data generation of minerals and rocks
=========================================================

Public API definition.
Only import stable, user-facing classes and modules here.
"""
# Last updated: 16.12.2025

# Version
__version__ = "1.1.2"

# -------------------------
# High-level Mineral Generators
# -------------------------
from .core.minerals.carbonates import Carbonates
from .core.minerals.oxides import Oxides
from .core.minerals.phyllosilicates import Phyllosilicates
from .core.minerals.sulfides import Sulfides
from .core.minerals.tectosilicates import Tectosilicates

# -------------------------
# High-level Rock Generators
# -------------------------
from .core.rocks.sedimentary import SedimentaryRocks
# from .core.rocks.igneous import IgneousRocks
# from .core.rocks.metamorphic import MetamorphicRocks

# -------------------------
# Public Submodules
# -------------------------
from .core import rocks, minerals, fluids, subsurface

# -------------------------
# Explicit public interface
# -------------------------
__all__ = [
    # Minerals
    "Carbonates",
    "Oxides",
    "Phyllosilicates",
    "Sulfides",
    "Tectosilicates",

    # Rocks
    "SedimentaryRocks",

    # Submodules
    "rocks",
    "minerals",
    "fluids",
    "subsurface",
]
