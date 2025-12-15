"""
GebPy – Synthetic data generation of minerals and rocks
=========================================================

Public API definition.
Only import stable, user-facing classes and modules here.
"""

# Version
__version__ = "0.1.0"

# -------------------------
# High-level Mineral Generators
# -------------------------
from gebpy.core.minerals.carbonates import Carbonates
from gebpy.core.minerals.oxides import Oxides
from gebpy.core.minerals.phyllosilicates import Phyllosilicates
from gebpy.core.minerals.sulfides import Sulfides
from gebpy.core.minerals.tectosilicates import Tectosilicates

# -------------------------
# High-level Rock Generators
# -------------------------
from gebpy.core.rocks.sedimentary import SedimentaryRocks
# from gebpy.core.rocks.igneous import IgneousRocks
# from gebpy.core.rocks.metamorphic import MetamorphicRocks

# -------------------------
# Public Submodules
# -------------------------
from gebpy.core import rocks, minerals, fluids, subsurface

# -------------------------
# Explicit public interface
# -------------------------
__all__ = [
    # Minerals
    "Phyllosilicates",
    "Tectosilicates",
    "Oxides",
    "Sulfides",

    # Rocks
    "SedimentaryRocks",

    # Submodules
    "rocks",
    "minerals",
    "fluids",
    "subsurface",
]
