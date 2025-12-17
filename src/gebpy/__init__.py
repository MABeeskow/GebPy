"""
GebPy – Synthetic data generation of minerals and rocks
=========================================================

Public API definition.
Only import stable, user-facing classes and modules here.
"""
# Last updated: 17.12.2025

# Version
__version__ = "1.1.5"

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
from .core.rocks.isotropic_rocks import IsotropicRocks
from .core.rocks.anisotropic_rocks import AnisotropicRocks
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
    "IsotropicRocks",
    "AnisotropicRocks",

    # Submodules
    "rocks",
    "minerals",
    "fluids",
    "subsurface",
]
