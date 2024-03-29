from setuptools import setup

setup(
	name="modules", 
	packages=["modules"]
)

setup(
	name="GebPy", 
	packages=["documents"],
	package_data={"documents": ["*.png", "dir/*"]}
)

import modules
from modules.gui_elements import SimpleElements as SE
from modules.sulfates import Sulfates
from modules.oxides import Oxides
from modules.sulfides import Sulfides
from modules.carbonates import Carbonates
from modules.halides import Halogenes
from modules.silicates import Tectosilicates, Phyllosilicates, Nesosilicates, Sorosilicates, Inosilicates
from modules.phosphates import Phosphates
from modules.siliciclastics import sandstone, shale
from modules.carbonates import limestone, dolomite
from modules.ore import Ores
from modules.igneous import Plutonic
from modules.evaporites import Evaporites
from modules.sequences import DataProcessing as DP
