from setuptools import setup, find_packages

with open("README.md", "r", encoding="latin-1") as f:
    description = f.read()

setup(
    name="gebpy",
    version="1.1.5",
    package_dir={"": "src"},
    packages=find_packages(where="src"),

    url="https://github.com/MABeeskow/GebPy",
    license="LGPL-3.0",
    author="Maximilian Alexander Beeskow",
    author_email="maximilian.beeskow@rwth-aachen.com",

    description=(
        "GebPy is a Python-based, open source tool for the generation of geological "
        "data of minerals, rocks and complete lithological sequences."
    ),

    keywords=[
        "geophysics", "geochemistry", "mineralogy",
        "petrology", "seismology", "stratigraphy",
        "geology", "geosciences"
    ],

    install_requires=[
        "numpy", "scipy", "pandas", "matplotlib"
    ],

    entry_points={
        "console_scripts": [
            "gebpy = gebpy.gebpy_app:gebpy"
        ]
    },

    include_package_data=True,
    package_data={"": ["lib/images/*.png"]},

    long_description=description,
    long_description_content_type="text/markdown",
)
