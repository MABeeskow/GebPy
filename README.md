[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

<a href="https://github.com/MABeeskow/GebPy">
  <img src="https://raw.githubusercontent.com/MABeeskow/GebPy/master/documents/readme_images/readme_header.png"
width="67%">
</a>

developed by Maximilian Alexander Beeskow

GebPy is an open-source, Python-based tool for the synthetic generation of geological data with focus on minerals, 
rocks and stratigraphy. The main assumption of GepPy is that all rock properties are determined by the mineral 
assemblage besides structural features. GebPy can be used for educational purposes, for example by the generation of 
mineralogical data that could be then further investigated in specific diagrams or in well log analysis courses after 
creating different stratigraphic sequences.\
GebPy is a one-man project and can still contain some technical bugs or some wrong assumptions, but it is driven by a 
huge motivation and has the goal to improve the knowledge and understanding of different geological aspects that are 
(actually) focused on mineralogy, rocks and stratigraphy.

Last updated: 02.12.2025

## 💭 Citing GebPy

Coming soon ...

## 🚀 Installation

GebPy can be installed by downloading this repository or - and we encourage this possibility - by the command 
<code>pip install gebpy</code>. If installed, it can be updated by the command <code>pip install gebpy --upgrade</code>.

## 💻 Resources & documentation

Coming soon ...

## 💎 Mineral data

Synthetic geophysical and compositional data are created based on the element properties that are determined by the 
mineral formula, and by the cell parameters which define the crystal lattice. In addition, the bulk and shear modulus 
is also taken from the literature or online databases. Based on this set of input parameters, GebPy calculates the 
* molar mass $M$, 
* molar and cell volume $V_m$ and $V$, 
* density $\varrho$,
* seismic velocities $v_P$ and $v_S$,
* elastic properties $E$ and $\nu$,
* electromagnetic properties $GR$, $PE$, $U$.

In addition, GebPy calculates the chemical composition of the mineral.

## 🪨 Rock data

Based on the sets of synthetic mineral data, GebPy calculates the bulk properties of a rock, including its chemical 
composition. The following geophysical properties are calculated by GebPy:

* bulk density $\varrho_B$,
* porosity $\phi$,
* seismic velocities $v_P$ and $v_S$,
* elastic properties $E$ and $\nu$,
* electromagnetic properties $GR$, $PE$, $U$.

The GebPy database contains currently more than 100 minerals that can be used for rock modeling. Many 
common rocks like sandstone, limestone, shale, many igneous and volcanic rocks but also some metamorphic and ore rocks 
are part of the database as well.

## 🏞️ Stratigraphic data

Since GebPy is able to generate synthetic rock data, it is also possible to stack them. Please be aware that GebPy does 
not consider effects that are related to structural geology, tectonics, fluid flow, geothermal heat, lithostatic 
pressure. 

## 📚 References

Physical, chemical and crystallographic parameters are mainly taken from websites like 
[webelements.com](https://webelements.com), [webminerals.com](https://webminerals.com), [mindat.org](https://mindat.org)
and [materialsproject.org](https://materialsproject.org). Physical, chemical and crystallographic equations and 
principles are generally taken from standard textbooks about these subjects.
However, a more detailed list will be added soon.