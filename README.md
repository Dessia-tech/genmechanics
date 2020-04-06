# Genmechanics

A General mechanics solver: Force and speed computing in rigid-body based mechanisms.

## Overview of modules

### Core 
This package beginned with a simple solving of mechanisms defined as parts linked by linkages (ball joints, revolute, prismatic, gear sets...) in a particular configuration. 

Speeds and forces are solved in all the mechanism, enabling to compute power losses and efficiency of mechanisms

![Sankey](https://github.com/Dessia-tech/genmechanics/raw/master/doc/source/images/sankey.png)

### Unidimensional
This module computes non-linear forces and linkage behaviors (unilateral contacts, non linear springs...) for mechanism with parts that have 1D motion.

![unidimensional](https://github.com/Dessia-tech/genmechanics/raw/master/doc/source/images/unidimensional_ballbearings.png)

### Dynamic positions
This is an update from core module where linkages positions are solved to find the mechanism configuration from some imposed linkages parameters.
Mechanism can be rendered with a babylonjs binding.

![crank_rod](https://github.com/Dessia-tech/genmechanics/raw/master/doc/source/images/crank_rod.png)


## Getting started
- install with pip: pip install genmachanics
- execute scripts from scripts folder on github: https://github.com/Dessia-tech/genmechanics/tree/master/scripts


