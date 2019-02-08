#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Springs test case

"""

import genmechanics.unidimensional as unidimensional
import matplotlib.pyplot as plt

l = 0.1
k1 = 100000
k2 = 110000
k3 = 180000
k4 = 150000
Fp = 500 # Preload force
nsteps = 50
F = -900
B = 0.020
d = 0.020
L = 0.3



shaft = unidimensional.Body(0, -d, name='Shaft')
ground = unidimensional.Body(0, d/2., name='Ground')
ir1 = unidimensional.Body(0, 0,  name='Inner ring1')
or1 = unidimensional.Body(0, d, name='Outer ring1')

bodies = [ground, shaft, ir1, or1]

nonlinear_linkages = []
nonlinear_linkages.append(unidimensional.CompressionSpring(ir1, or1, k1, -1e-5, 'bearing 2'))
nonlinear_linkages.append(unidimensional.CompressionSpring(or1, ir1, k1, -1e-5, 'bearing 3'))

linear_linkages = []

nonlinear_linkages.append(unidimensional.UnilateralContact(shaft, ir1, 0, name='Inner rings contact 1-2'))
nonlinear_linkages.append(unidimensional.UnilateralContact(or1, ground, 0, name='Outer rings contact 1-2'))

p1 = unidimensional.Load(shaft, 500)

id1 = unidimensional.ImposedDisplacement(ground, 0.)
imposed_displacements = [id1]

loads = [p1]

sm = unidimensional.UnidimensionalModel(bodies, linear_linkages, nonlinear_linkages, loads,
                         imposed_displacements)

result = sm.Solve(500)

#for result in results:
result.Plot(intensity_factor=1e-5)