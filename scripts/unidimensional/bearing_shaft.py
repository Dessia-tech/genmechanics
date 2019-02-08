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
ir1 = unidimensional.Body(-B, 0,  name='Inner ring1')
or1 = unidimensional.Body(-B, d, name='Outer ring1')
ir2 = unidimensional.Body(B, 0, name='Inner ring2')
or2 = unidimensional.Body(B, d, name='Outer ring2')
ir3 = unidimensional.Body(2*B, 0, name='Inner ring3')
or3 = unidimensional.Body(2*B, d, name='Outer ring4')
ir4 = unidimensional.Body(L, 0, name='Inner ring4')
or4 = unidimensional.Body(L, d, name='Outer ring3')

bodies = [ground, shaft, ir1, or1, ir2, or2, ir3, or3, ir4, or4]

#s1 = unidimensional.Spring(ir1, or1, 1, 0)
s11 = unidimensional.CompressionSpring(ir1, or1, k1, -1e-5, 'bearing 2')
s12 = unidimensional.CompressionSpring(or1, ir1, k1, -1e-5, 'bearing 3')
#s4 = unidimensional.Spring(ir4, or4, 1, 1e-2)
s41 = unidimensional.CompressionSpring(ir4, or4, k1, -1e-4, 'bearing 2')
s42 = unidimensional.CompressionSpring(or4, ir4, k1, -1e-4, 'bearing 3')
#pl41 = unidimensional.Load(or4, 1)
#pl42 = unidimensional.Load(ir4, 1)

linear_linkages = []

c1 = unidimensional.UnilateralContact(ir1, ir2, 2*B, name='Inner rings contact 1-2')
c2 = unidimensional.UnilateralContact(or1, or2, 2*B, name='Outer rings contact 1-2')
c3 = unidimensional.UnilateralContact(ground, or1, - B, name='Left outer ring contact')
c4 = unidimensional.UnilateralContact(shaft, ir1, - B, name='Left inner ring contact')


c5 = unidimensional.UnilateralContact(or2, or3, B, name='Outer rings contact 2-3')
c6 = unidimensional.UnilateralContact(ir2, ir3, B, name='Inner rings contact 3-4')

c7 = unidimensional.UnilateralContact(or3, ground, -2*B, name='Right inner ring contact')
c8 = unidimensional.UnilateralContact(ir3, shaft, -2*B, name='Right inner ring contact')

s2 = unidimensional.CompressionSpring(ir2, or2, k1, -1e-4, 'bearing 2')
s3 = unidimensional.CompressionSpring(or3, ir3, k1, -1e-4, 'bearing 3')

c9 = unidimensional.UnilateralContact(or4, ground, -L, name='Left outer ring contact')
c10 = unidimensional.UnilateralContact(shaft, ir4, L, name='Left inner ring contact')

nonlinear_linkages = [c1, c2, c3, c4, c5, c6, c7, c8, s2, s3, c9, c10, s11, s12, s41, s42]

p1 = unidimensional.Load(ir1, 500)
pl1 = unidimensional.Load(or3, -Fp)
pl2 = unidimensional.Load(ir3, Fp)


id1 = unidimensional.ImposedDisplacement(ground, 0.)
imposed_displacements = [id1]


#nonlinear_imposed_displacements = []
l1 = unidimensional.Load(shaft, F)
loads = [p1, l1, pl1, pl2]

sm = unidimensional.UnidimensionalModel(bodies, linear_linkages, nonlinear_linkages, loads,
                         imposed_displacements)

result = sm.Solve(500)

#for result in results:
result.Plot(intensity_factor=1e-5)