#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

unidimensional test case

"""

import genmechanics.unidimensional as unidimensional

l = 0.1
k = 100000
F = -500
B = 0.020
d = 0.020

ground = unidimensional.Body(B/2, B/2, name='Ground')
ir1 = unidimensional.Body(0, 0,  name='Inner ring1')
or1 = unidimensional.Body(0, d, name='Outer ring1')
ir2 = unidimensional.Body(B, 0, name='Inner ring2')
or2 = unidimensional.Body(B, d, name='Outer ring2')

bodies = [ground, ir1, ir2, or1, or2]

s1 = unidimensional.Spring(ir1, or1, k, 0)
s2 = unidimensional.Spring(ir2, or2, 2*k, 0)

linear_linkages = [s1, s2]

c1 = unidimensional.UnilateralContact(ir1, ir2, B, name='Inner rings contact')
c2 = unidimensional.UnilateralContact(or1, or2, B, name='Outer rings contact')
c3 = unidimensional.UnilateralContact(ground, ir1, -B/2, name='Left inner ring contact')
c4 = unidimensional.UnilateralContact(ground, or1, -B/2, name='Left outer ring contact')
c5 = unidimensional.UnilateralContact(ir2, ground, -B/2, name='Right inner ring contact')

nonlinear_linkages = [c1, c2, c3, c5]

l1 = unidimensional.Load(or2, F)

loads = [l1]

id1 = unidimensional.ImposedDisplacement(ground, 0.)
imposed_displacements = [id1]


#nonlinear_imposed_displacements = []

sm = unidimensional.UnidimensionalModel(bodies, linear_linkages, nonlinear_linkages, loads,
                         imposed_displacements)
#A = sm.LinearSolve()
result = sm.Solve()

#for result in results:
result.Plot(intensity_factor=1e-5)