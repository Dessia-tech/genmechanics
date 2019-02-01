#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Springs test case

"""

import genmechanics.springs as springs

l = 0.1
k = 10000
F = 100

b1 = springs.Body(0)
b2 = springs.Body(l)
b3 = springs.Body(2*l)

bodies = [b1, b2, b3]

s1 = springs.Spring(b1, b2, k, 0)
s2 = springs.Spring(b2, b3, k, 0)

linear_linkages = [s1, s2]
nonlinear_linkages = []

l1 = springs.Load(b3, F)

loads = [l1]


id1 = springs.ImposedDisplacement(b1, 0.)
imposed_displacements = [id1]

nonlinear_imposed_displacements = []

sm = springs.SpringModel(bodies, linear_linkages, nonlinear_linkages, loads,
                    imposed_displacements, nonlinear_imposed_displacements)

result = sm.SolveLinear()

result.Plot()