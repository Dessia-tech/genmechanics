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
F = 900
B = 0.020
d = 0.020

ground = unidimensional.Body(3*B/2., B/2, name='Ground')
ir1 = unidimensional.Body(0, 0,  name='Inner ring1')
or1 = unidimensional.Body(0, d, name='Outer ring1')
ir2 = unidimensional.Body(B, 0, name='Inner ring2')
or2 = unidimensional.Body(B, d, name='Outer ring2')
ir3 = unidimensional.Body(2*B, 0, name='Inner ring3')
or3 = unidimensional.Body(2*B, d, name='Outer ring4')
ir4 = unidimensional.Body(3*B, 0, name='Inner ring4')
or4 = unidimensional.Body(3*B, 3*d/4., name='Outer ring3')

bodies = [ground, ir1, or1, ir2, or2, ir3, or3, ir4, or4]

s1 = unidimensional.Spring(ir1, or1, k1, 0)
s2 = unidimensional.Spring(ir2, or2, k2, 0)

linear_linkages = [s1, s2]

c1 = unidimensional.UnilateralContact(ir1, ir2, B, name='Inner rings contact 1-2')
c2 = unidimensional.UnilateralContact(or1, or2, B, name='Outer rings contact 1-2')
c3 = unidimensional.UnilateralContact(ground, ir1, -3*B/2., name='Left inner ring contact')
c4 = unidimensional.UnilateralContact(ground, or1, -3*B/2., name='Left outer ring contact')


c5 = unidimensional.UnilateralContact(or2, or3, B, name='Outer rings contact 2-3')
c6 = unidimensional.UnilateralContact(ir3, ir4, B, name='Inner rings contact 3-4')

c7 = unidimensional.UnilateralContact(or4, ground, -3*B/2., name='Right inner ring contact')

s3 = unidimensional.CompressionSpring(or3, ir3, 2*k3, 0, 'bearing 3')
s4 = unidimensional.CompressionSpring(ir4, or4, 2*k4, 0, 'bearing 4')

nonlinear_linkages = [c1, c2, c3, c4, c5, c6, c7, s3, s4]


pl1 = unidimensional.Load(or3, -Fp)
pl2 = unidimensional.Load(ir3, Fp)


id1 = unidimensional.ImposedDisplacement(ground, 0.)
imposed_displacements = [id1]


#nonlinear_imposed_displacements = []
l1 = unidimensional.Load(ir1, F)
loads = [l1, pl1, pl2]

sm = unidimensional.UnidimensionalModel(bodies, linear_linkages, nonlinear_linkages, loads,
                         imposed_displacements)

strain_bearing1 = []
strain_bearing2 = []
strain_bearing3 = []
strain_bearing4 = []
load_values = []
for i in range(nsteps):
    Fi = i*F/float(nsteps+1)
    load_values.append(Fi)
    l1.value = Fi
    result = sm.Solve()
    
    # b1
    u1 = result.displacements[s1.body1]
    u2 = result.displacements[s1.body2]
    print(s1.Strains((u1, u2)))
    
    if u2 - u1 < 0:
        strain_bearing1.append(-k1*(u2-u1))
    else:
        strain_bearing1.append(0.)

    # b2
    u1 = result.displacements[s2.body1]
    u2 = result.displacements[s2.body2]
    if u2 - u1 < 0:
        strain_bearing2.append(-k2*(u2-u1))
    else:
        strain_bearing2.append(0.)

    # b3
    u1 = result.displacements[s3.body1]
    u2 = result.displacements[s3.body2]
    if u2 - u1 < 0:
        strain_bearing3.append(-k3*(u2-u1))
    else:
        strain_bearing3.append(0.)

    # b4
    u1 = result.displacements[s4.body1]
    u2 = result.displacements[s4.body2]
    if u2 - u1 < 0:
        strain_bearing4.append(-k4*(u2-u1))
    else:
        strain_bearing4.append(0.)


result.Plot(intensity_factor=1e-5)
plt.figure()
plt.plot(load_values, strain_bearing1, label='bearing1')
plt.plot(load_values, strain_bearing2, label='bearing2')
plt.plot(load_values, strain_bearing3, label='bearing3')
plt.plot(load_values, strain_bearing4, label='bearing4')
plt.legend()