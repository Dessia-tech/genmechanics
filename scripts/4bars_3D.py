#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4 bars mechanism
"""

import volmdlr as vm
from genmechanics.dynamic_positions import Part, RevoluteLinkage, MovingMechanism

ground = Part('ground')
rod1 = Part('rod 1', interest_points=[vm.o3D])
rod2 = Part('rod 2', interest_points=[vm.o3D])
rod3 = Part('rod 3')

l1 = 0.15
l2 = 0.18
l3 = 0.16
h = 0.17


axis1 = vm.Vector3D((0.9, 0.2, 0))
axis4 = vm.Vector3D((-0.2, 0.4, 0))

axis1.Normalize()
axis4.Normalize()

linkage1 = RevoluteLinkage(ground, vm.Point3D((-0.5*h, -0.3*h, 0)), vm.Basis3D.from_two_vectors(axis1, vm.y3D),
                           rod1,  vm.Point3D((-0.5*l1, 0, 0)), vm.YZX, name='l1')


linkage2 = RevoluteLinkage(rod2, vm.Point3D((-0.5*l2,0, 0)), vm.YZX,
                           rod1, vm.Point3D((0.5*l1, 0, 0)), vm.YZX, name='l2')
#
#linkage3 = RevoluteLinkage(rod3, vm.Point3D((-0.5*l3, 0, 0,)), vm.z3D,
#                           rod2, vm.Point3D((0.5*l2,0, 0)), name='l3')
#
#linkage4 = RevoluteLinkage(ground, vm.Point3D((0.5*h, 0, 0)), axis4,
#                           rod3, vm.Point3D((0.5*l3, 0, 0, )), name='l4')

#mechanism = MovingMechanism([linkage1, linkage2, linkage3, linkage4], ground, name='3 bar')
mechanism = MovingMechanism([linkage1, linkage2], ground, name='3 bar')
#
configuration = mechanism.solve_configurations({0: [0.01, 0.05, 0.15, 0.2, 0.25,
                                                    0.3, 0.35, 0.4, 0.45, 0.5, 0.6,
                                                    0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6,
                                                    1.8, 2, 2.2, 2.4, 2.6, 2.8, 3,
                                                    3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4,
                                                    4.6, 4.8, 5, 5.2, 5.4, 5.6],
                                                1: [-0.01, -0.05, -0.15, -0.2, -0.25,
                                                    -0.3, -0.35, -0.4, -0.45, -0.5, -0.6,
                                                    -0.7, -0.8, -0.9, -1, -1.2, -1.4, -1.6,
                                                    -1.8, -2, -2.2, -2.4, -2.6, -2.8, -3,
                                                    -3.2, -3.4, -3.6, -3.8, -4, -4.2, -4.4,
                                                    -4.6, -4.8, -5, -5.2, -5.4, -5.6]
                                                })

configuration.plot2D()

#configuration.plot_kinematic_parameters(linkage1, 0, linkage4, 0)

#manual_configuration = MechanismConfiguration(mechanism, [3.1415/4, 3.1415/8, -3.1415/3])
#manual_configuration.plot2D()

#f = mechanism.part_frame(rod2, [3.1415/8, -3.1415/8, -3.1415/10])

#reduced_mechanism = MovingMechanism([l1, l2, l3], ground, '3 bar')
#manual_configuration_reduced_mechanism = MechanismConfiguration(reduced_mechanism, [0.37, 0.15, -0.04])
#manual_configuration_reduced_mechanism.plot2D()
configuration.babylonjs(plot_frames=True, plot_trajectories=True)