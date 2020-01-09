#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
# import numpy as npy

import volmdlr as vm
from genmechanics.dynamic_positions import Part, RevoluteLinkage, MovingMechanism, MechanismConfigurations

l = 0.15
#h = 0.32

ground = Part('Ground')
arm = Part('Arm', interest_points=[vm.Point3D((0., 0., 1))])

arm_linkage = RevoluteLinkage(ground, vm.Point3D((-0.07, 0.1, 0.12)), vm.XYZ.Rotation(vm.Z3D, 1.),
                           arm, vm.Point3D((0, 0, 0.)), vm.XYZ,
                           name='arm linkage')

mechanism = MovingMechanism([arm_linkage],
                            ground, name='Crank_rod_mechanism')

manual_configuration = MechanismConfigurations(mechanism, None,
                                               [[0.],
                                                [0.05],
                                                [0.1],
                                                [0.15],
                                                [0.2],
                                                [0.25],
                                                [0.3]])


manual_configuration.babylonjs(plot_frames=True, plot_trajectories=True,
                               plot_instant_rotation_axis=True)

#configurations = mechanism.solve_configurations({0: npy.arange(0, 2*3.14, 0.1)})
#for configuration in configurations:
#    configuration.plot2D(x=vm.y3D, y=vm.z3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.y3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)

#configuration.plot_kinematic_parameters(crank_ground, 0, piston_ground, 0)

#configuration.babylonjs(plot_frames=True)
