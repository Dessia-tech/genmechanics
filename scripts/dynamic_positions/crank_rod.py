#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:54:16 2019

@author: masfaraud
"""
import numpy as npy

import volmdlr as vm
from genmechanics.dynamic_positions import Part, RevoluteLinkage, MovingMechanism, PrismaticLinkage, MechanismConfigurations


l = 0.15
h = 0.32
H = 0.4

ground = Part('Ground')
crank = Part('Crank')
rod = Part('Rod')
piston = Part('piston')

crank_ground = RevoluteLinkage(ground, vm.Point3D((0, 0, 0)), vm.XYZ,
                               crank, vm.Point3D((0, 0, 0)), vm.XYZ,
                               name='crank_ground linkage')

crank_rod = RevoluteLinkage(crank, vm.Point3D((0, l, 0)), vm.XYZ,
                            rod, vm.Point3D((0, 0, 0)), vm.XYZ,
                            name='crank rod linkage')

rod_piston = RevoluteLinkage(rod, vm.Point3D((0, h, 0,)), vm.XYZ,
                             piston, vm.Point3D((-0.1 , 0, 0)), vm.YZX,
                             name='rod_piston linkage')

piston_ground = PrismaticLinkage(ground, vm.Point3D((0, 0, H)), vm.ZXY,
                                 piston, vm.Point3D((0., 0, 0, )), vm.XYZ,
                                 name='piston_ground linkage')

mechanism = MovingMechanism([crank_ground, crank_rod, rod_piston, piston_ground],
                            ground, name='Crank_rod_mechanism')

# manual_configuration = MechanismConfigurations(mechanism, [[0.1, 0, 0]])
# manual_configuration.plot2D(x=vm.y3D, y=vm.z3D)
# manual_configuration.babylonjs()


for initial_configuration in mechanism.find_configurations({0: 0.}, 10, number_starts=10):
    frame_piston = mechanism.part_frame(piston, initial_configuration)
    if frame_piston.origin[2] > 0:
        break

configuration = mechanism.solve_from_initial_configuration(initial_configuration, {0: npy.arange(0.1, 0.5*3.14, 0.1)})
# configuration.plot2D(x=vm.y3D, y=vm.z3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.y3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)

configuration.plot_kinematic_parameters(crank_ground, 0, piston_ground, 0)


configuration.babylonjs(plot_frames=True, plot_instant_rotation_axis=True)
