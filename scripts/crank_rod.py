#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:54:16 2019

@author: masfaraud
"""
import numpy as npy

import volmdlr as vm
from genmechanics.dynamic_positions import Part, RevoluteLinkage, MovingMechanism, SlidingRevoluteLinkage, PrismaticLinkage, MechanismConfigurations


l = 0.15
h = 0.32

ground = Part('Ground')
crank = Part('Crank')
rod = Part('Rod')
piston = Part('piston', interest_points=[vm.Point3D((0., 0., 0.1))])

crank_ground = RevoluteLinkage(ground, vm.Point3D((0, 0, 0)), vm.x3D,
                               crank, vm.Point3D((0, 0, 0)),
                               name='crank_ground linkage')

crank_rod = RevoluteLinkage(crank, vm.Point3D((0, l, 0)), vm.x3D,
                            rod, vm.Point3D((0, 0, 0)),
                            name='crank rod linkage')

rod_piston = RevoluteLinkage(rod, vm.Point3D((0., h, 0,)), vm.x3D,
                             piston, vm.Point3D((0 , 0, 0)),
                             name='rod_piston linkage')

piston_ground = PrismaticLinkage(ground, vm.Point3D((0, 0, 0.)), vm.z3D,
                                 piston, vm.Point3D((0., 0, 0, )),
                                 name='piston_ground linkage')

mechanism = MovingMechanism([crank_ground, crank_rod, rod_piston, piston_ground],
                            ground, name='Crank_rod_mechanism')

#manual_configuration = MechanismConfigurations(mechanism, [[3.1415/4, 3.1415/8, 0.2, -0.3]])
#manual_configuration.plot2D(x=vm.y3D, y=vm.z3D)


configuration = mechanism.solve_configurations({0: npy.arange(0, 2*3.14, 0.1)})
configuration.plot2D(x=vm.y3D, y=vm.z3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.y3D, plot_frames=False)
#configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)

configuration.plot_kinematic_parameters(crank_ground, 0, piston_ground, 0)

configuration.babylonjs(plot_frames=True)
