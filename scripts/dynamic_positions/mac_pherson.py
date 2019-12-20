#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

simple mac pherson exemple

"""
import numpy as npy
import volmdlr as vm
from genmechanics.dynamic_positions import Part, RevoluteLinkage, MovingMechanism, BallLinkage, SlidingRevoluteLinkage, PrismaticLinkage


Z_UBJ = 0.4
Y_UBJ = 0.08

HUB_HEIGHT = 0.20

TRIANGLE_WIDTH = 0.3
TRIANGLE_BALLJOINT_X = -0.11

Y_CONTACT_TYRE = 0.4
Z_CONTACT_TYRE = -0.25


car_body = Part('Car body')
triangle = Part('Triangle')
strut = Part('Strut')
hub = Part('Wheel hub')


triangle_body_joint = RevoluteLinkage(car_body, vm.Point3D((0, 0, 0)), vm.x3D, triangle, vm.Point3D((0, 0, 0)), 'triangle body joint')
upper_ball_joint = BallLinkage(car_body, vm.Point3D((0, Y_UBJ, Z_UBJ)), strut, vm.Point3D((0, 0, 0)), 'strut body joint')
triangle_ball_joint = BallLinkage(triangle, vm.Point3D((TRIANGLE_BALLJOINT_X, TRIANGLE_WIDTH, 0,)), hub, vm.Point3D((HUB_HEIGHT ,0, 0)), 'triangle_ball_joint')

#strut_axis = -vm.z3D.Rotation(vm.o3D, vm.x3D, theta)
strut_joint = PrismaticLinkage(strut, vm.Point3D((0.2, 0, 0, )), vm.x3D, hub, vm.Point3D((0, 0, 0)), 'strut_joint')

mechanism = MovingMechanism([triangle_body_joint, upper_ball_joint, triangle_ball_joint, strut_joint], car_body, 'Mac Pherson')

configuration = mechanism.solve_configurations({0: [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.05, 0.1, 0.15, 0.3, 0.4, 0.5, 0.6]})

# Determining contact point at curb weight
hub_frame_curb_weight = mechanism.part_frame(hub, configuration.steps[6])
hub.interest_points.append(hub_frame_curb_weight.NewCoordinates(vm.Point3D((0, Y_CONTACT_TYRE, 0))))
hub.interest_points.append(hub_frame_curb_weight.NewCoordinates(vm.Point3D((0, Y_CONTACT_TYRE, Z_CONTACT_TYRE))))

configuration.plot2D(x=vm.y3D, y=vm.z3D, plot_frames=False)
configuration.plot2D(x=vm.x3D, y=vm.z3D, plot_frames=False)

configuration.plot_kinematic_parameters(triangle_body_joint, 0, strut_joint, 0)

configuration.plot2D_trajectory(hub.interest_points[0], hub, car_body)
configuration.plot_trajectory(hub.interest_points[0], hub, car_body)
configuration.plot_trajectory(hub.interest_points[1], hub, car_body)

configuration.babylonjs(plot_frames=False, plot_trajectories=True)
