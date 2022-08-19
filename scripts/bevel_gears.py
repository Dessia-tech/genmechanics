# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:21:07 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads
import numpy as npy
import genmechanics.geometry as geometry
import math
import matplotlib.pyplot as plt
import scipy.linalg as linalg

# plt.figure()
L = 0.2
e1 = 0.1
e2 = 0.07
C = 300
w = 300
r_gear_1 = 0.2
r_gear_2 = 0.4
shaft_axe_1 = npy.array([1, 0, 0])
shaft_axe_2 = npy.array([0, 1, 0])

alpha_gs1 = 18/360*2*3.1415
beta_gs1 = 25/360*2*3.1415

ground = genmechanics.Part('ground')
shaft1 = genmechanics.Part('shaft1')
shaft2 = genmechanics.Part('shaft2')

radial_shaft_vector = npy.cross(shaft_axe_1, shaft_axe_2)
transversal_shaft_1_vector = npy.cross(radial_shaft_vector, shaft_axe_1)
transversal_shaft_2_vector = npy.cross(radial_shaft_vector, shaft_axe_2)
p1a = npy.array([0, 0, 0])
p1b = (L-0.1) * shaft_axe_1
p2a = L * shaft_axe_1 + e1 * shaft_axe_2 + r_gear_2*transversal_shaft_2_vector
p2b = p2a + L * shaft_axe_2

egs1 = geometry.direction_2_euler(shaft_axe_1, transversal_shaft_1_vector)

egs2 = geometry.direction_2_euler(shaft_axe_2, transversal_shaft_2_vector)
bearing1a = linkages.FrictionlessBallLinkage(ground, shaft1, p1a, egs1, 'bearing 1a')
bearing1b = linkages.FrictionlessLinearAnnularLinkage(ground, shaft1, p1b, egs1, 'bearing 1b')
bearing2a = linkages.FrictionlessBallLinkage(ground, shaft2, p2a, egs2, 'bearing 2a')
bearing2b = linkages.FrictionlessLinearAnnularLinkage(ground, shaft2, p2b, egs2, 'bearing 2b')




center_gearing = L*shaft_axe_1 + r_gear_1*transversal_shaft_1_vector
radial_vector = [0, 1, 0]
axial_vector = [-1, 0, 0]
transversal_vector = npy.cross(axial_vector, radial_vector)
pitch_angle_gear_part_1 = math.radians(30)
mean_spiral_angle = math.radians(20)
normal_pressure_angle = math.atan(math.tan(alpha_gs1) * math.cos(mean_spiral_angle))
gearset12 = linkages.FrictionlessBevelGearLinkage(shaft1, shaft2, center_gearing,
                                                  radial_vector_gear_part_1=radial_vector,
                                                  axial_vector_gear_part_1=axial_vector,
                                                  pressure_angle=alpha_gs1,
                                                  mean_spiral_angle=mean_spiral_angle,
                                                  pitch_angle_gear_part_1=pitch_angle_gear_part_1,
                                                  name='Gear set 1')
#gearset23=linkages.GearSetLinkage(shaft2,shaft3,[3*L/2,e1*r1+(e1+e2)*(1-r1),0],dir23)

imposed_speeds = [(bearing1a, 0, w)]

load1 = loads.KnownLoad(shaft1, (-L/4)*shaft_axe_1, [0, 0, 0], [0, 0, 0], [C, 0, 0], 'input torque')
load2 = loads.SimpleUnknownLoad(shaft2, p2b + 0.1 * shaft_axe_2, egs2, [], [0], 'output torque')

mech = genmechanics.Mechanism([bearing1a, bearing1b, bearing2a, bearing2b, gearset12],
                               ground, imposed_speeds, [load1], [load2])

#r=mech.StaticAnalysis()

for l,r in mech.static_results.items():
    for d,v in r.items(): 
        print(l.name,d,v)
    
    
print('Cth: ',-r_gear_1/(1-r_gear_1)*C)

Ft = C/r_gear_1
Fr = Ft * (math.tan(normal_pressure_angle) * math.cos(pitch_angle_gear_part_1) -
           math.sin(mean_spiral_angle) * math.sin(pitch_angle_gear_part_1))/math.cos(mean_spiral_angle)

Fa = Ft * (math.tan(normal_pressure_angle) * math.sin(pitch_angle_gear_part_1) +
           math.sin(mean_spiral_angle) * math.cos(pitch_angle_gear_part_1))/math.cos(mean_spiral_angle)

euler_angle = gearset12.euler_angles

transfert_matrix_gear = geometry.euler_2_transfer_matrix(*euler_angle)

F_gear_transfert = npy.dot(transfert_matrix_gear, npy.array([-Ft, -Fa, -Fr]))
print(transversal_vector)
print([Ft, Fa, Fr])
print(F_gear_transfert)
print(euler_angle)

bearing2a_result_force = mech.static_results[bearing2a]
bearing2b_result_force = mech.static_results[bearing2b]
bearing2a_force = npy.array([bearing2a_result_force[0], bearing2a_result_force[1], bearing2a_result_force[2]])
bearing2b_force = npy.array([0, bearing2b_result_force[0], bearing2b_result_force[1]])

transfert_matrix_bearing = geometry.euler_2_transfer_matrix(*egs2)
F_bearing_2a_transfert = npy.dot(transfert_matrix_bearing, bearing2a_force)

F_bearing_2b_transfert = npy.dot(transfert_matrix_bearing, bearing2b_force)
print(F_bearing_2a_transfert)
print(F_bearing_2b_transfert)
for i in range(3):
    result_F = F_gear_transfert[i] + F_bearing_2a_transfert[i] + F_bearing_2b_transfert[i]
    print('Test F sur axe '+ str(i))
    print(result_F)