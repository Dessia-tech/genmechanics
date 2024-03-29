# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:14:34 2016

@author: Steven Masfaraud
"""

import numpy as npy
from math import cos, sin, tan, atan
import genmechanics.geometry as geometry
from dessia_common.core import DessiaObject


class Linkage:
    _non_serializable_attributes = ['position', 'euler_angles',
                                    'static_matrix1', 'static_matrix2',
                                    'static_behavior_occurence_matrix',
                                    'static_behavior_nonlinear_eq_indices',
                                    'static_behavior_linear_eq',
                                    'static_behavior_nonlinear_eq',
                                    'static_require_kinematic', 'P', 'kinematic_matrix']

    def __init__(self,
                 part1, part2, position,
                 euler_angles,
                 static_matrix1, static_matrix2,
                 static_behavior_occurence_matrix,
                 static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,
                 static_behavior_nonlinear_eq,
                 static_require_kinematic, name=''):
        self.part1 = part1
        self.part2 = part2
        self.position = npy.array(position)
        self.euler_angles = euler_angles
        self.name = name

        self.static_matrix1 = static_matrix1
        self.static_matrix2 = static_matrix2
        self.static_behavior_occurence_matrix = static_behavior_occurence_matrix

        self.static_behavior_nonlinear_eq_indices = static_behavior_nonlinear_eq_indices
        self.static_behavior_linear_eq = static_behavior_linear_eq
        self.static_behavior_nonlinear_eq = static_behavior_nonlinear_eq
        self.static_require_kinematic = static_require_kinematic

        self.n_static_unknowns = static_matrix1.shape[1]

        self.P = geometry.euler_2_transfer_matrix(*self.euler_angles)

    def babylon(self, length, forces, torques):
        xa, za, ya = self.P[:, 0]
        #        theta=math.acos(z/self.width)
        #        phi=math.atan(y/x)
        x, z, y = self.position
        s = """
        var sphere = BABYLON.Mesh.CreateSphere("{} center", 8, {}, scene);
        sphere.position=new BABYLON.Vector3({},{},{});
        var lineX = BABYLON.Mesh.CreateDashedLines("{} axis", [
            new BABYLON.Vector3({}, {}, {}),
            new BABYLON.Vector3({}, {}, {})
        ],0.05,0.05,10, scene);
        """.format(self.name, length / 30, x, y, z, self.name, x - 0.5 * xa, y - 0.5 * ya, z - 0.5 * za, x + 0.5 * xa,
                   y + 0.5 * ya, z + 0.5 * za)

        if forces is not None:
            # print(forces)
            # for i,fi in enumerate(forces):
            s += """var lineF = BABYLON.Mesh.CreateLines("{} axis", [
            new BABYLON.Vector3({}, {}, {}),
            new BABYLON.Vector3({}, {}, {})
        ], scene);
            lineF.color=new BABYLON.Color3(1,0,0);""".format('force {}'.format(self.name), x, y, z, x + forces[0],
                                                             y + forces[1], z + forces[2])

        return s


class HolonomicLinkage(Linkage):
    holonomic = True

    def __init__(self, part1, part2, position, euler_angles, static_matrix1, static_matrix2,
                 static_behavior_occurence_matrix, static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq, static_behavior_nonlinear_eq,
                 kinematic_matrix, static_require_kinematic=False, name=''):
        Linkage.__init__(self, part1, part2, position, euler_angles, static_matrix1, static_matrix2,
                         static_behavior_occurence_matrix, static_behavior_nonlinear_eq_indices,
                         static_behavior_linear_eq, static_behavior_nonlinear_eq,
                         static_require_kinematic, name)
        self.kinematic_matrix = kinematic_matrix
        self.n_kinematic_unknowns = kinematic_matrix.shape[1]


class NonHolonomicLinkage(Linkage):
    holonomic = False

    def __init__(self,
                 part1, part2, position, euler_angles,
                 static_matrix1, static_matrix2,
                 static_behavior_occurence_matrix,
                 static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq, static_behavior_nonlinear_eq,
                 kinematic_directions,
                 static_require_kinematic=False, name=''):
        Linkage.__init__(self, part1, part2, position, euler_angles, static_matrix1, static_matrix2,
                         static_behavior_occurence_matrix, static_behavior_nonlinear_eq_indices,
                         static_behavior_linear_eq, static_behavior_nonlinear_eq,
                         static_require_kinematic, name)
        self.kinematic_directions = kinematic_directions


class FrictionlessRevoluteLinkage(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, name='Frictionless Revolute Linkage'):
        static_matrix2 = npy.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 0],
                                    [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([])
        static_behavior_nonlinear_eq_indices = []
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = []
        kinematic_matrix = npy.array([[1], [0], [0], [0], [0], [0]])
        static_require_kinematic = False
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices, static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq, kinematic_matrix, static_require_kinematic, name)
        # print(self.to_dict())
        # print('EZA')


class RevoluteLinkage(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, Ca, Cr, Cw, name='Revolute Linkage'):
        self.Ca = Ca
        self.Cr = Cr
        self.Cw = Cw
        static_matrix2 = npy.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                                    [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([[1, 1, 1, 1, 0, 0]])
        static_behavior_nonlinear_eq_indices = [0]
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = [lambda x, w, v: abs(
            w[0]) / w[0] * (Ca * abs(x[0]) + Cr * (x[1] ** 2 + x[2] ** 2) ** 0.5 + Cw * w[0]) + x[3] if w[0] != 0 else
            x[3]]
        kinematic_matrix = npy.array([[1], [0], [0], [0], [0], [0]])
        static_require_kinematic = True
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices, static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq, kinematic_matrix,
                                  static_require_kinematic, name)

    def change_coefficients(self, Ca, Cr, Cw):
        self.Ca = Ca
        self.Cr = Cr
        self.Cw = Cw
        self.static_behavior_nonlinear_eq = [lambda x, w, v: abs(
            w[0]) / w[0] * (Ca * abs(x[0]) + Cr * (x[1] ** 2 + x[2] ** 2) ** 0.5 + Cw * w[0]) + x[3] if w[0] != 0 else
            x[3]]


# class CylindricalLinkage(Linkage):
#    def __init__(self,part1,part2,position,euler_angles,name=''):
#        static_matrix=npy.array([[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]])
#        static_linear_eq_indices=[True,True,True,True,True]
#        static_linear_eq=npy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
#        nonlinear_static_eq=[]
#        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_linear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)
#
# class PrismaticLinkage(Linkage):
#    def __init__(self,part1,part2,position,euler_angles,name=''):
#        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_linear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)
#

class FrictionlessBallLinkage(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, name=''):
        static_matrix2 = npy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([])
        static_behavior_nonlinear_eq_indices = []
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = []
        kinematic_matrix = npy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        static_require_kinematic = False
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices, static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq, kinematic_matrix, static_require_kinematic, name)


#


class FrictionlessLinearAnnularLinkage(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, name='Frictionless Linear Annular Linkage'):
        static_matrix2 = npy.array([[0, 0], [1, 0], [0, 1], [0, 0], [0, 0], [0, 0]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([])
        static_behavior_nonlinear_eq_indices = []
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = []
        kinematic_matrix = npy.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                      [0, 0, 0, 0], [0, 0, 0, 0]])
        static_require_kinematic = False
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices,
                                  static_behavior_linear_eq, static_behavior_nonlinear_eq,
                                  kinematic_matrix, static_require_kinematic, name)


class BallLinkage(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, Ca, Cr, Cw, name='Ball Linkage'):
        self.Ca = Ca
        self.Cr = Cr
        self.Cw = Cw
        static_matrix2 = npy.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                    [0, 0, 0, 0], [0, 0, 0, 0]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([[1, 1, 1, 1]])
        static_behavior_nonlinear_eq_indices = [0]
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = [lambda x, w, v: abs(
            w[0]) / w[0] * (Ca * abs(x[0]) + Cr * (x[1] ** 2 + x[2] ** 2) ** 0.5 + Cw * w[0]) + x[3] if w[0] != 0 else
            x[3]]
        kinematic_matrix = npy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        static_require_kinematic = True
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices, static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq, kinematic_matrix,
                                  static_require_kinematic, name)

    def ChangeCoefficients(self, Ca, Cr, Cw):
        self.Ca = Ca
        self.Cr = Cr
        self.Cw = Cw
        self.static_behavior_nonlinear_eq = [lambda x, w, v: abs(
            w[0]) / w[0] * (Ca * abs(x[0]) + Cr * (x[1] ** 2 + x[2] ** 2) ** 0.5 + Cw * w[0]) + x[3] if w[0] != 0 else
            x[3]]


class LinearAnnularLinkage(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, Cr, Cw, name='Linear Annular Linkage'):
        self.Cr = Cr
        self.Cw = Cw
        static_matrix2 = npy.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0], [0, 0, 0]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([[1, 1, 1]])
        static_behavior_nonlinear_eq_indices = [0]
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = [lambda x, w, v: abs(
            w[0]) / w[0] * (Cr * (w[0] ** 2 + x[1] ** 2) ** 0.5 + Cw * w[0]) + x[2] if w[0] != 0 else x[2]]
        kinematic_matrix = npy.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                      [0, 0, 0, 0], [0, 0, 0, 0]])
        static_require_kinematic = True
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices, static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq, kinematic_matrix,
                                  static_require_kinematic, name)

    def change_coefficients(self, Cr, Cw):
        self.Cr = Cr
        self.Cw = Cw
        self.static_behavior_nonlinear_eq = [lambda x, w, v: abs(
            w[0]) / w[0] * (Cr * (w[0] ** 2 + x[1] ** 2) ** 0.5 + Cw * w[0]) + x[2] if w[0] != 0 else x[2]]


# class AxialStop(HolonomicLinkage):
#    def __init__(self,part1,part2,position,euler_angles,Ca,Cw,name='Axial Stop'):
#        self.Ca=Ca
#        self.Cw=Cw
#        static_matrix2=npy.array([[1,0],[0,0],[0,0],[0,1],[0,0],[0,0]])
#        static_matrix1=-static_matrix2
#        static_behavior_occurence_matrix=npy.array([[1,1]])
#        static_behavior_nonlinear_eq_indices=[0]
#        static_behavior_linear_eq=npy.array([])
#        static_behavior_nonlinear_eq=[lambda x,w,v:abs(w[0])/w[0]*Ca*abs(x[0])+x[1] if w[0]!=0 else x[1]]
#        kinematic_matrix=npy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
#        static_require_kinematic=True
#        HolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
#                                  static_matrix1,static_matrix2,static_behavior_occurence_matrix,
#                                  static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,
#                                  static_behavior_nonlinear_eq,kinematic_matrix,
#                                  static_require_kinematic,name)
#
#
#    def ChangeCoefficients(self,Ca,Cw):
#        self.Ca=Ca
#        self.Cw=Cw
#        self.static_behavior_nonlinear_eq=[lambda x,w,v:abs(w[0])/w[0]*Ca*abs(x[0])+x[1] if w[0]!=0 else x[1]]
#

class RotationalStop(HolonomicLinkage):
    def __init__(self, part1, part2, position, euler_angles, name='Rotational Stop'):
        static_matrix2 = npy.array([[0], [0], [0], [1], [0], [0]])
        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([])
        static_behavior_nonlinear_eq_indices = []
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = []
        kinematic_matrix = npy.array([[0, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0],
                                      [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
        static_require_kinematic = True
        HolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                  static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices, static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq, kinematic_matrix,
                                  static_require_kinematic, name)


class FrictionlessGearSetLinkage(NonHolonomicLinkage):
    """
    :param pressure_angle: pressure angle
    :param helix_angle: helix angle
    """

    def __init__(self, part1, part2, position, radial_vector, axial_vector,
                 pressure_angle, helix_angle, name='Gear Set Linkage'):
        self.pressure_angle = pressure_angle  # pressure angle
        self.helix_angle = helix_angle

        transversal_vector = npy.cross(axial_vector, radial_vector)
        euler_angles = geometry.direction_2_euler(transversal_vector, axial_vector)
        static_matrix2 = npy.array([[1, 0], [tan(helix_angle), 0], [0, -1], [0, 0], [0, 0], [0, 0]])

        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([[1, 1]])
        static_behavior_nonlinear_eq_indices = [0]
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = [lambda x, w, v: tan(pressure_angle) * abs(x[0]) / cos(helix_angle) + x[1]]
        directions = [npy.array([1, 0, 0])]
        static_require_kinematic = True
        NonHolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                     static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                     static_behavior_nonlinear_eq_indices,
                                     static_behavior_linear_eq, static_behavior_nonlinear_eq,
                                     directions, static_require_kinematic, name)

    def change_parameters(self, pressure_angle, helix_angle):
        self.pressure_angle = pressure_angle
        self.helix_angle = helix_angle
        self.static_matrix2 = npy.array([[1, 0], [tan(helix_angle), 0], [0, -1], [0, 0], [0, 0], [0, 0]])
        self.static_matrix1 = - self.static_matrix2
        self.static_behavior_nonlinear_eq = [lambda x, w, v: tan(pressure_angle) * abs(x[0]) / cos(helix_angle) + x[1]]


class FrictionlessGearSetLinkage(NonHolonomicLinkage):
    """
    :param pressure_angle: pressure angle
    :param helix_angle: helix angle
    """

    def __init__(self, part1, part2, position, radial_vector, axial_vector,
                 pressure_angle, helix_angle, name='Gear Set Linkage'):
        self.pressure_angle = pressure_angle  # pressure angle
        self.helix_angle = helix_angle

        transversal_vector = npy.cross(axial_vector, radial_vector)
        euler_angles = geometry.direction_2_euler(transversal_vector, axial_vector)
        static_matrix2 = npy.array([[1, 0], [tan(helix_angle), 0], [0, -1], [0, 0], [0, 0], [0, 0]])

        static_matrix1 = -static_matrix2
        static_behavior_occurence_matrix = npy.array([[1, 1]])
        static_behavior_nonlinear_eq_indices = [0]
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = [lambda x, w, v: tan(pressure_angle) * abs(x[0]) * cos(helix_angle) + x[1]]
        directions = [npy.array([1, 0, 0])]
        static_require_kinematic = True
        NonHolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                     static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                     static_behavior_nonlinear_eq_indices,
                                     static_behavior_linear_eq, static_behavior_nonlinear_eq,
                                     directions, static_require_kinematic, name)

    def change_parameters(self, pressure_angle, helix_angle):
        self.pressure_angle = pressure_angle
        self.helix_angle = helix_angle
        self.static_matrix2 = npy.array([[1, 0], [tan(helix_angle), 0], [0, -1], [0, 0], [0, 0], [0, 0]])
        self.static_matrix1 = - self.static_matrix2
        self.static_behavior_nonlinear_eq = [lambda x, w, v: tan(pressure_angle) * abs(x[0]) * cos(helix_angle) + x[1]]


class GearSetLinkage(NonHolonomicLinkage):
    """
    :param alpha: pressure angle
    :param beta: helix angle
    """

    def __init__(self, part1, part2, position, radial_vector, axial_vector, pressure_angle, helix_angle,
                 Cf, Cv, name='Gear Set Linkage'):
        self.pressure_angle = pressure_angle  # pressure angle
        self.helix_angle = helix_angle  # hélix angle
        self.Cf = Cf  # force coefficient
        self.Cv = Cv  # Speed coefficient

        transversal_vector = npy.cross(axial_vector, radial_vector)
        euler_angles = geometry.direction_2_euler(transversal_vector, axial_vector)

        static_matrix2 = npy.array([[1, 0, 0], [tan(helix_angle), 0, 0], [0, 0, -1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        static_matrix1 = npy.array([[0, 1, 0], [0, tan(helix_angle), 0], [0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        static_behavior_occurence_matrix = npy.array([[1, 1, 1], [1, 1, 0]])
        static_behavior_nonlinear_eq_indices = [0, 1]
        static_behavior_linear_eq = npy.array([])
        #        static_behavior_nonlinear_eq=[lambda x,w,v:abs(sin(beta)*cos(alpha)*max(abs(x[0]),abs(x[1])))+x[2],
        #                                      lambda x,w,v: x[1]-x[0]*(Cf*(1+sin(beta)**2*cos(alpha)**2)**0.5-1)+Cv*abs(v[0])
        # if v[0]*x[0]>0 else
        # x[0]-x[1]*(Cf*(1+sin(beta)**2*cos(alpha)**2)**0.5-1)+Cv*abs(v[0])]
        static_behavior_nonlinear_eq = self.update_behavior()

        directions = [npy.array([1, 0, 0])]
        static_require_kinematic = True
        NonHolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                     static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                     static_behavior_nonlinear_eq_indices,
                                     static_behavior_linear_eq, static_behavior_nonlinear_eq,
                                     directions, static_require_kinematic, name)

    def update_behavior(self):
        f1 = lambda x, w, v: tan(self.pressure_angle) * max(abs(x[1]), abs(x[0])) * cos(self.helix_angle) + x[2]

        def f2(x, w, v):
            if x[0] == 0.:
                if v[0] * x[1] > 0:
                    r = x[1] - (x[0] * (self.Cf * (1 + (tan(self.pressure_angle) ** 2 + sin(self.helix_angle)
                                                        ** 2) / cos(self.helix_angle)) ** 0.5 - 1) + self.Cv * abs(
                        v[0]))
                else:
                    r = x[0] - (x[1] * (self.Cf * (1 + (tan(self.pressure_angle) ** 2 + sin(self.helix_angle)
                                                        ** 2) / cos(self.helix_angle)) ** 0.5 - 1) + self.Cv * abs(
                        v[0]))
            else:
                if v[0] * x[0] < 0:
                    r = x[1] - (x[0] * (self.Cf * (1 + (tan(self.pressure_angle) ** 2 + sin(self.helix_angle)
                                                        ** 2) / cos(self.helix_angle)) ** 0.5 - 1) + self.Cv * abs(
                        v[0]))
                else:
                    r = x[0] - (x[1] * (self.Cf * (1 + (tan(self.pressure_angle) ** 2 + sin(self.helix_angle)
                                                        ** 2) / cos(self.helix_angle)) ** 0.5 - 1) + self.Cv * abs(
                        v[0]))
            return r

        return [f1, f2]

    def change_coefficients(self, Cf, Cv):
        self.Cf = Cf
        self.Cv = Cv

        self.static_behavior_nonlinear_eq = self.update_behavior()

    def change_parameters(self, pressure_angle, helix_angle):
        self.pressure_angle = pressure_angle
        self.helix_angle = helix_angle
        self.static_matrix2 = npy.array([[1, 0, 0], [tan(helix_angle), 0, 0], [
            0, 0, -1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.static_matrix1 = npy.array([[0, 1, 0], [0, tan(helix_angle), 0], [
            0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]])

        self.static_behavior_nonlinear_eq = self.update_behavior()


class FrictionlessBevelGearLinkage(NonHolonomicLinkage):
    """
    :param pressure_angle: pressure angle
    :param helix_angle: helix angle
    """

    def __init__(self, part1, part2, position, radial_vector_gear_part_1, axial_vector_gear_part_1,
                 pressure_angle, mean_spiral_angle, pitch_angle_gear_part_1, hand_of_spiral_gear_part_1='LH',
                 name='Bevel Gear Linkage'):
        self.pressure_angle = pressure_angle  # pressure angle
        self.mean_spiral_angle = mean_spiral_angle
        self.pitch_angle_gear_part_1 = pitch_angle_gear_part_1
        self.hand_of_spiral_gear_part_1 = hand_of_spiral_gear_part_1
        self.radial_vector_gear_part_1 = radial_vector_gear_part_1
        self.axial_vector_gear_part_1 = axial_vector_gear_part_1
        normal_pressure_angle = atan(tan(pressure_angle) * cos(mean_spiral_angle))
        if hand_of_spiral_gear_part_1 == 'LH':
            hand_of_spiral_coefficient = 1
        else:
            hand_of_spiral_coefficient = -1

        transversal_vector = npy.cross(axial_vector_gear_part_1, radial_vector_gear_part_1)

        euler_angles = geometry.direction_2_euler(transversal_vector, axial_vector_gear_part_1)

        static_matrix1 = npy.array([[1, 0, 0],
                                    [0, tan(pitch_angle_gear_part_1), 1 / tan(pitch_angle_gear_part_1)],
                                    [0, 1, -1],
                                    [0, 0, 0], [0, 0, 0], [0, 0, 0]])

        static_matrix2 = -static_matrix1
        static_behavior_occurence_matrix = npy.array([[1, 1, 0], [1, 0, 1]])
        static_behavior_nonlinear_eq_indices = [0, 1]
        static_behavior_linear_eq = npy.array([])
        static_behavior_nonlinear_eq = [lambda x, w, v: (tan(normal_pressure_angle) *
                                                         cos(pitch_angle_gear_part_1)) *
                                        abs(x[0]) / cos(mean_spiral_angle) - x[1],
                                        lambda x, w, v: x[0] * hand_of_spiral_coefficient *
                                        sin(mean_spiral_angle) * sin(pitch_angle_gear_part_1) - x[1] *
                                        cos(mean_spiral_angle)]
        directions = [npy.array([1, 0, 0])]
        static_require_kinematic = True
        NonHolonomicLinkage.__init__(self, part1, part2, position, euler_angles,
                                     static_matrix1, static_matrix2, static_behavior_occurence_matrix,
                                     static_behavior_nonlinear_eq_indices,
                                     static_behavior_linear_eq, static_behavior_nonlinear_eq,
                                     directions, static_require_kinematic, name)

    def change_parameters(self, pressure_angle, mean_spiral_angle, pitch_angle_gear_part_1,
                          hand_of_spiral_gear_part_1='LH'):
        self.pressure_angle = pressure_angle
        self.mean_spiral_angle = mean_spiral_angle
        self.pitch_angle_gear_part_1 = pitch_angle_gear_part_1
        self.hand_of_spiral_gear_part_1 = hand_of_spiral_gear_part_1
        self.static_matrix1 = npy.array([[1, 0, 0],
                                         [0, tan(pitch_angle_gear_part_1), 1 / tan(pitch_angle_gear_part_1)],
                                         [0, 1, -1],
                                         [0, 0, 0], [0, 0, 0], [0, 0, 0]])

        self.static_matrix2 = -self.static_matrix1
        normal_pressure_angle = atan(tan(pressure_angle) * cos(mean_spiral_angle))
        if hand_of_spiral_gear_part_1 == 'LH':
            hand_of_spiral_coefficient = 1
        else:
            hand_of_spiral_coefficient = -1
        self.static_behavior_nonlinear_eq = [lambda x, w, v: (tan(normal_pressure_angle) *
                                                              cos(pitch_angle_gear_part_1)) *
                                             abs(x[0]) / cos(mean_spiral_angle) - x[1],
                                             lambda x, w, v: x[0] * hand_of_spiral_coefficient *
                                             sin(mean_spiral_angle) * sin(pitch_angle_gear_part_1) - x[
                                                                 1] *
                                             cos(mean_spiral_angle)]

# class FrictionlessGearSetLinkage(NonHolonomicLinkage):
#     """
#     :param alpha: pressure angle
#     :param beta: helix angle
#     """
#     def __init__(self,part1,part2,position,euler_angles,alpha,beta,name='Gear Set Linkage'):
#         self.alpha=alpha# pressure angle
#         self.beta=beta# hélix angle
#
#         static_matrix2=npy.array([[cos(beta)*cos(alpha),0],[sin(beta),0],[0,-1],[0,0],[0,0],[0,0]])
#         static_matrix1=-static_matrix2
#         static_behavior_occurence_matrix=npy.array([[1,1,1]])
#         static_behavior_nonlinear_eq_indices=[0]
#         static_behavior_linear_eq=npy.array([])
#         static_behavior_nonlinear_eq=[lambda x,w,v:abs(sin(beta)*cos(alpha)*abs(x[0]))+x[1]]
#         directions=[npy.array([1,0,0])]
#         static_require_kinematic=True
#         NonHolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
#                                      static_matrix1,static_matrix2,static_behavior_occurence_matrix,
#                                      static_behavior_nonlinear_eq_indices,
#                                      static_behavior_linear_eq,static_behavior_nonlinear_eq,
#                                      directions,static_require_kinematic,name)

#    def ChangeCoefficients(self,Cf,Cv):
#        self.Cf=Cf
#        self.Cv=Cv
#        self.static_behavior_nonlinear_eq=[lambda x,w,v:abs(sin(self.beta)*cos(self.alpha)*max(abs(x[0]),abs(x[1])))+x[2],
#                                      lambda x,w,v: x[1]-x[0]*(Cf*(1+sin(self.beta)**2*cos(self.alpha)**2)**0.5-1)+Cv*abs(v[0])
# if v[0]*x[0]>0 else
# x[0]-x[1]*(Cf*(1+sin(self.beta)**2*cos(self.alpha)**2)**0.5-1)+Cv*abs(v[0])]

#    def ChangeParameters(self,alpha,beta):
#        self.alpha=alpha
#        self.beta=beta
#        self.static_matrix2=npy.array([[cos(beta)*cos(alpha),0,0],[sin(beta),0,0],[0,0,-1],[0,0,0],[0,0,0],[0,0,0]])
#        self.static_matrix1=npy.array([[0,cos(beta)*cos(alpha),0],[0,sin(beta),0],[0,0,1],[0,0,0],[0,0,0],[0,0,0]])
#        self.static_behavior_nonlinear_eq=[lambda x,w,v:abs(sin(beta)*cos(alpha)*max(abs(x[0]),abs(x[1])))+x[2],
#                                      lambda x,w,v: x[1]-x[0]*(self.Cf*(1+sin(beta)**2*cos(alpha)**2)**0.5-1)+self.Cv*abs(v[0])
# if v[0]*x[0]>0 else
# x[0]-x[1]*(self.Cf*(1+sin(beta)**2*cos(alpha)**2)**0.5-1)+self.Cv*abs(v[0])]
