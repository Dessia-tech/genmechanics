# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:14:34 2016

@author: Steven Masfaraud
"""

import numpy as npy
from math import cos,sin

class Linkage:
    def __init__(self,part1,part2,position,euler_angles,static_matrix1,static_matrix2,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 static_require_kinematic,name=''):
        self.part1=part1
        self.part2=part2
        self.position=position
        self.euler_angles=euler_angles
        self.name=name        
        
        self.static_matrix1=static_matrix1
        self.static_matrix2=static_matrix2
        self.static_behavior_occurence_matrix=static_behavior_occurence_matrix
        
        self.static_behavior_nonlinear_eq_indices=static_behavior_nonlinear_eq_indices
        self.static_behavior_linear_eq=static_behavior_linear_eq
        self.static_behavior_nonlinear_eq=static_behavior_nonlinear_eq
        self.static_require_kinematic=static_require_kinematic
        
        self.n_static_unknowns=static_matrix1.shape[1]
        

        

class HolonomicLinkage(Linkage):
    holonomic=True
    def __init__(self,part1,part2,position,euler_angles,static_matrix1,static_matrix2,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 kinematic_matrix,static_require_kinematic=False,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix1,static_matrix2,
                         static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                         static_behavior_linear_eq,static_behavior_nonlinear_eq,
                         static_require_kinematic,name)
        self.kinematic_matrix=kinematic_matrix
        self.n_kinematic_unknowns=kinematic_matrix.shape[1]
        
    

    def LocalSpeeds(self,q):
        """
        :returns a (2,3) shape numpy arraywith speeds in local coordinate system. First row are rotatings speeds,
            the second linear speeds
        """
        return npy.dot(self.kinematic_matrix,q)

        
class NonHolonomicLinkage(Linkage):
    holonomic=False
    def __init__(self,part1,part2,position,euler_angles,static_matrix1,static_matrix2,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 kinematic_directions,static_require_kinematic=False,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix1,static_matrix2,
                         static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                         static_behavior_linear_eq,static_behavior_nonlinear_eq,
                         static_require_kinematic,name)
        self.kinematic_directions=kinematic_directions
        

        
class FrictionlessRevoluteLinkage(HolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix2=npy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
        static_matrix1=-static_matrix2
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        kinematic_matrix=npy.array([[1],[0],[0],[0],[0],[0]])
        static_require_kinematic=False
        HolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                  static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq,kinematic_matrix,static_require_kinematic,name)

#class CylindricalLinkage(Linkage):
#    def __init__(self,part1,part2,position,euler_angles,name=''):
#        static_matrix=npy.array([[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]])
#        static_linear_eq_indices=[True,True,True,True,True]
#        static_linear_eq=npy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
#        nonlinear_static_eq=[]
#        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_linear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)
#
#class PrismaticLinkage(Linkage):
#    def __init__(self,part1,part2,position,euler_angles,name=''):
#        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_linear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)
#

class FrictionlessBallLinkage(HolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix2=npy.array([[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,0],[0,0,0]])
        static_matrix1=-static_matrix2
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        kinematic_matrix=npy.array([[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,0],[0,0,0]])
        static_require_kinematic=False
        HolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                  static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq,kinematic_matrix,static_require_kinematic,name)
#
class FrictionlessLinearAnnularLinkage(HolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix2=npy.array([[0,0],[1,0],[0,1],[0,0],[0,0],[0,0]])
        static_matrix1=-static_matrix2
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        kinematic_matrix=npy.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0],[0,0,0,0]])
        static_require_kinematic=False
        HolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                  static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices,
                                  static_behavior_linear_eq,static_behavior_nonlinear_eq,
                                  kinematic_matrix,static_require_kinematic,name)

class FrictionlessGearSetLinkage(NonHolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,alpha,beta,name=''):
        self.alpha=alpha
        self.beta=beta
        static_matrix2=npy.array([[cos(beta)*cos(alpha),0],[sin(alpha),0],[0,-1],[0,0],[0,0],[0,0]])
        static_matrix1=-static_matrix2
        static_behavior_occurence_matrix=npy.array([[1,1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[lambda x:abs(sin(alpha)*cos(beta)*x[0])-x[1]]
        directions=[npy.array([1,0,0])]
        static_require_kinematic=False
        NonHolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                     static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                     static_behavior_nonlinear_eq_indices,
                                     static_behavior_linear_eq,static_behavior_nonlinear_eq,
                                     directions,static_require_kinematic,name)

class BallLinkage(HolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,Ca,Cr,name=''):
        self.Ca=Ca
        self.Cr=Cr
        static_matrix2=npy.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0],[0,0,0,0]])
        static_matrix1=-static_matrix2
        static_behavior_occurence_matrix=npy.array([[1,1,1,1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[lambda x,v:abs(v[0])/v[0]*(Ca*abs(x[0])+Cr*(x[1]**2+x[2]**2)**0.5)+x[3] if v[0]!=0 else 0.]
        kinematic_matrix=npy.array([[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,0],[0,0,0]])
        static_require_kinematic=True
        HolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                  static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq,kinematic_matrix,
                                  static_require_kinematic,name)

class LinearAnnularLinkage(HolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,Cr,name=''):
        self.Cr=Cr
        static_matrix2=npy.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,0]])
        static_matrix1=-static_matrix2
        static_behavior_occurence_matrix=npy.array([[1,1,1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[lambda x,v:abs(v[0])/v[0]*(Cr*(x[0]**2+x[1]**2)**0.5)+x[2] if v[0]!=0 else 0.]
        kinematic_matrix=npy.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0],[0,0,0,0]])
        static_require_kinematic=True
        HolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                  static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                  static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,
                                  static_behavior_nonlinear_eq,kinematic_matrix,
                                  static_require_kinematic,name)

class GearSetLinkage(NonHolonomicLinkage):
    def __init__(self,part1,part2,position,euler_angles,alpha,beta,cf,name=''):
        self.alpha=alpha
        self.beta=beta
        self.cf=cf# force coefficient
        static_matrix2=npy.array([[cos(beta)*cos(alpha),0,0],[sin(alpha),0,0],[0,0,-1],[0,0,0],[0,0,0],[0,0,0]])
        static_matrix1=npy.array([[0,cos(beta)*cos(alpha),0],[0,sin(alpha),0],[0,0,1],[0,0,0],[0,0,0],[0,0,0]])
        static_behavior_occurence_matrix=npy.array([[1,1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[lambda x:abs(sin(alpha)*cos(beta))*max(x[0],x[1])+x[2],
                                      lambda x,v: x[1]-x[0]*(self.cf*(1+sin(alpha)**2*cos(alpha)**2)**0.5-1)
                                      if v[3]*x[0]>0 else x[0]-x[1]*(self.cf*(1+sin(alpha)**2*cos(alpha)**2)**0.5-1)]
        directions=[npy.array([1,0,0])]
        static_require_kinematic=False
        NonHolonomicLinkage.__init__(self,part1,part2,position,euler_angles,
                                     static_matrix1,static_matrix2,static_behavior_occurence_matrix,
                                     static_behavior_nonlinear_eq_indices,
                                     static_behavior_linear_eq,static_behavior_nonlinear_eq,
                                     directions,static_require_kinematic,name)
        
