# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:14:34 2016

@author: Steven Masfaraud
"""

from genmechanics import Linkage
import numpy as npy

class FrictionlessRevoluteLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix=npy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)

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
class FrictionlessBallLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix=npy.array([[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,0],[0,0,0]])
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)
#
class FrictionlessLinearAnnularLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix=npy.array([[0,0],[1,0],[0,1],[0,0],[0,0],[0,0]])
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)

class FrictionlessGearSetLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        static_matrix=npy.array([[1],[0],[0],[0],[0],[0]])
        static_behavior_occurence_matrix=npy.array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[]
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)

class BallLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,Ca,Cr,name=''):
        static_matrix=npy.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0],[0,0,0,0]])
        static_behavior_occurence_matrix=npy.array([[1,1,1,1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[lambda x:Ca*x[0]+Cr*(x[1]**2+x[2]**2)**0.5-x[3]]
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)

class LinearAnnularLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,Cr,name=''):
        static_matrix=npy.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0,0,0]])
        static_behavior_occurence_matrix=npy.array([[1,1,1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=npy.array([])
        static_behavior_nonlinear_eq=[lambda x:Cr*(x[0]**2+x[1]**2)**0.5-x[2]]
        Linkage.__init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name)
