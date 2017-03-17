#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:15:09 2017


"""

from numpy import array,zeros
import genmechanics.geometry as geometry

class KnownLoad:
    def __init__(self,part,position,euler_angles,forces,torques,name=''):
        self.part=part
        self.position=position
        self.euler_angles=euler_angles
        self.forces=array(forces)
        self.torques=array(torques)
        self.name=name
        
        self.P=geometry.Euler2TransferMatrix(*self.euler_angles) 
        
class UnknownLoad:
    """
    :param force_directions: a list of directions for force (0,1,2)
    :param torque_directions: a list of directions for torque (0,1,2)
    """
    def __init__(self,part,position,euler_angles,static_matrix,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 static_require_kinematic,name=''):
        
        self.part=part
        self.position=position
        self.euler_angles=euler_angles

        self.static_matrix=static_matrix

        self.static_behavior_occurence_matrix=static_behavior_occurence_matrix
        
        self.static_behavior_nonlinear_eq_indices=static_behavior_nonlinear_eq_indices
        self.static_behavior_linear_eq=static_behavior_linear_eq
        self.static_behavior_nonlinear_eq=static_behavior_nonlinear_eq
        self.static_require_kinematic=static_require_kinematic

        self.name=name
#        print(euler_angles)
        self.P=geometry.Euler2TransferMatrix(*self.euler_angles) 
        self.n_static_unknowns=self.static_matrix.shape[1]


class SimpleUnknownLoad(UnknownLoad):
    """
    :param force_directions: a list of directions for force (0,1,2)
    :param torque_directions: a list of directions for torque (0,1,2)
    """
    def __init__(self,part,position,euler_angles,force_directions,torque_directions,name=''):
        
#        print(euler_angles)

        self.force_directions=force_directions
        self.torque_directions=torque_directions
        lfd=len(force_directions)
        ltd=len(torque_directions)
        static_matrix=zeros((6,lfd+ltd))
        for i,k in enumerate(force_directions):
            static_matrix[k,i]=1
        for i,k in enumerate(torque_directions):
            static_matrix[k+3,i+lfd]=1
            
        static_behavior_occurence_matrix=array([])
        static_behavior_occurence_matrix=array([])
        static_behavior_nonlinear_eq_indices=[]
        static_behavior_linear_eq=array([])
        static_behavior_nonlinear_eq=[]

        static_require_kinematic=False

        UnknownLoad.__init__(self,part,position,euler_angles,static_matrix,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 static_require_kinematic,name)
        

class SplashLoad(UnknownLoad):
    """
    Creates a splash force linked to the rotationnal speed around
    """
    def __init__(self,part,position,euler_angles,area,radius,Cl,Ct,Rl,mu,name):
        self.Cl=Cl
        self.Ct=Ct
        self.Rl=Rl# limit reynods number
        self.radius=radius
        self.mu=mu
        self.area=area
        
        static_matrix=zeros((6,1))
        static_matrix[3,0]=1# Resistant torque on X
        static_behavior_occurence_matrix=array([[1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=array([])
        static_behavior_nonlinear_eq=[lambda x,w,v:x[0]+self.Cl*self.area*w[0]**2/abs(w[0]) if self.radius**2*w[0]/self.mu<self.Rl else x[0]+abs(w[0])/w[0]*self.Ct*self.area]
        static_require_kinematic=True
        
        UnknownLoad.__init__(self,part,position,euler_angles,static_matrix,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 static_require_kinematic,name)        
        
        
        