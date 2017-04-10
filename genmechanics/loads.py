#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:15:09 2017


"""

from numpy import array,zeros
import genmechanics.geometry as geometry

import math

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
        

class GearSplashLoad(UnknownLoad):
    """
    Creates a splash force linked to the rotationnal speed around local X axis
    :param Rel: Reynolds limit
    :param d: caracteristic length. Default value: radius
    :param h: height of wet surface of gear
    """
    def __init__(self,part,position,euler_angles,h,radius,Cl,Ct,Rel,nu,d=None,name='Splash load'):
        self.Cl=Cl
        self.Ct=Ct
        self.Rel=Rel# limit reynods number
        self.radius=radius
        self.nu=nu
        
        if d==None:
            d=radius            
        self.d=d
            
        self.h=h

        if h<0:
            self.area=0
        elif h<self.radius:
            hp=self.radius-h
            self.area=self.radius**2*math.acos(hp/self.radius)-hp*math.sqrt(self.radius**2-hp**2)
        elif h<2*self.radius:
            h2=2*self.radius-h
            hp=self.radius-h2
            self.area=math.pi*self.radius**2-(self.radius**2*math.acos(hp/self.radius)-hp*math.sqrt(self.radius**2-hp**2))
        else:
            self.area=math.pi*self.radius**2
            
        
        self.nuSR=self.area*self.nu*self.radius
        self.wl=self.Rel*self.nu/self.radius/self.d
        
        static_matrix=zeros((6,1))
        static_matrix[3,0]=1# Resistant torque on X
        static_behavior_occurence_matrix=array([[1]])
        static_behavior_nonlinear_eq_indices=[0]
        static_behavior_linear_eq=array([])
        static_behavior_nonlinear_eq=[lambda x,w,v:
            x[0]+self.Cl*self.nuSR*w[0]
            if abs(w[0])<self.wl
            else x[0]+self.nuSR*(self.Ct*w[0]+(self.Cl-self.Ct)*self.wl*w[0]/abs(w[0]))]
        static_require_kinematic=True
        
        UnknownLoad.__init__(self,part,position,euler_angles,static_matrix,
                 static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,
                 static_behavior_linear_eq,static_behavior_nonlinear_eq,
                 static_require_kinematic,name)        
        
        
    def ChangeCoefficients(self,Cl,Ct,Rel):
        self.Cl=Cl
        self.Ct=Ct
        self.Rel=Rel# limit reynods number
        self.wl=self.Rel*self.nu/self.radius/self.d
        self.static_behavior_nonlinear_eq=[lambda x,w,v:
            x[0]+self.Cl*self.nuSR*w[0]
            if abs(w[0])<self.wl
            else x[0]+self.nuSR*(self.Ct*w[0]+(self.Cl-self.Ct)*self.wl*w[0]/abs(w[0]))]
