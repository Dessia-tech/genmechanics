#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:15:09 2017


"""

from numpy import array,zeros
import genmechanics.geometry as geometry

class KnownMechanicalLoad:
    def __init__(self,part,position,euler_angles,forces,torques,name=''):
        self.part=part
        self.position=position
        self.euler_angles=euler_angles
        self.forces=array(forces)
        self.torques=array(torques)
        self.name=name
        
        self.P=geometry.Euler2TransferMatrix(*self.euler_angles) 
        

class UnknownMechanicalLoad:
    """
    :param force_directions: a list of directions for force (0,1,2)
    :param torque_directions: a list of directions for torque (0,1,2)
    """
    def __init__(self,part,position,euler_angles,force_directions,torque_directions,name=''):
        
        self.part=part
        self.position=position
        self.euler_angles=euler_angles
        self.force_directions=force_directions
        self.torque_directions=torque_directions
#        self.force_matrix=npy.array([[force[0],0,0],[0,force[1],0],[0,0,force[2]]])
#        self.torque_matrix=npy.array([[torque[0],0,0],[0,torque[1],0],[0,0,torque[2]]])
        lfd=len(force_directions)
        ltd=len(torque_directions)
        self.static_matrix=zeros((6,lfd+ltd))
        for i,k in enumerate(force_directions):
            self.static_matrix[k,i]=1
        for i,k in enumerate(torque_directions):
            self.static_matrix[k+3,i+lfd]=1
        self.name=name
        
        self.P=geometry.Euler2TransferMatrix(*self.euler_angles) 
        self.n_static_unknowns=self.static_matrix.shape[1]

#class SplashLoad(UnknownMechanicalLoad):
#    """
#    Creates a splash force linked to the rotationnal speed around
#    """
#    def __init__(self,part,position,euler_angles,force_directions,torque_directions):
#        
#        UnknownMechanicalLoad.__init__(self,part,position,euler_angles,force_directions,torque_directions,name='')
    