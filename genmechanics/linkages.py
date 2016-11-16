# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:14:34 2016

@author: steven
"""

from genmechanics import Linkage
import numpy as npy

class RevoluteLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,name)
        
        self.n_static_unknows=5
        self.static_unknows=npy.identity(6)
        self.static_unknows[3,3]=0
        self.n_kinematic_unknows=1
        self.kinematic_unknows=npy.identity(1)
        self.static_matrix=npy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])

class CylindricalLinkage(Linkage):
    def __init__(self,position,euler_angles,name=''):
        Linkage.__init__(self,position,euler_angles,name)

class PrismaticLinkage(Linkage):
    def __init__(self,position,euler_angles,name=''):
        Linkage.__init__(self,position,euler_angles,name)
