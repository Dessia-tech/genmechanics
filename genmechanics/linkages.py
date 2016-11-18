# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:14:34 2016

@author: steven
"""

from genmechanics import Linkage
import numpy as npy

class RevoluteLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,[0],[0,1,2,4,5],name)

class CylindricalLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,[0,3],[1,2,4,5],name)

class PrismaticLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,[3],[1,2,3,4,5],name)

class BallLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,[0,1,2],[0,1,2],name)

class LinearAnnularLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,[0,1,2,3],[1,2],name)

class GearSetLinkage(Linkage):
    def __init__(self,part1,part2,position,euler_angles,name=''):
        Linkage.__init__(self,part1,part2,position,euler_angles,[],[0],name)
