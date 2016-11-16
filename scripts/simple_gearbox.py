# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:21:07 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages

L=0.2
F=2000

ground=genmechanics.Part('ground')
shaft=genmechanics.Part('shaft')

bearing1=linkages.RevoluteLinkage(ground,shaft,[L,0,0],[0,0,0])
#bearing2=linkages.RevoluteLinkage(ground,shaft,[L,0,0],[0,0,0])
mech=genmechanics.Mechanism([bearing1],ground)

load1=genmechanics.KnownMechanicalLoad(shaft,[L/2,0,0],[0,0,0],[F,0,0],[300,0,0])
load2=genmechanics.UnknownMechanicalLoad(shaft,[L/2,0,0],[0,0,0],[0,0,0],[1,0,0])

r=mech.StaticSolve([load1],[load2])

print(r)