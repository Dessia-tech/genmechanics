# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 17:31:34 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads
import numpy as npy
import scipy.linalg as linalg

F=100
L=0.1

ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')

bearing=genmechanics.linkages.RevoluteLinkage(ground,shaft1,[0,0,0],[0,0,0],'bearing')

mech=genmechanics.Mechanism([bearing],ground)

load1=loads.KnownMechanicalLoad(shaft1,[0,L,0],[0,0,0],[0,0,F],[0,0,0],'load')
load2=loads.UnknownMechanicalLoad(shaft1,[0,0,0],[0,0,0],[],[0],'output_torque')

r=mech.StaticAnalysis([load1],[load2])

for (l,d),v in r.items():
    print(l.name,d,v)