#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 10:50:15 2016

@author: steven
"""


import genmechanics
import genmechanics.linkages as linkages
import numpy as npy

L=0.2
e1=0.05
C=300
r1=0.8
W=300

Ca=0.0008
Cr=0.0006

ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')
shaft2=genmechanics.Part('shaft2')

p1a=npy.array([0,0,0])
p1b=npy.array([L,0,0])
p2a=npy.array([0,e1,0])
p2b=npy.array([L,e1,0])

pgs1=0.5*(p1a+p1b)*r1+(1-r1)*0.5*(p2a+p2b)

dgs1=npy.cross(p1b-p1a,p2a-p1a)
egs1=genmechanics.geometry.Direction2Euler(*dgs1)

bearing1=linkages.FrictionlessRevoluteLinkage(ground,shaft1,p1a,[0,0,0],'bearing1a')
#bearing1=linkages.LinearAnnularLinkage(ground,shaft1,p1b,[0,0,0],Cr,'bearing1b')
bearing2=linkages.FrictionlessRevoluteLinkage(ground,shaft2,p2a,[0,0,0],'bearing2a')
#bearing2b=linkages.LinearAnnularLinkage(ground,shaft2,p2b,[0,0,0],Cr,'bearing2b')

gearset12=linkages.FrictionlessGearSetLinkage(shaft1,shaft2,pgs1,egs1,'Gear set 2')

load1=genmechanics.KnownMechanicalLoad(shaft1,[-L/4,0,0],[0,0,0],[0,0,0],[C,0,0],'input torque')
load2=genmechanics.UnknownMechanicalLoad(shaft2,[L/2,0,0],[0,0,0],[],[0],'output torque')
imposed_speeds=[(bearing1,0,W)]

mech=genmechanics.Mechanism([bearing1,bearing2,gearset12],ground,imposed_speeds,[load1],[load2])


for l,lv in mech.static_results.items():
    for d,v in lv.items():
        print(l.name,d,v)
print('=================')
for l,lv in mech.kinematic_results.items():
    for d,v in lv.items():
        print(l.name,d,v)
        
        
print('wth: ',(1-r1)/r1*W)
