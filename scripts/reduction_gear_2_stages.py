# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:21:07 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages
import numpy as npy

L=0.2
e1=0.05
e2=0.07
C=300
r1=0.8
r2=0.7

Ca=0.0008
Cr=0.0006

ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')
shaft2=genmechanics.Part('shaft2')
shaft3=genmechanics.Part('shaft3')

p1a=npy.array([0,0,0])
p1b=npy.array([L,0,0])
p2a=npy.array([0,e1,0])
p2b=npy.array([L,e1,0])
p3a=npy.array([0,e1+e2,e2])
p3b=npy.array([L,e1+e2,e2])

pgs1=0.5*(p1a+p1b)*r1+(1-r1)*0.5*(p2a+p2b)
pgs2=0.5*(p2a+p2b)*r2+(1-r2)*0.5*(p3a+p3b)

dgs1=npy.cross(p1b-p1a,p2a-p1a)
egs1=genmechanics.geometry.Direction2Euler(*dgs1)
dgs2=npy.cross(p3b-p2a,p3a-p2a)
egs2=genmechanics.geometry.Direction2Euler(*dgs2)

bearing1a=linkages.BallLinkage(ground,shaft1,p1a,[0,0,0],Ca,Cr,'bearing1a')
bearing1b=linkages.LinearAnnularLinkage(ground,shaft1,p1b,[0,0,0],Cr,'bearing1b')
bearing2a=linkages.BallLinkage(ground,shaft2,p2a,[0,0,0],Ca,Cr,'bearing2a')
bearing2b=linkages.LinearAnnularLinkage(ground,shaft2,p2b,[0,0,0],Cr,'bearing2b')
bearing3a=linkages.BallLinkage(ground,shaft3,p3a,[0,0,0],Ca,Cr,'bearing3a')
bearing3b=linkages.LinearAnnularLinkage(ground,shaft3,p3b,[0,0,0],Cr,'bearing3b')


gearset12=linkages.FrictionlessGearSetLinkage(shaft1,shaft2,pgs1,egs1,'Gear set 2')
gearset23=linkages.FrictionlessGearSetLinkage(shaft2,shaft3,pgs2,egs2,'Gear set 2')

mech=genmechanics.Mechanism([bearing1a,bearing1b,bearing2a,bearing2b,bearing3a,bearing3b,gearset12,gearset23],ground)

load1=genmechanics.KnownMechanicalLoad(shaft1,[-L/4,0,0],[0,0,0],[0,0,0],[C,0,0],'input torque')
load2=genmechanics.UnknownMechanicalLoad(shaft3,[L/2,0,0],[0,0,0],[],[0],'output torque')


r=mech.StaticAnalysis([load1],[load2])
#print(r1)
#for i in range(100):
#    r2=mech.StaticAnalysis([load1],[load2])
#    if not npy.allclose(r1,r2):
#        print('!!')
#    else:
#        print('??')
        

for (l,d),v in r.items():
    print(l.name,d,v)

print('Cth: ',r1/(1-r1)*r2/(1-r2)*C)
