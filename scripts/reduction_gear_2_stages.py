# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:21:07 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads

import numpy as npy
import scipy.linalg as linalg

Lb=0.2# bearing width
Lgs=0.045# Gearset width


C=300
w=300
r1=0.048
r2=0.047
y2=0.078
z2=-0.026
y3=0.203
z3=0


Ca=0.0008
Cr=0.0006
Cf=0.01
Cwb=0.00001# Speed coeff for bearings
Cvgs=0.0001# Speed coeff for gear sets

alpha_gs1=20/360*2*3.1415
beta_gs1=25/360*2*3.1415
alpha_gs2=20/360*2*3.1415
beta_gs2=-25/360*2*3.1415


ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')
shaft2=genmechanics.Part('shaft2')
shaft3=genmechanics.Part('shaft3')

p1a=npy.array([0.,0.,0.])
p1b=npy.array([Lb+Lgs,0.,0.])
p2a=npy.array([0,y2,z2])
p2b=npy.array([Lb+2*Lgs,y2,z2])
p3a=npy.array([0,y3,z3])
p3b=npy.array([Lb+2*Lgs,y3,z3])

pgs11=npy.array([0.5*(Lgs+Lb),0,0])
pgs12=npy.array([0.5*(Lgs+Lb),y2,z2])
pgs22=npy.array([0.5*(3*Lgs+Lb),y2,z2])
pgs23=npy.array([0.5*(3*Lgs+Lb),y3,z3])

pgs1=pgs11+(pgs12-pgs11)/linalg.norm(pgs12-pgs11)*r1
pgs2=pgs22+(pgs23-pgs22)/linalg.norm(pgs23-pgs22)*r2

dir_axis=npy.array([1,0,0])

dgs1=npy.cross(p1b-p1a,p2a-p1a)
egs1=genmechanics.geometry.Direction2Euler(dgs1,dir_axis)
dgs2=npy.cross(p3b-p2a,p3a-p2a)
egs2=genmechanics.geometry.Direction2Euler(dgs2,dir_axis)

bearing1a=linkages.BallLinkage(ground,shaft1,p1a,[0,0,0],Ca,Cr,Cwb,'bearing1a')
bearing1b=linkages.LinearAnnularLinkage(ground,shaft1,p1b,[0,0,0],Cr,Cwb,'bearing1b')
bearing2a=linkages.BallLinkage(ground,shaft2,p2a,[0,0,0],Ca,Cr,Cwb,'bearing2a')
bearing2b=linkages.LinearAnnularLinkage(ground,shaft2,p2b,[0,0,0],Cr,Cwb,'bearing2b')
bearing3a=linkages.BallLinkage(ground,shaft3,p3a,[0,0,0],Ca,Cr,Cwb,'bearing3a')
bearing3b=linkages.LinearAnnularLinkage(ground,shaft3,p3b,[0,0,0],Cr,Cwb,'bearing3b')


gearset12=linkages.GearSetLinkage(shaft1,shaft2,pgs1,egs1,alpha_gs1,beta_gs1,Cf,Cvgs,'Gear set 1')
gearset23=linkages.GearSetLinkage(shaft2,shaft3,pgs2,egs2,alpha_gs2,beta_gs2,Cf,Cvgs,'Gear set 2')

load1=loads.KnownLoad(shaft1,[-Lb/2,0,0],[0,0,0],[0,0,0],[C,0,0],'input torque')
load2=loads.SimpleUnknownLoad(shaft3,[2*(Lgs+Lb),y2,z2],[0,0,0],[],[0],'output torque')

imposed_speeds=[(bearing1a,0,w)]

mech=genmechanics.Mechanism([bearing1a,bearing1b,bearing2a,bearing2b,bearing3a,bearing3b,gearset12,gearset23],ground,imposed_speeds,[load1],[load2])


for l,lv in mech.static_results.items():
    for d,v in lv.items():
        print(l.name,d,v)

#print('Cth: ',r1/(1-r1)*r2/(1-r2)*C)# False now

for l,lv in mech.kinematic_results.items():
    for d,v in lv.items():
        print(l.name,d,v)

#print('wth: ',(1-r1)/r1*(1-r2)/r2*w) # False now

mech.GlobalSankey()

