# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:21:07 2016

@author: steven
"""

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads
import numpy as npy
import matplotlib.pyplot as plt
import scipy.linalg as linalg

# plt.figure()
L=0.2
e1=0.05
e2=0.07
C=300
w=300
r1=0.3

Ca=0.000
Cr=0.000
Cf=0.01
Cwb=0.# Speed coeff for bearings
Cvgs=0.0001# Speed coeff for gear sets

alpha_gs1=18/360*2*3.1415
beta_gs1=25/360*2*3.1415

ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')
shaft2=genmechanics.Part('shaft2')

p1a=npy.array([0,0,0])
p1b=npy.array([L,0,0])
p2a=npy.array([0,0,e1])
p2b=npy.array([L,0,e1])


pgs1 = 0.5*(p1a+p1b)*r1+(1-r1)*0.5*(p2a+p2b)
print(pgs1)

dir_axis=npy.array([0,0,1])
dgs1=npy.cross(p1b-p1a, p2a-p1a)
dgs1 = [0, 1, 0]
egs1=[0,0,0]
bearing1a=linkages.BallLinkage(ground,shaft1,p1a,egs1,Ca,Cr,Cwb,'bearing 1a')
bearing1b=linkages.LinearAnnularLinkage(ground,shaft1,p1b,egs1,Cr,Cwb,'bearing 1b')
bearing2a=linkages.BallLinkage(ground,shaft2,p2a,egs1,Ca,Cr,Cwb,'bearing 2a')
bearing2b=linkages.LinearAnnularLinkage(ground,shaft2,p2b,egs1,Cr,Cwb,'bearing 2b')

dir_axis=npy.array([1,0,0])
dgs1=npy.cross(p1b-p1a, p2a-p1a)
dgs1 = [0, 1, 0]
egs1=genmechanics.geometry.Direction2Euler(dgs1,[0,0,1])
gearset12=linkages.GearSetLinkage(shaft1, shaft2, pgs1, radial_vector=p2a-p1a, axial_vector=dir_axis,
                                  pressure_angle=alpha_gs1, helix_angle=beta_gs1,Cf=Cf, Cv=Cvgs, name='Gear set 1')
#gearset23=linkages.GearSetLinkage(shaft2,shaft3,[3*L/2,e1*r1+(e1+e2)*(1-r1),0],dir23)

imposed_speeds=[(bearing1a,0,w)]

load1=loads.KnownLoad(shaft1,[-L/4,0,0],[0,0,0],[0,0,0],[C,0,0],'input torque')
load2=loads.SimpleUnknownLoad(shaft2,[3*L/2,0,0],[0,0,0],[],[0],'output torque')

mech=genmechanics.Mechanism([bearing1a, bearing1b, bearing2a, bearing2b, gearset12],
                            ground, imposed_speeds, [load1],[load2])

#r=mech.StaticAnalysis()

for l,r in mech.static_results.items():
    for d,v in r.items(): 
        print(l.name,d,v)
    
    
print('Cth: ',-r1/(1-r1)*C)


mech.GlobalSankey()