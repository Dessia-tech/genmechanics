#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 09:39:37 2016

@author: steven
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 10:50:15 2016

@author: steven
"""


import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads
import numpy as npy

L=2.2
C=300
Fa=5000
Fr=2500
W=300

Ca=0.0008
Cr=0.0006

ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')

p1=npy.array([0,0,0])
p2=npy.array([L,0,0])

bearing1=linkages.BallLinkage(ground,shaft1,p1,[0,0,0],Cr,Ca,'bearing1a')
#bearing1=linkages.LinearAnnularLinkage(ground,shaft1,p1b,[0,0,0],Cr,'bearing1b')
bearing2=linkages.LinearAnnularLinkage(ground,shaft1,p2,[0,0,0],Cr,'bearing2a')
#bearing2b=linkages.LinearAnnularLinkage(ground,shaft2,p2b,[0,0,0],Cr,'bearing2b')


load1=loads.KnownMechanicalLoad(shaft1,[-L/8,0,0],[0,0,0],[Fa,Fr,0],[C,0,0],'input torque and loads')
load2=loads.UnknownMechanicalLoad(shaft1,[9*L/8,0,0],[0,0,0],[],[0],'output torque')
imposed_speeds=[(bearing1,0,W)]

mech=genmechanics.Mechanism([bearing1,bearing2],ground,imposed_speeds,[load1],[load2])

#K=mech.static_results['K']
for l,lv in mech.static_results.items():
    for d,v in lv.items():
        print(l.name,d,v)
print('=================')
for l,lv in mech.kinematic_results.items():
    for d,v in lv.items():
        print(l.name,d,v)
        
        