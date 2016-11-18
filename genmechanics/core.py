# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:06:30 2016

@author: steven
"""

import numpy as npy
import networkx as nx
from genmechanics import geometry
from scipy import linalg

class Part:
    def __init__(self,name=''):
        self.name=name
        
class Linkage:
    def __init__(self,part1,part2,position,euler_angles,kinematic_unknowns,static_unknowns,name=''):
        self.part1=part1
        self.part2=part2
        self.position=position
        self.euler_angles=euler_angles
        self.name=name        
        
        self.n_kinematic_unknowns=len(kinematic_unknowns)
        self.n_static_unknowns=len(static_unknowns)
        self.kinematic_matrix=npy.zeros((6,len(kinematic_unknowns)))
        for i,k in enumerate(kinematic_unknowns):
            self.kinematic_matrix[k,i]=1
        self.static_matrix=npy.zeros((6,len(static_unknowns)))
        for i,k in enumerate(static_unknowns):
            self.static_matrix[k,i]=1

        
class KnownMechanicalLoad:
    def __init__(self,part,position,euler_angles,force,torque,name=''):
        self.part=part
        self.position=position
        self.euler_angles=euler_angles
        self.force=force
        self.torque=torque
        self.force_matrix=npy.array([[force[0],0,0],[0,force[1],0],[0,0,force[2]]])
        self.torque_matrix=npy.array([[torque[0],0,0],[0,torque[1],0],[0,0,torque[2]]])
        self.name=name

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
        self.static_matrix=npy.zeros((6,lfd+ltd))
        for i,k in enumerate(force_directions):
            self.static_matrix[k,i]=1
        for i,k in enumerate(torque_directions):
            self.static_matrix[k+3,i+lfd]=1
        self.name=name
        
class Mechanism:
    def __init__(self,linkages,ground,name=''):
        self.linkages=linkages
        self.ground=ground
        self.name=name

    def _get_parts(self):
        parts=[]
        for linkage in self.linkages:
            for part in [linkage.part1,linkage.part2]:
                if not part in parts:
                    if part!=self.ground:
                        parts.append(part)
        return parts

    parts=property(_get_parts)

    def _get_graph(self):
        G=nx.Graph()
        G.add_nodes_from(self.parts)
        for linkage in self.linkages:
            G.add_node(linkage)
            G.add_edge(linkage,linkage.part1)
            G.add_edge(linkage,linkage.part2)
        return G

    graph=property(_get_graph)

        
    def StaticAnalysis(self,known_loads,unknown_loads):
        # Static parametrisation
        self.sdof={}
        self.n_sdof=0
        for linkage in self.linkages:
#            for i in range(linkage.n_static_unknows):
            self.sdof[linkage]=list(range(self.n_sdof,self.n_sdof+linkage.n_static_unknowns))
            self.n_sdof+=linkage.n_static_unknowns
        for load in unknown_loads:
            load_unknowns=len(load.force_directions)+len(load.torque_directions)
            self.sdof[load]=list(range(self.n_sdof,self.n_sdof+load_unknowns))
            self.n_sdof+=load_unknowns
            
        # Loading sorting by part
        uloads_parts={}
        for load in unknown_loads:
            try:
                uloads_parts[load.part].append(load)
            except:
                uloads_parts[load.part]=[load]

        loads_parts={}
        for load in known_loads:
            try:
                loads_parts[load.part].append(load)
            except:
                loads_parts[load.part]=[load]
                
#        print(loads_parts,uloads_parts)
                
        
        # Static assembly
        lparts=len(self.parts)
        K=npy.zeros((6*lparts,self.n_sdof))
        F=npy.zeros((6*lparts,1))
        for ip,part in enumerate(self.parts):
            # linkage contribution to the LHS
            for linkage in self.graph[part].keys():
                P=geometry.Euler2TransferMatrix(*linkage.euler_angles)            
                u=linkage.position
#                uprime=npy.dot(npy.transpose(P),u)
                uprime=u
                L=geometry.CrossProductMatrix(uprime)
                J1=npy.dot(P,linkage.static_matrix[:3,:])
                J2=npy.dot(L,npy.dot(P,linkage.static_matrix[:3,:]))+npy.dot(P,linkage.static_matrix[3:,:])
                J=npy.vstack([J1,J2])
#                print(J)
                if linkage.part1 is part:
                    side=1
                else:
                    side=-1
                for indof,ndof in enumerate(self.sdof[linkage]):
                    K[ip*6:(ip+1)*6,ndof]+=side*J[:,indof]

            # Unknowns loads contribution to the LHS
            try:
                uloads=uloads_parts[part]
            except:
                uloads=[]                
            for load in uloads:
                P=geometry.Euler2TransferMatrix(*load.euler_angles)            
                u=load.position
#                uprime=npy.dot(npy.transpose(P),u)
                uprime=u
                L=geometry.CrossProductMatrix(uprime)
                J1=npy.dot(P,load.static_matrix[:3,:])
                J2=npy.dot(L,npy.dot(P,load.static_matrix[:3,:]))+npy.dot(P,load.static_matrix[3:,:])
                J=npy.vstack([J1,J2])
                for indof,ndof in enumerate(self.sdof[load]):
                    K[ip*6:(ip+1)*6,ndof]+=J[:,indof]

            # knowns loads contribution to the RHS
            try:
                loads=loads_parts[part]
            except:
                loads=[]                
            for load in loads:
                P=geometry.Euler2TransferMatrix(*load.euler_angles)            
                u=load.position
#                uprime=npy.dot(npy.transpose(P),u)
                uprime=u
                L=geometry.CrossProductMatrix(uprime)
                F1=npy.dot(P,load.force)
                F2=npy.dot(L,npy.dot(P,load.force))+npy.dot(P,load.torque)
                Fe=npy.vstack([F1.reshape((3,1)),F2.reshape((3,1))])
#                print(F)
#                print(F1.shape,F2.shape)
#                print(Fe)
#                for indof,ndof in enumerate(self.sdof[load]):
                F[ip*6:(ip+1)*6,:]+=Fe
#        F=F[~npy.all(K==0,axis=1)]
#        K=K[~npy.all(K==0,axis=1)]
#        print(K.shape,F.shape)
#        import matplotlib.pyplot as plt
#        plt.figure()
#        plt.pcolor(K)
#        plt.gca().invert_yaxis()
#        plt.colorbar()
#        print(K)
        q=linalg.solve(K,F)
#        q=Kplus=linalg.pinv(K)
#        q=npy.dot(Kplus,F)
        results={}
        for link,dofs in self.sdof.items():
            for idof,dof in enumerate(dofs):
                results[link,idof]=q[dof,0]
        return results
        
        