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
    def __init__(self,part1,part2,position,euler_angles,name=''):
        self.part1=part1
        self.part2=part2
        self.position=position
        self.euler_angles=euler_angles
        self.name=name        

        
class KnownMechanicalLoad:
    def __init__(self,part,position,euler_angles,force,torque):
        self.part=part
        self.position=position
        self.euler_angles=euler_angles
        self.force=force
        self.torque=torque
        self.force_matrix=npy.array([[force[0],0,0],[0,force[1],0],[0,0,force[2]]])
        self.torque_matrix=npy.array([[torque[0],0,0],[0,torque[1],0],[0,0,torque[2]]])

class UnknownMechanicalLoad:
    def __init__(self,part,position,euler_angles,force,torque):
        self.part=part
        self.position=position
        self.euler_angles=euler_angles
        self.force=force
        self.torque=torque
        self.force_matrix=npy.array([[force[0],0,0],[0,force[1],0],[0,0,force[2]]])
        self.torque_matrix=npy.array([[torque[0],0,0],[0,torque[1],0],[0,0,torque[2]]])
        
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

        
    def StaticSolve(self,known_loads,unknown_loads):
        # Static parametrisation
        self.sdof={}
        self.n_sdof=0
        for linkage in self.linkages:
#            for i in range(linkage.n_static_unknows):
            self.sdof[linkage]=list(range(self.n_sdof,self.n_sdof+linkage.n_static_unknows))
            self.n_sdof+=linkage.n_static_unknows
        for load in unknown_loads:
            load_unknowns=npy.count_nonzero(load.force)+npy.count_nonzero(load.torque)
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
            # linkage contribution
            for linkage in self.graph[part].keys():
                P=geometry.Euler2TransferMatrix(*linkage.euler_angles)            
                u=linkage.position
                uprime=npy.dot(npy.transpose(P),u)
                L=geometry.CrossProductMatrix(uprime)
#                print(L)
                J1=npy.dot(P,linkage.static_matrix[:3,:])
                J2=npy.dot(L,npy.dot(P,linkage.static_matrix[:3,:]))+npy.dot(P,linkage.static_matrix[3:,:])
                J=npy.vstack([J1,J2])
                if linkage.part1 is part:
                    side=1
                else:
                    side=-1
                for indof,ndof in enumerate(self.sdof[linkage]):
                    K[ip*6:(ip+1)*6,ndof]+=side*J[:,indof]

            try:
                uloads=uloads_parts[part]
            except:
                uloads=[]                
            for load in uloads:
                P=geometry.Euler2TransferMatrix(*load.euler_angles)            
                u=load.position
                uprime=npy.dot(npy.transpose(P),u)
                L=geometry.CrossProductMatrix(uprime)
                J1=npy.dot(P,load.force_matrix)
                J2=npy.dot(L,npy.dot(P,load.force_matrix))+npy.dot(P,load.torque_matrix)
                J=npy.vstack([J1,J2])
                for indof,ndof in enumerate(self.sdof[load]):
                    K[ip*6:(ip+1)*6,ndof]+=J[:,indof]

            try:
                loads=loads_parts[part]
            except:
                loads=[]                
            for load in loads:
                P=geometry.Euler2TransferMatrix(*load.euler_angles)            
                u=load.position
                uprime=npy.dot(npy.transpose(P),u)
                L=geometry.CrossProductMatrix(uprime)
                F1=npy.dot(P,load.force)
                F2=npy.dot(L,npy.dot(P,load.force))+npy.dot(P,load.torque)
                Fe=npy.vstack([F1.reshape((3,1)),F2.reshape((3,1))])
#                print(F)
#                print(F1.shape,F2.shape)
#                print(Fe)
#                for indof,ndof in enumerate(self.sdof[load]):
                F[ip*6:(ip+1)*6,:]+=Fe
                
        q=linalg.solve(K,F)
        results={}
        for link,dofs in self.sdof.items():
            for idof,dof in enumerate(dofs):
                results[link,idof]=q[dof,0]
        return results
        
        