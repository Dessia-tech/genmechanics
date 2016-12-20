# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:06:30 2016

@author: steven
"""

import numpy as npy
import networkx as nx
from genmechanics import geometry,tools
from scipy import linalg
from scipy.optimize import fsolve

class Part:
    def __init__(self,name=''):
        self.name=name
        
class Linkage:
    def __init__(self,part1,part2,position,euler_angles,static_matrix,static_behavior_occurence_matrix,static_behavior_nonlinear_eq_indices,static_behavior_linear_eq,static_behavior_nonlinear_eq,name=''):
        self.part1=part1
        self.part2=part2
        self.position=position
        self.euler_angles=euler_angles
        self.name=name        
        
        self.static_matrix=static_matrix
        self.static_behavior_occurence_matrix=static_behavior_occurence_matrix
        
        self.static_behavior_nonlinear_eq_indices=static_behavior_nonlinear_eq_indices
        self.static_behavior_linear_eq=static_behavior_linear_eq
        self.static_behavior_nonlinear_eq=static_behavior_nonlinear_eq
        
        self.n_static_unknowns=static_matrix.shape[1]

        
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
                
        
        lparts=len(self.parts)
        M=npy.zeros((6*lparts,self.n_sdof))
        K=npy.zeros((6*lparts,self.n_sdof))
        F=npy.zeros((6*lparts,1))
        q=npy.zeros((self.n_sdof,1))
        nonlinear_eq={}

        # Occurence matrix assembly
        for ip,part in enumerate(self.parts):
            # linkage contribution to the LHS
            for linkage in self.graph[part].keys():
                P=geometry.Euler2TransferMatrix(*linkage.euler_angles)            
                u=linkage.position
#                uprime=npy.dot(npy.transpose(P),u)
                uprime=u
                L=geometry.CrossProductMatrix(uprime)
                Me1=npy.abs(npy.dot(P,linkage.static_matrix[:3,:]))>1e-10
                Me2=npy.abs(npy.dot(L,npy.dot(P,linkage.static_matrix[:3,:]))+npy.dot(P,linkage.static_matrix[3:,:]))>1e-10
#                print(Me1,Me2)
                Me=npy.vstack([Me1,Me2])
                for indof,ndof in enumerate(self.sdof[linkage]):
                    M[ip*6:(ip+1)*6,ndof]+=Me[:,indof]

                Ke1=npy.dot(P,linkage.static_matrix[:3,:])
                Ke2=npy.dot(L,npy.dot(P,linkage.static_matrix[:3,:]))+npy.dot(P,linkage.static_matrix[3:,:])
                Ke=npy.vstack([Ke1,Ke2])
                for indof,ndof in enumerate(self.sdof[linkage]):
                    K[ip*6:(ip+1)*6,ndof]+=Ke[:,indof]

                    
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
                Me1=npy.abs(npy.dot(P,load.static_matrix[:3,:]))>1e-10
                Me2=npy.abs(npy.dot(L,npy.dot(P,load.static_matrix[:3,:]))+npy.dot(P,load.static_matrix[3:,:]))>1e-10
                Me=npy.vstack([Me1,Me2])
#                print(Me)
                for indof,ndof in enumerate(self.sdof[load]):
                    M[ip*6:(ip+1)*6,ndof]+=Me[:,indof]
                
                Ke1=npy.dot(P,load.static_matrix[:3,:])
                Ke2=npy.dot(L,npy.dot(P,load.static_matrix[:3,:]))+npy.dot(P,load.static_matrix[3:,:])
                Ke=npy.vstack([Ke1,Ke2])
                for indof,ndof in enumerate(self.sdof[load]):
                    K[ip*6:(ip+1)*6,ndof]+=Ke[:,indof]      
                    
                    
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
                F[ip*6:(ip+1)*6,:]+=Fe
        
        neq=6*lparts
        neq_linear=neq
        indices_r=list(range(neq))
        
        # behavior equations of linkages
        for linkage in self.linkages:
            neq_linkage=linkage.static_behavior_occurence_matrix.shape[0]
#            print(neq_linkage)
            if neq_linkage>0:
#                print('iii')
                # Adding to occurence matrix
                Me=npy.zeros((neq_linkage,self.n_sdof))
#                print(linkage.static_behavior_occurence_matrix)
#                print(self.sdof[linkage])
                for indof,ndof in enumerate(self.sdof[linkage]):
                    Me[:,ndof]+=linkage.static_behavior_occurence_matrix[:,indof]    
#                print(Me)
                M=npy.vstack([M,Me])

                # Adding linear equations to System matrix K
                neq_linear_linkage=linkage.static_behavior_linear_eq.shape[0]
#                print(n)
                if neq_linear_linkage>0:
                    Ke=npy.zeros((neq_linear_linkage,self.n_sdof))
                    for indof,ndof in enumerate(self.sdof[linkage]):
                        Ke[neq_linear:neq_linear+6,ndof]+=linkage.static_behavior_linear_eq[:,indof]   
#                    print(Ke)
                    K=npy.vstack([K,Ke])
#                    print(neq_linear_linkage)
                    indices_r.extend(range(neq_linear,neq_linear+neq_linear_linkage))

                # Collecting non linear equations
                for i,eq in zip(linkage.static_behavior_nonlinear_eq_indices,linkage.static_behavior_nonlinear_eq):
                    nonlinear_eq[neq+i]=eq

                        
                # Updating counters
                neq+=neq_linkage
                neq_linear+=neq_linear_linkage
                

        
#        print(K)
#        return M
        
#        for i in range(30):        
        solvable,solvable_var,resolution_order=tools.EquationsSystemAnalysis(M,None)
        print(solvable)
        print(F)
#        print(resolution_order)
#        print(nonlinear_eq)
        for eqs,variables in resolution_order:
            print(eqs,variables)
            linear=True
            linear_eqs=[]
            for eq in eqs:
                try:
                    nonlinear_eq[eq]
                    linear=False
                    print('!!!!!')
                except KeyError:
                    linear_eqs.append(eq)            
            if linear:
                eqs_r=npy.array([indices_r[eq] for eq in eqs])
#                print('eqr',eqs,)
                other_vars=npy.array([i for i in range(self.n_sdof) if i not in variables])
                Kr=K[eqs_r[:,None],npy.array(variables)]
                Fr=F[eqs_r,:]-npy.dot(K[eqs_r[:,None],other_vars],q[other_vars,:])
#                print(Kr,Fr)
                q[variables,0]=linalg.solve(Kr,Fr)
            else:
                nl_eqs=[]
                other_vars=npy.array([i for i in range(self.n_sdof) if i not in variables])
                print('ov: ',other_vars)
                for eq in eqs:
                    try:
#                        ind_var=[variables.index(i) for i in variables if M[eq,i]]
                        f1=nonlinear_eq[eq]
                        vars_func=[i for i in range(self.n_sdof) if M[eq,i]]
                        def f2(x,f1=f1,vars_func=vars_func,variables=variables,q=q):
                            
                            x2=[]
                            for variable in vars_func:
                                try:
                                    x2.append(x[variables.index(variable)])
                                except ValueError:
                                    x2.append(q[variable,0])
#                            print('xx2',x,x2)
                            return f1(x2)
                        print(f2)
                        nl_eqs.append(f2)
                    except KeyError:
                        # lambdification of linear equations
                        print('Fnl: ',F[indices_r[eq],0])
                        f2=lambda x,indices_r=indices_r,K=K,F=F,eq=eq,q=q,other_vars=other_vars:(npy.dot(K[indices_r[eq],variables],x)-F[indices_r[eq],0]+npy.dot(K[indices_r[eq],other_vars],q[other_vars,:]))
                        nl_eqs.append(f2)
                        print('f2nl: ',f2(npy.zeros(len(variables))),f2)
#                print(nl_eqs)
                print('fj')
#                import inspect
                for fj in nl_eqs:
                    print(fj,fj(npy.zeros(len(variables))))
                f=lambda x:[fi(x) for fi in nl_eqs]
#                print(f(npy.random.random(len(variables))))
                xs=fsolve(f,npy.zeros(len(variables)))
#                print('xfx',xs,f(xs))
                q[variables,0]=xs
#                print(q[variables,0])
#            print(q[variables,0])
#            print(q)

#        return f
#        if solvable:
#            print(solvable_var)
#            print(resolution_order)

#        # Static assembly
#        lparts=len(self.parts)
#        K=npy.zeros((6*lparts,self.n_sdof))
#        F=npy.zeros((6*lparts,1))
#        for ip,part in enumerate(self.parts):
#            # linkage contribution to the LHS
#            for linkage in self.graph[part].keys():
#                P=geometry.Euler2TransferMatrix(*linkage.euler_angles)            
#                u=linkage.position
##                uprime=npy.dot(npy.transpose(P),u)
#                uprime=u
#                L=geometry.CrossProductMatrix(uprime)
#                J1=npy.dot(P,linkage.static_matrix[:3,:])
#                J2=npy.dot(L,npy.dot(P,linkage.static_matrix[:3,:]))+npy.dot(P,linkage.static_matrix[3:,:])
#                J=npy.vstack([J1,J2])
##                print(J)
#                if linkage.part1 is part:
#                    side=1
#                else:
#                    side=-1
#                for indof,ndof in enumerate(self.sdof[linkage]):
#                    K[ip*6:(ip+1)*6,ndof]+=side*J[:,indof]
#
#            # Unknowns loads contribution to the LHS
#            try:
#                uloads=uloads_parts[part]
#            except:
#                uloads=[]                
#            for load in uloads:
#                P=geometry.Euler2TransferMatrix(*load.euler_angles)            
#                u=load.position
##                uprime=npy.dot(npy.transpose(P),u)
#                uprime=u
#                L=geometry.CrossProductMatrix(uprime)
#                J1=npy.dot(P,load.static_matrix[:3,:])
#                J2=npy.dot(L,npy.dot(P,load.static_matrix[:3,:]))+npy.dot(P,load.static_matrix[3:,:])
#                J=npy.vstack([J1,J2])
#                for indof,ndof in enumerate(self.sdof[load]):
#                    K[ip*6:(ip+1)*6,ndof]+=J[:,indof]
#
#            # knowns loads contribution to the RHS
#            try:
#                loads=loads_parts[part]
#            except:
#                loads=[]                
#            for load in loads:
#                P=geometry.Euler2TransferMatrix(*load.euler_angles)            
#                u=load.position
##                uprime=npy.dot(npy.transpose(P),u)
#                uprime=u
#                L=geometry.CrossProductMatrix(uprime)
#                F1=npy.dot(P,load.force)
#                F2=npy.dot(L,npy.dot(P,load.force))+npy.dot(P,load.torque)
#                Fe=npy.vstack([F1.reshape((3,1)),F2.reshape((3,1))])
##                print(F)
##                print(F1.shape,F2.shape)
##                print(Fe)
##                for indof,ndof in enumerate(self.sdof[load]):
#                F[ip*6:(ip+1)*6,:]+=Fe
##        F=F[~npy.all(K==0,axis=1)]
##        K=K[~npy.all(K==0,axis=1)]
##        print(K.shape,F.shape)
##        import matplotlib.pyplot as plt
##        plt.figure()
##        plt.pcolor(K)
##        plt.gca().invert_yaxis()
##        plt.colorbar()
##        print(K)
#        q=linalg.solve(K,F)
##        q=Kplus=linalg.pinv(K)
##        q=npy.dot(Kplus,F)
        results={}
        for link,dofs in self.sdof.items():
            for idof,dof in enumerate(dofs):
                results[link,idof]=q[dof,0]
        return results
        
        