# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 20:13:02 2016

@author: steven
"""

import networkx as nx

def EquationsSystemAnalysis(Mo,vars_to_solve,overconstrain_stop=True):
    """
        Analyse a free equations system given by its ocurence matrix.
        :return: False (if system is unsolvable) if overconstrain_stop==True
        else, returns True, the solvable variables and the resolution order
    """
    
    if vars_to_solve==None:
        vars_to_solve=range(Mo.shape[1])
    
    
    neq,nvar=Mo.shape
    G=nx.Graph()
    Gp=nx.DiGraph()
    pos={}
    for i in range(nvar):
        G.add_node('v'+str(i),bipartite=0)
        Gp.add_node('v'+str(i),bipartite=0)
        pos['v'+str(i)]=[i,0]
    
    for i in range(neq):
        G.add_node('e'+str(i),bipartite=1)
        Gp.add_node('e'+str(i),bipartite=1)
        pos['e'+str(i)]=[i,1]
        for j in range(nvar):
            if Mo[i,j]==1:
                G.add_edge('e'+str(i),'v'+str(j))
                Gp.add_edge('e'+str(i),'v'+str(j))
            
        
    
#    plt.figure()
#    pos=nx.spring_layout(G)
#    nx.draw(G,pos)
#    nx.draw_networkx_labels(G,pos)
               
    for Gi in (G.subgraph(c).copy() for c in nx.connected_components(G)):
        M=nx.bipartite.maximum_matching(Gi)
    #    print('M',M,len(M))
        
        for n1,n2 in M.items():
    #        print(n1,n2)
            Gp.add_edge(n1,n2)    
    
#    pos=nx.spring_layout(Gp)
#    plt.figure()
#    nx.draw(Gp,pos)
#    nx.draw_networkx_labels(Gp,pos)
    
    sinks=[]
    sources=[]
    for node in Gp.nodes():
        if Gp.out_degree(node)==0:
            sinks.append(node)
        elif Gp.in_degree(node)==0:
            sources.append(node)
    
    G2=sources[:]
    for node in sources:
        for node2 in nx.descendants(Gp,node):
            if node2 not in G2:
                G2.append(node2)
    
    if overconstrain_stop:
        if G2!=[]:
#            print('G2 (sur-contraint): ',G2)
#            eG2=[int(elem[1:]) for elem in G2 if elem[0]=='e']
#            vG2=[int(elem[1:]) for elem in G2 if elem[0]=='v']
#            print(Mo[eG2,vG2])
            return (False,[],None)

    
    G3=sinks[:]
    for node in sinks:
        for node2 in nx.ancestors(Gp,node):
            if node2 not in G3:
                G3.append(node2)

#    print(sources)
#    print(sinks)
    
#    print('G3 (sous-contraint): ',G3)

    solvable_vars=[]
    for var in vars_to_solve: 
        if not 'v'+str(var) in G2+G3:
            solvable_vars.append(var)
            
        
    G1=G.copy()
    G1.remove_nodes_from(G2+G3)
#    print('G1: ',G1.nodes())

#    pos=nx.spring_layout(G1)
#    plt.figure()
#    nx.draw(G1,pos)
#    nx.draw_networkx_labels(G1,pos)
    
#    M1=nx.bipartite.maximum_matching(G1)
    G1p=nx.DiGraph()

    G1p.add_nodes_from(G1.nodes())
    for e in G1.edges():
        # equation vers variable
        if e[0][0]=='v':
            G1p.add_edge(e[0],e[1])        
        else:
            G1p.add_edge(e[1],e[0])   

    
    for G1i in (G1.subgraph(c).copy() for c in nx.connected_components(G1)):#nx.connected_component_subgraphs(G1):
        M1=nx.bipartite.maximum_matching(G1i)
    

        for n1,n2 in M1.items():
    #        print(n1,n2)
            if n1[0]=='e':
                G1p.add_edge(n1,n2)        
            else:
                G1p.add_edge(n2,n1)
            
        
    scc=list(nx.strongly_connected_components(G1p))
#    pos=nx.spring_layout(G1p)
#    plt.figure()
#    nx.draw(G1p,pos)
#    nx.draw_networkx_labels(G1p,pos)
    
    if scc!=[]:
        C=nx.condensation(G1p,scc)
#        print('scc: ',scc)
#        plt.figure()
#        pos=nx.spectral_layout(C)
#        nx.draw(C,pos)
#        nx.draw_networkx_labels(C,pos)
        # on cherche les indices des blocs contenant chacune des variables
#        print(list(C.nodes()))
        isc_vars=[]
        for isc,sc in enumerate(scc):
#            print(sc)
            for var in solvable_vars:
                if 'v'+str(var) in sc:
                    isc_vars.append(isc)
                    break
#                print(isc)
        ancetres_vars=isc_vars[:]
        
        for isc_var in isc_vars:                                
            for ancetre in nx.ancestors(C,isc_var):
                if ancetre not in ancetres_vars:
                    ancetres_vars.append(ancetre)
                
#        ancetres_var.add(isc_var)
#        print('ancetres: ',ancetres_var)
#        ordre_sc=nx.topological_sort_recursive(C,reverse=False)
#        print(ancetres_vars)
        ordre_sc=[sc for sc in nx.topological_sort(C) if sc in ancetres_vars]
#        print(ordre_sc)
        ordre_ev=[]
        for isc in ordre_sc:            
            evs=sorted(scc[isc])# liste d'équations et de variables triées pour être séparées

            levs=int(len(evs)/2)
#            print(evs,levs,levs/2)
            ordre_ev.append(([int(e[1:]) for e in evs[0:levs]],[int(v[1:]) for v in evs[levs:]]))
            
#        print(ordre_ev)        
        return (True,solvable_vars,ordre_ev)
        
    return (False,[],None)