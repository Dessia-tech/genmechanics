# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 20:13:02 2016

@author: steven
"""

import networkx as nx


def construction_sources_and_sinks(graphe):
    sinks = []
    sources = []
    for node in graphe.nodes():
        if graphe.out_degree(node) == 0:
            sinks.append(node)
        elif graphe.in_degree(node) == 0:
            sources.append(node)

    g2 = sources[:]

    for node in sources:
        for node2 in nx.descendants(graphe, node):
            if node2 not in g2:
                g2.append(node2)

    g3 = sinks[:]
    for node in sinks:
        for node2 in nx.ancestors(graphe, node):
            if node2 not in g3:
                g3.append(node2)
    return g2, g3


def construction_graphes_system_analysis(Mo):
    neq, nvar = Mo.shape
    g = nx.Graph()
    gp = nx.DiGraph()
    pos = {}
    for i in range(nvar):
        g.add_node('v' + str(i), bipartite=0)
        gp.add_node('v' + str(i), bipartite=0)
        pos['v' + str(i)] = [i, 0]

    for i in range(neq):
        g.add_node('e' + str(i), bipartite=1)
        gp.add_node('e' + str(i), bipartite=1)
        pos['e' + str(i)] = [i, 1]
        for j in range(nvar):
            if Mo[i, j] == 1:
                g.add_edge('e' + str(i), 'v' + str(j))
                gp.add_edge('e' + str(i), 'v' + str(j))

    for gi in (g.subgraph(c).copy() for c in nx.connected_components(g)):
        M = nx.bipartite.maximum_matching(gi)
        #    print('M',M,len(M))

        for n1, n2 in M.items():
            # print(n1,n2)
            gp.add_edge(n1, n2)

    return g, gp


def calcul_order_ev(graph, strongly_connected_components, solvable_vars, ):
    c = nx.condensation(graph, strongly_connected_components)
    isc_vars = []
    for isc, sc in enumerate(strongly_connected_components):
        # print(sc)
        for var in solvable_vars:
            if 'v' + str(var) in sc:
                isc_vars.append(isc)
                break

    ancetres_vars = isc_vars[:]

    for isc_var in isc_vars:
        for ancetre in nx.ancestors(c, isc_var):
            if ancetre not in ancetres_vars:
                ancetres_vars.append(ancetre)

    ordre_sc = [sc for sc in nx.topological_sort(c) if sc in ancetres_vars]
    ordre_ev = []
    for isc in ordre_sc:
        evs = sorted(strongly_connected_components[isc])  # liste d'équations et de variables triées pour être séparées

        levs = int(len(evs) / 2)
        ordre_ev.append(([int(e[1:]) for e in evs[0:levs]], [int(v[1:]) for v in evs[levs:]]))

    return ordre_ev


def equations_system_analysis(Mo, vars_to_solve, overconstrain_stop=True):
    """
        Analyse a free equations system given by its ocurence matrix.
        :return: False (if system is unsolvable) if overconstrain_stop==True
        else, returns True, the solvable variables and the resolution order
    """

    if vars_to_solve is None:
        vars_to_solve = range(Mo.shape[1])

    g, gp, = construction_graphes_system_analysis(Mo)

    g2, g3 = construction_sources_and_sinks(gp)

    if overconstrain_stop:
        if g2:

            return (False, [], None)

    solvable_vars = []
    for var in vars_to_solve:
        if not 'v' + str(var) in g2 + g3:
            solvable_vars.append(var)

    g1 = g.copy()
    g1.remove_nodes_from(g2 + g3)

    g1p = nx.DiGraph()

    g1p.add_nodes_from(g1.nodes())
    for e in g1.edges():
        # equation vers variable
        if e[0][0] == 'v':
            g1p.add_edge(e[0], e[1])
        else:
            g1p.add_edge(e[1], e[0])

    for g1i in (g1.subgraph(c).copy() for c in nx.connected_components(g1)):  # nx.connected_component_subgraphs(G1):
        M1 = nx.bipartite.maximum_matching(g1i)

        for n1, n2 in M1.items():
            # print(n1,n2)
            if n1[0] == 'e':
                g1p.add_edge(n1, n2)
            else:
                g1p.add_edge(n2, n1)

    scc = list(nx.strongly_connected_components(g1p))

    if scc:
        ordre_ev = calcul_order_ev(g1p, scc, solvable_vars)
        return (True, solvable_vars, ordre_ev)

    return (False, [], None)
