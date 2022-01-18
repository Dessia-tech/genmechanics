# -*- coding: utf-8 -*-
"""

"""

from itertools import combinations
import numpy as npy
import networkx as nx
from genmechanics import geometry, tools
from scipy import linalg
from scipy.optimize import fsolve

from dessia_common import DessiaObject

import volmdlr as vm
import volmdlr.edges as edges

from genmechanics.templates import babylon_template

import webbrowser
import os


class ModelError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return 'Model Error: ' + self.message


class Part(DessiaObject):
    _eq_is_data_eq = False

    def __init__(self, name='', interest_points=None):
        if interest_points is None:
            self.interest_points = []
        else:
            self.interest_points = interest_points

        self.name = name

    @classmethod
    def wireframe_lines(cls, points):
        if len(points) == 2:
            return [edges.LineSegment3D(points[0], points[1])]
        else:
            lines = []
            full_graph = nx.Graph()
            full_graph.add_nodes_from(points)
            points.append(vm.Point3D.mean_point(points))
            for point1, point2 in combinations(points, 2):
                full_graph.add_edge(point1, point2, weight=point1.point_distance(point2))

            wireframe_graph = nx.minimum_spanning_tree(full_graph)
            for point1, point2 in wireframe_graph.edges():
                lines.append(vm.LineSegment3D(point1, point2))

            return lines


class Mechanism:
    _eq_is_data_eq = False

    def __init__(self,
                 linkages,
                 ground,
                 imposed_speeds,
                 known_static_loads,
                 unknown_static_loads,
                 name=''):

        self.linkages = linkages
        self.ground = ground
        self.name = name
        self.imposed_speeds = imposed_speeds
        self.known_static_loads = known_static_loads
        self.unknown_static_loads = unknown_static_loads

        self._utd_kinematic_results = False
        self._utd_static_results = False

        self._holonomic_paths = {}

    def _get_parts(self):
        parts = []
        for linkage in self.linkages:
            for part in [linkage.part1, linkage.part2]:
                if part not in parts:
                    if part != self.ground:
                        parts.append(part)
        return parts

    parts = property(_get_parts)

    def _get_graph(self):
        g = nx.Graph()
        g.add_nodes_from(self.parts)
        for linkage in self.linkages:
            g.add_node(linkage)
            g.add_edge(linkage, linkage.part1)
            g.add_edge(linkage, linkage.part2)
        return g

    graph = property(_get_graph)

    def _get_holonomic_graph(self):
        g = nx.Graph()
        g.add_nodes_from(self.parts)
        for linkage in self.linkages:
            if linkage.holonomic:
                g.add_node(linkage)
                g.add_edge(linkage, linkage.part1)
                g.add_edge(linkage, linkage.part2)
        return g

    holonomic_graph = property(_get_holonomic_graph)

    def _settings_path(self, part1, part2):
        if (part1, part2) in self._settings_paths:
            return self._settings_paths[part1, part2]
        elif (part2, part1) in self._settings_paths:
            path = [(p2, linkage, not linkage_side, p1)
                    for (p1, linkage, linkage_side, p2)
                    in self._settings_paths[part2, part1][::-1]]
            self._settings_paths[part1, part2] = path
            return self._settings_paths[part1, part2]
        else:
            path = []
            raw_path = list(nx.shortest_path(self.settings_graph, part1, part2))
#            print('rp', raw_path)
            for part1, linkage, part2 in zip(raw_path[:-2:2], raw_path[1::2], raw_path[2::2] + [part2]):
                path.append((part1, linkage, linkage.part1 == part1, part2))

            self._settings_paths[part1, part2] = path
            return path

    def settings_graph(self):
        graph = self.holonomic_graph.copy()
        deleted_linkages = []
        for cycle in nx.cycle_basis(graph):
            ground_distance = [(l, len(nx.shortest_path(graph, l, self.ground)))
                               for l in cycle
                               if l in self.linkages
                               and l not in deleted_linkages
                               and not l.positions_require_kinematic_parameters
                               ]
#            print(ground_distance)
            linkage_to_delete = max(ground_distance, key=lambda x: x[1])[0]
#            print(linkage_to_delete)
            deleted_linkages.append(linkage_to_delete)
            graph.remove_node(linkage_to_delete)

        self.linkages_kinematic_setting = [l for l in self.linkages if l not in deleted_linkages]
        self.settings_graph = graph

    def part_linkages(self):
        part_linkages = {}
        for linkage in self.linkages:
            for part in (linkage.part1, linkage.part2):
                if part not in part_linkages:
                    part_linkages[part] = [linkage]
                else:
                    part_linkages[part].append(linkage)
        return part_linkages

    def plot_graph(self):
        s = """<html>
        <head>
        <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.20.0/vis.min.js"></script>
        <link href="https://cdnjs.cloudflare.com/ajax/libs/vis/4.20.0/vis.min.css" rel="stylesheet" type="text/css" />

        <style type="text/css">
            #mynetwork {
                border: 1px solid lightgray;
            }
        </style>
    </head>
    <body>
    <div id="mynetwork"></div>

    <script type="text/javascript">
    var nodes = new vis.DataSet([\n"""
        index = {}
        for ipart, part in enumerate(self.parts + [self.ground]):
            index[part] = ipart
            s += "{{id: {}, label: '{}'}},\n".format(ipart, part.name)
#        s+=']);\n'
        n = len(self.parts) + 1
#        index[self.ground]=n
#        n+=1
        for il, linkage in enumerate(self.linkages):
            index[linkage] = n + il
            s += "{{id: {}, label: '{}'}},\n".format(n + il, linkage.name)
        s += ']);\n'

        s += "var edges = new vis.DataSet(["
        for linkage in self.linkages:
            s += '{{from: {}, to: {}}},\n'.format(index[linkage], index[linkage.part1])
            s += '{{from: {}, to: {}}},\n'.format(index[linkage], index[linkage.part2])
        s += ']);'

        s += """
    // create a network
    var container = document.getElementById('mynetwork');

    // provide the data in the vis format
    var data = {
        nodes: nodes,
        edges: edges
    };
    var options = {};

    // initialize your network!
    var network = new vis.Network(container, data, options);
</script>
</body>
</html>"""

        with open('gm_graph_viz.html', 'w') as file:
            file.write(s)

        webbrowser.open('file://' + os.path.realpath('gm_graph_viz.html'))

    def draw_power_graph(self, return_graph=False):
        """
        Draw graph with tulip
        """
        import matplotlib.pyplot as plt

        # Creating tulip graph
        G = nx.Graph()
#        nodes={}
#        edges={}
        widths = []
        labels = {}
        edges = []
        for part in self.parts:
            G.add_node(part)

            labels[part] = part.name
#        nodes[self.ground]=G.addNode()
        for linkage in self.linkages:
            G.add_node(linkage)
            labels[linkage] = linkage.name
#            for part in [linkage.part1,linkage.part2]:
            G.add_edge(linkage, linkage.part1)
            widths.append(abs(self.transmitted_linkage_power(linkage, 0)))
            edges.append((linkage, linkage.part1))
            G.add_edge(linkage, linkage.part2)
            widths.append(abs(self.transmitted_linkage_power(linkage, 1)))
            edges.append((linkage, linkage.part2))

        for load in self.unknown_static_loads + self.known_static_loads:
            G.add_node(load)
            G.add_edge(load, load.part)
            widths.append(abs(self.load_power(load)))
            labels[load] = load.name
            edges.append((load, load.part))
        max_widths = max(widths)
#        print(widths)
        widths = [6 * w / max_widths for w in widths]
#        edges[linkage,part]=e
        if not return_graph:
            plt.figure()
            pos = nx.spring_layout(G)
            nx.draw_networkx_nodes(G, pos, nodelist=self.linkages, node_color='grey')
            nx.draw_networkx_nodes(G, pos, nodelist=self.parts)
            nx.draw_networkx_nodes(G, pos, nodelist=self.unknown_static_loads, node_color='red')
            nx.draw_networkx_nodes(G, pos, nodelist=self.known_static_loads, node_color='green')
            nx.draw_networkx_nodes(G, pos, nodelist=self.parts, node_color='cyan')
            nx.draw_networkx_labels(G, pos, labels)
            nx.draw_networkx_edges(G, pos, edges, width=widths, edge_color='blue')
            nx.draw_networkx_edges(G, pos)
        else:
            return G
#        nx.draw_networkx_labels(G,pos,labels)

    def change_imposed_speeds(self, imposed_speeds):
        self.imposed_speeds = imposed_speeds
        self._utd_kinematic_results = False
        self._utd_static_results = False  # Change of speeds might change resistant forces: update needed

    def change_loads(self, known_static_loads, unknown_static_loads):
        self.known_static_loads = known_static_loads
        self.unknown_static_loads = unknown_static_loads
        self._utd_static_results = False

    def speeds(self, position, part_ref, part):
        """
        Speeds from point belonging to part with part_ref as reference
        """
        q = self.kinematic_vector  # force kdof computation
        path = nx.shortest_path(self.holonomic_graph, part_ref, part)
        V = npy.zeros((3, self.n_kdof))
        W = npy.zeros((3, self.n_kdof))
        for il, linkage2 in enumerate(path):
            if not linkage2.__class__.__name__ == 'Part':
                try:
                    if path[il + 1] == linkage2.part2:
                        side = 1
                    else:
                        side = -1
                except IndexError:
                    # linkage is last element of list
                    if path[il - 1] == linkage2.part1:
                        side = 1
                    else:
                        side = -1
                # It's really a linkage
                P = geometry.euler_2_transfer_matrix(*linkage2.euler_angles)
                u = linkage2.position
                uprime = u - position
                L = geometry.cross_product_matrix(uprime)
                Ve = npy.dot(L, npy.dot(P, linkage2.kinematic_matrix[:3, :])
                             ) + npy.dot(P, linkage2.kinematic_matrix[3:, :])
                We = npy.dot(P, linkage2.kinematic_matrix[:3, :])
                for indof, ndof in enumerate(self.kdof[linkage2]):
                    V[:, ndof] += side * Ve[:, indof]
                    W[:, ndof] += side * We[:, indof]

        return npy.hstack((npy.dot(W, q).flatten(), npy.dot(V, q).flatten()))

    def local_linkage_speeds(self, linkage):
        """
        :returns: relative speed of linkage in local linkage coordinate system
        """

        try:
            # holonomic linkage
            r = self.kinematic_results[linkage]  # results
            vr = npy.array([r[i] for i in range(linkage.n_kinematic_unknowns)])  # vector of results
            return npy.dot(linkage.kinematic_matrix, vr).flatten()
        except:
            s = self.Speed(linkage.position, self.ground, linkage.part1)
            v = npy.dot(linkage.P, s[:3])
            w = npy.dot(linkage.P, s[3:])
#            return npy.dot(linkage.P,s)
            return npy.hstack((w, v)).flatten()

    def local_linkage_forces(self, linkage, num_part):
        """
        :returns :Local forces of linkage on part given by its number (0 for part1, 1 for part2)
        """
        r = self.static_results[linkage]  # results
        vr = npy.array([r[i] for i in range(linkage.n_static_unknowns)])  # vector of results
        if num_part == 0:
            return npy.dot(linkage.static_matrix1, vr)
        else:

            return npy.dot(linkage.static_matrix2, vr)

    def global_linkage_forces(self, linkage, num_part):
        P = geometry.euler_2_transfer_matrix(*linkage.euler_angles)
        lf = self.local_linkage_forces(linkage, num_part)
        F = npy.zeros(6)
        F[:3] = npy.dot(P, lf[:3])
        F[3:] = npy.dot(P, lf[3:])
        return F

    def local_load_forces(self, load):
        try:
            F = npy.zeros(6)
            F[:3] = load.forces
            F[3:] = load.torques
            return F
        except AttributeError:
            r = self.static_results[load]  # results
            vr = npy.array([r[i] for i in range(load.n_static_unknowns)])  # vector of results
            return npy.dot(load.static_matrix, vr)

    def global_load_forces(self, load):
        f = self.local_load_forces(load)
#        print(f)
        F = npy.zeros(6)
        F[:3] = npy.dot(load.P, f[:3]).flatten()
#        print(npy.dot(load.P,f[3:]).flatten())
        F[3:] = npy.dot(load.P, f[3:]).flatten()
        return F

    def linkage_power_losses(self, linkage):
        if linkage.holonomic:
            tf = self.local_linkage_forces(linkage, 0)
            ft = npy.hstack((tf[3:], tf[:3]))
            return npy.dot(self.local_linkage_speeds(linkage), ft)
        else:
            d = self.global_linkage_forces(linkage, 0) + self.global_linkage_forces(linkage, 1)
            df = d[:3]
            dt = d[3:]
            s = -self.speeds(linkage.position, self.ground, linkage.part1)
            w = s[:3]
            v = s[3:]
#            print(df,dt,w,v)
            return npy.dot(df, v) + npy.dot(dt, w)

    def load_power(self, load):
        s = self.speeds(load.position, self.ground, load.part)
        f = self.global_load_forces(load)
        P = npy.dot(s[3:], f[:3])
        P += npy.dot(f[3:], s[:3])
        return P

    def transmitted_linkage_power(self, linkage, num_part):
        if num_part == 1:
            s = self.speeds(linkage.position, self.ground, linkage.part2)
        else:
            s = self.speeds(linkage.position, self.ground, linkage.part1)
        f = self.global_linkage_forces(linkage, num_part)
        P = npy.dot(s[3:], f[:3])
        P += npy.dot(f[3:], s[:3])
        return P

    def _get_static_results(self):
        if not self._utd_static_results:
            self._static_results = self._static_analysis()
            self._utd_static_results = True
        return self._static_results

    static_results = property(_get_static_results)

    def _get_kinematic_results(self):
        if not self._utd_kinematic_results:
            self._kinematic_results, self._kinematic_vector = self._kinematic_analysis()
            self._utd_kinematic_results = True
        return self._kinematic_results

    kinematic_results = property(_get_kinematic_results)

    def _get_kinematic_vector(self):
        if not self._utd_kinematic_results:
            self._kinematic_results, self._kinematic_vector = self._kinematic_analysis()
            self._utd_kinematic_results = True
        return self._kinematic_vector

    kinematic_vector = property(_get_kinematic_vector)

    def behavior_equation_of_unknows_loads_static(self, M, K, nonlinear_eq, neq, neq_linear, indices_r):

        for load in self.unknown_static_loads:
            neq_load = load.static_behavior_occurence_matrix.shape[0]
            if neq_load > 0:
                # Adding to occurence matrix
                Me = npy.zeros((neq_load, self.n_sdof))
                for indof, ndof in enumerate(self.sdof[load]):
                    Me[:, ndof] += load.static_behavior_occurence_matrix[:, indof]
                M = npy.vstack([M, Me])

                # Adding linear equations to System matrix K
                neq_linear_load = load.static_behavior_linear_eq.shape[0]
                if neq_linear_load > 0:
                    Ke = npy.zeros((neq_linear_load, self.n_sdof))
                    for indof, ndof in enumerate(self.sdof[load]):
                        Ke[neq_linear:neq_linear + 6, ndof] += load.static_behavior_linear_eq[:, indof]
                    K = npy.vstack([K, Ke])

                    indices_r.extend(range(neq_linear, neq_linear + neq_linear_load))

                # Collecting non linear equations

                if load.static_require_kinematic:
                    # Absolute speed in this case in local coordinate system
                    s = self.speeds(load.position, self.ground, load.part)
                    w = npy.dot(load.P.T, s[:3])
                    v = npy.dot(load.P.T, s[3:])
                else:
                    w, v = 0, 0
                for i, fct in zip(load.static_behavior_nonlinear_eq_indices, load.static_behavior_nonlinear_eq):
                    nonlinear_eq[neq + i] = lambda x, v=v, w=w, fct=fct: fct(x, w, v)

                # Updating counters
                neq += neq_load
                neq_linear += neq_linear_load
        return M, K

    def behavior_equation_of_linkages_static(self, M, K, nonlinear_eq, neq, neq_linear, indices_r):

        for linkage in self.linkages:
            neq_linkage = linkage.static_behavior_occurence_matrix.shape[0]
            if neq_linkage > 0:
                # Adding to occurence matrix
                Me = npy.zeros((neq_linkage, self.n_sdof))
                for indof, ndof in enumerate(self.sdof[linkage]):
                    Me[:, ndof] += linkage.static_behavior_occurence_matrix[:, indof]
                M = npy.vstack([M, Me])

                # Adding linear equations to System matrix K
                neq_linear_linkage = linkage.static_behavior_linear_eq.shape[0]
                if neq_linear_linkage > 0:
                    Ke = npy.zeros((neq_linear_linkage, self.n_sdof))
                    for indof, ndof in enumerate(self.sdof[linkage]):
                        Ke[neq_linear:neq_linear + 6, ndof] += linkage.static_behavior_linear_eq[:, indof]
                    K = npy.vstack([K, Ke])
                    indices_r.extend(range(neq_linear, neq_linear + neq_linear_linkage))

                # Collecting non linear equations

                if linkage.holonomic:
                    if linkage.static_require_kinematic:
                        # Relatives speed in this case
                        ls = self.local_linkage_speeds(linkage)
                        w = ls[:3]
                        v = ls[3:]
                    else:
                        w, v = 0, 0
                    for i, fct in zip(linkage.static_behavior_nonlinear_eq_indices,
                                      linkage.static_behavior_nonlinear_eq):
                        nonlinear_eq[neq + i] = lambda x, v=v, w=w, fct=fct: fct(x, w, v)
                else:
                    if linkage.static_require_kinematic:
                        # Absolute speed in this case in local coordinate system
                        s = self.speeds(linkage.position, self.ground, linkage.part1)
                        w = npy.dot(linkage.P.T, s[:3])
                        v = npy.dot(linkage.P.T, s[3:])
                    else:
                        w, v = 0, 0

                    for i, fct in zip(linkage.static_behavior_nonlinear_eq_indices,
                                      linkage.static_behavior_nonlinear_eq):
                        nonlinear_eq[neq + i] = lambda x, v=v, w=w, fct=fct: fct(x, w, v)

                # Updating counters
                neq += neq_linkage
                neq_linear += neq_linear_linkage

        return M, K

    def occurence_matrix_assembly_static(self, M, K, F, uloads_parts, loads_parts):

        for ip, part in enumerate(self.parts):
            # linkage contribution to the LHS
            for linkage in self.graph[part].keys():
                P = geometry.euler_2_transfer_matrix(*linkage.euler_angles)
                u = linkage.position
                uprime = u
                L = geometry.cross_product_matrix(uprime)
                if part == linkage.part1:
                    static_matrix = linkage.static_matrix1
                else:
                    static_matrix = linkage.static_matrix2
                Me1 = npy.abs(npy.dot(P, static_matrix[:3, :])) > 1e-10
                Me2 = npy.abs(npy.dot(L, npy.dot(P, static_matrix[:3, :])) + npy.dot(P, static_matrix[3:, :])) > 1e-10

                Me = npy.vstack([Me1, Me2])
                for indof, ndof in enumerate(self.sdof[linkage]):
                    M[ip * 6:(ip + 1) * 6, ndof] += Me[:, indof]

                Ke1 = npy.dot(P, static_matrix[:3, :])
                Ke2 = npy.dot(L, npy.dot(P, static_matrix[:3, :])) + npy.dot(P, static_matrix[3:, :])

                Ke = npy.vstack([Ke1, Ke2])

                for indof, ndof in enumerate(self.sdof[linkage]):
                    K[ip * 6:(ip + 1) * 6, ndof] += Ke[:, indof]

            # Unknowns loads contribution to the LHS
            try:
                uloads = uloads_parts[part]
            except:
                uloads = []
            for load in uloads:
                P = geometry.euler_2_transfer_matrix(*load.euler_angles)
                u = load.position
                uprime = u
                L = geometry.cross_product_matrix(uprime)
                Me1 = npy.abs(npy.dot(P, load.static_matrix[:3, :])) > 1e-10
                Me2 = npy.abs(
                    npy.dot(L, npy.dot(P, load.static_matrix[:3, :])) + npy.dot(P, load.static_matrix[3:, :])) > 1e-10
                Me = npy.vstack([Me1, Me2])
                for indof, ndof in enumerate(self.sdof[load]):
                    M[ip * 6:(ip + 1) * 6, ndof] += Me[:, indof]

                Ke1 = npy.dot(P, load.static_matrix[:3, :])
                Ke2 = npy.dot(L, npy.dot(P, load.static_matrix[:3, :])) + npy.dot(P, load.static_matrix[3:, :])
                Ke = npy.vstack([Ke1, Ke2])
                for indof, ndof in enumerate(self.sdof[load]):
                    K[ip * 6:(ip + 1) * 6, ndof] = Ke[:, indof]  # minus because of sum is in LHS

            # knowns loads contribution to the RHS
            try:
                loads = loads_parts[part]
            except:
                loads = []
            for load in loads:
                P = geometry.euler_2_transfer_matrix(*load.euler_angles)
                u = load.position
                uprime = u
                L = geometry.cross_product_matrix(uprime)
                F1 = npy.dot(P, load.forces)
                F2 = npy.dot(L, npy.dot(P, load.forces)) + npy.dot(P, load.torques)
                Fe = npy.hstack([F1.reshape(3), F2.reshape(3)])
                F[ip * 6:(ip + 1) * 6] -= Fe  # minus because of sum is in LHS

        return M, K, F

    def construction_matrix_m_k_f_and_non_linear_eq_static(self, uloads_parts, loads_parts):

        lparts = len(self.parts)
        M = npy.zeros((6 * lparts, self.n_sdof))
        K = npy.zeros((6 * lparts, self.n_sdof))
        F = npy.zeros(6 * lparts)
        nonlinear_eq = {}

        # Occurence matrix assembly
        M, K, F = self.occurence_matrix_assembly_static(M, K, F, uloads_parts, loads_parts)
        neq = 6 * lparts
        neq_linear = neq
        indices_r = list(range(neq))

        # behavior equations of linkages
        M, K = self.behavior_equation_of_linkages_static(M, K, nonlinear_eq, neq, neq_linear, indices_r)

        # behavior equations of unknowns loads
        M, K = self.behavior_equation_of_unknows_loads_static(M, K, nonlinear_eq, neq, neq_linear, indices_r)

        return M, K, F, nonlinear_eq, indices_r

    def construction_matrix_q_static(self, M, K, F, nonlinear_eq, indices_r, resolution_order):
        q = npy.zeros(self.n_sdof)

        for eqs, variables in resolution_order:
            #            print(eqs,variables)
            linear = True
            linear_eqs = []
            for eq in eqs:
                try:
                    nonlinear_eq[eq]
                    linear = False
                except KeyError:
                    linear_eqs.append(eq)
            if linear:
                eqs_r = npy.array([indices_r[eq] for eq in eqs])
                other_vars = npy.array([i for i in range(self.n_sdof) if i not in variables])
                Kr = K[eqs_r[:, None], npy.array(variables)]
                Fr = F[eqs_r] - npy.dot(K[eqs_r[:, None], other_vars], q[other_vars])

                q[variables] = linalg.solve(Kr, Fr)

            else:

                nl_eqs = []
                other_vars = npy.array([i for i in range(self.n_sdof) if i not in variables])
                for eq in eqs:
                    try:
                        f1 = nonlinear_eq[eq]
                        vars_func = [i for i in range(self.n_sdof) if M[eq, i]]

                        def f2(x, f1=f1, vars_func=vars_func, variables=variables, q=q):
                            x2 = []
                            for variable in vars_func:
                                try:
                                    x2.append(x[variables.index(variable)])
                                except ValueError:
                                    x2.append(q[variable])
                            return f1(x2)

                        nl_eqs.append(f2)
                    except KeyError:
                        # lambdification of linear equations

                        f2 = lambda x, indices_r=indices_r, K=K, F=F, eq=eq, q=q, other_vars=other_vars: npy.dot(
                            K[indices_r[eq], variables], x) - F[indices_r[eq]] + npy.dot(K[indices_r[eq], other_vars],
                                                                                         q[other_vars])

                        nl_eqs.append(f2)
                f = lambda x: [fi(x) for fi in nl_eqs]
                xs = fsolve(f, npy.zeros(len(variables)), full_output=0)

                if npy.sum(npy.abs(f(xs))) > 1e-4:
                    raise ModelError('No convergence of nonlinear phenomena solving' + str(npy.sum(npy.abs(f(xs)))))

                q[variables] = xs
        return q

    def _static_analysis(self):
        # Static parametrisation
        self.sdof = {}
        self.n_sdof = 0
        kinematic_analysis_required = False
        for linkage in self.linkages:
            self.sdof[linkage] = list(range(self.n_sdof, self.n_sdof + linkage.n_static_unknowns))
            self.n_sdof += linkage.n_static_unknowns
            if linkage.static_require_kinematic:
                kinematic_analysis_required = True

        uloads_parts = {}
        for load in self.unknown_static_loads:
            load_unknowns = load.static_matrix.shape[1]
            self.sdof[load] = list(range(self.n_sdof, self.n_sdof + load_unknowns))
            self.n_sdof += load_unknowns
            if load.static_require_kinematic:
                kinematic_analysis_required = True

            # Loading sorting by part
            try:
                uloads_parts[load.part].append(load)
            except:
                uloads_parts[load.part] = [load]

        # Force kinematic computation if required
        if kinematic_analysis_required:
            self.kinematic_results

        # Knowns loads sorting by parts
        loads_parts = {}
        for load in self.known_static_loads:
            try:
                loads_parts[load.part].append(load)
            except:
                loads_parts[load.part] = [load]

        lparts = len(self.parts)
        M, K, F, nonlinear_eq, indices_r = self.construction_matrix_m_k_f_and_non_linear_eq_static(
            uloads_parts, loads_parts)

        solvable, solvable_var, resolution_order = tools.equations_system_analysis(M, None)
#        print(resolution_order)
        if not solvable:
            raise ModelError('Overconstrained system')

        q = self.construction_matrix_q_static(M, K, F, nonlinear_eq, indices_r, resolution_order)
        results = {}
        for link, dofs in self.sdof.items():
            rlink = {}
            for idof, dof in enumerate(dofs):
                rlink[idof] = q[dof]
            results[link] = rlink

        return results

    def construction_matrix_m_k_f_kinematic(self, neq, loops, nhl, ):
        ll = len(loops)
        ieq = 6 * ll
        K = npy.zeros((neq, self.n_kdof))
        F = npy.zeros((neq, 1))
        M = npy.zeros((neq, self.n_kdof))

        for il, loop in enumerate(loops):
            for ilk, linkage in enumerate(loop):
                if not linkage.__class__.__name__ == 'Part':
                    try:
                        if loop[ilk + 1] == linkage.part2:
                            side = 1
                        else:
                            side = -1
                    except IndexError:
                        # linkage is last element of list
                        if loop[ilk - 1] == linkage.part1:
                            side = 1
                        else:
                            side = -1
                    # Linkage
                    P = geometry.euler_2_transfer_matrix(*linkage.euler_angles)
                    u = linkage.position
                    uprime = u
                    L = geometry.cross_product_matrix(uprime)
                    Me1 = npy.abs(npy.dot(P, linkage.kinematic_matrix[:3, :])) > 1e-10
                    Me2 = npy.abs(npy.dot(L, npy.dot(P, linkage.kinematic_matrix[:3, :])) + npy.dot(P,
                                                                                                    linkage.kinematic_matrix[
                                                                                                        3:, :])) > 1e-10
                    Me = npy.vstack([Me1, Me2])
                    for indof, ndof in enumerate(self.kdof[linkage]):
                        M[6 * il:6 * il + 6, ndof] += Me[:, indof]

                    Ke1 = npy.dot(P, linkage.kinematic_matrix[:3, :])
                    Ke2 = npy.dot(L, npy.dot(P, linkage.kinematic_matrix[:3, :])) + npy.dot(P,
                                                                                            linkage.kinematic_matrix[3:,
                                                                                                                     :])
                    Ke = side * npy.vstack([Ke1, Ke2])
                    for indof, ndof in enumerate(self.kdof[linkage]):
                        K[6 * il:6 * il + 6, ndof] += Ke[:, indof]

        # Non holonomic equations
        for linkage in nhl:
            # Speed computation
            try:
                path = nx.shortest_path(self.holonomic_graph,
                                        linkage.part1,
                                        linkage.part2)
            except nx.NetworkXNoPath:
                raise ModelError(
                    'No path between {} and {} for linkage {} of type {}'.format(linkage.part1.name, linkage.part2.name,
                                                                                 linkage.name,
                                                                                 linkage.__class__.__name__))

            V = npy.zeros((3, self.n_kdof))
            for il, linkage2 in enumerate(path):
                if not linkage2.__class__.__name__ == 'Part':
                    try:
                        if path[il + 1] == linkage2.part2:
                            side = 1
                        else:
                            side = -1
                    except IndexError:
                        # linkage is last element of list
                        if path[il - 1] == linkage2.part1:
                            side = 1
                        else:
                            side = -1
                    # It's really a linkage
                    P = geometry.euler_2_transfer_matrix(*linkage2.euler_angles)
                    u = linkage2.position
                    uprime = u - linkage.position
                    L = geometry.cross_product_matrix(uprime)
                    Ve = npy.dot(L, npy.dot(P, linkage2.kinematic_matrix[:3, :])) + npy.dot(P,
                                                                                            linkage2.kinematic_matrix[
                                                                                                3:, :])
                    for indof, ndof in enumerate(self.kdof[linkage2]):
                        V[:, ndof] += side * Ve[:, indof]

            P = geometry.euler_2_transfer_matrix(*linkage.euler_angles)
            for direction in linkage.kinematic_directions:
                K[ieq, :] = npy.dot(npy.dot(P, direction), V)
                ieq += 1

        # Imposed speeds equations
        for linkage, index, speed in self.imposed_speeds:
            K[ieq, self.kdof[linkage][index]] = 1
            F[ieq, 0] = speed
            ieq += 1

        return M, K, F

    def _kinematic_analysis(self):
        # Kinematic setting
        self.n_kdof = 0
        self.kdof = {}
        loops = nx.cycle_basis(self.holonomic_graph)
        ll = len(loops)
        ieq = 6 * ll
        neq = ieq + len(self.imposed_speeds)
        nhl = []
        for linkage in self.linkages:
            try:
                # Holonomic linkage
                ndofl = linkage.kinematic_matrix.shape[1]
                for dofl in range(ndofl):
                    try:
                        self.kdof[linkage].append(self.n_kdof + dofl)
                    except KeyError:
                        self.kdof[linkage] = [self.n_kdof + dofl]
                self.n_kdof += ndofl
            except AttributeError:
                # Non holonomic linkage
                nhl.append(linkage)
                neq += len(linkage.kinematic_directions)

        q = npy.zeros((self.n_kdof, 1))
        M, K, F = self.construction_matrix_m_k_f_kinematic(neq, loops, nhl)

        # deducing M from K for last lines

        M[6 * ll:, :] = npy.abs(K[6 * ll:, :]) > 1e-10

        solvable, solvable_var, resolution_order = tools.equations_system_analysis(M, None)

        if solvable:
            for eqs, variables in resolution_order:
                eqs = npy.array(eqs)
                other_vars = npy.array([i for i in range(self.n_kdof) if i not in variables])
                Kr = K[eqs[:, None], npy.array(variables)]
                Fr = F[eqs, :] - npy.dot(K[eqs[:, None], other_vars], q[other_vars, :])
                q[variables, :] = linalg.solve(Kr, Fr)
            results = {}
            for link, dofs in self.kdof.items():
                rlink = {}
                for idof, dof in enumerate(dofs):
                    rlink[idof] = q[dof, 0]
                results[link] = rlink
            return results, q
        else:
            raise ModelError('Overconstrained system')

    def global_sankey(self):
        from matplotlib.sankey import Sankey
        flows = []
        orientations = []
        labels = []

        for load in self.known_static_loads + self.unknown_static_loads:
            pl = self.load_power(load)
#            print(pl)
            flows.append(pl)
            labels.append(load.name)
            if (load.__class__.__name__ == 'SimpleUnknownLoad') | (load.__class__.__name__ == 'KnownLoad'):
                orientations.append(0)
            else:
                orientations.append(-1)

        for linkage in self.linkages:
            pl = self.linkage_power_losses(linkage)
#            print(pl)
            if pl != 0.:
                flows.append(-pl)
                orientations.append(-1)
                labels.append(linkage.name)
#            pl.append(0.5*l)

        l = max([abs(f) for f in flows])
#        pl=[0.1*l]*len(flows)

        sankey = Sankey(unit='W', scale=1 / l)

        sankey.add(flows=flows,
                   orientations=orientations, labels=labels,)

        sankey.finish()

    def plot(self, u=1):
        points = []
        for linkage in self.linkages:
            points.append(vm.Point3D(linkage.position[0], linkage.position[1], linkage.position[2]))
            list_result = linkage.position + u * linkage.P[:, 0]
            points.append(vm.edges.Line3D(vm.Point3D(linkage.position[0], linkage.position[1], linkage.position[2]),
                                          vm.Point3D(list_result[0], list_result[1], list_result[2])))
        mdl = vm.core.VolumeModel(points)
        mdl.plot(equal_aspect=False)

    def scene_caracteristics(self):
        min_vect = self.linkages[0].position.copy()
        max_vect = self.linkages[0].position.copy()
        center = self.linkages[0].position.copy()
        n = 1
        for linkage in self.linkages[1:] + self.known_static_loads + self.unknown_static_loads:
            # print(linkage)
            for i, (xmin, xmax, xi) in enumerate(zip(min_vect, max_vect, linkage.position)):
                # print(i,xmin,xmax,xi)
                if xi < xmin:
                    min_vect[i] = xi
                    # print('min',min_vect)
                if xi > xmax:
                    max_vect[i] = xi
                    # print('max',max_vect)
            # print(linkage,linkage.position)
            center += linkage.position
            n += 1

        center = center / n
        # print(min_vect,max_vect)
        max_length = linalg.norm(min_vect - max_vect)
        # print(center,max_length)
        return center, max_length

    def babylon_script(self, forces=True):

        center, length = self.scene_caracteristics()

        if forces:
            max_force = 0.
            max_torque = 0
            for linkage in self.linkages:
                f = self.global_linkage_forces(linkage, 0)[0:3]
                t = self.global_linkage_forces(linkage, 0)[3:]
                max_force = max(max_force, max(f))
                max_torque = max(max_torque, max(t))

            for load in self.known_static_loads + self.unknown_static_loads:
                f = self.global_load_forces(load)[0:3]
                t = self.global_load_forces(load)[3:]
                max_force = max(max_force, max(f))
                max_torque = max(max_torque, max(t))

#            print(max_force,max_torque)

        linkages_strings = []
        for linkage in self.linkages:
            try:
                if forces:
                    # print([f/max_force/length/4 for f in self.GlobalLoadForces(load)[0:3]])
                    linkages_strings.append(
                        linkage.babylon(length, [f / max_force * length / 4 for f in self.global_linkage_forces(
                            linkage, 0)[0:3]], [f / max_torque * length / 4 for f in self.global_linkage_forces(
                                linkage, 0)[3:]]))
                else:
                    linkages_strings.append(linkage.babylon(length, None, None))

            except AttributeError:
                print('error')
                pass

        return babylon_template.substitute(linkages_strings=linkages_strings,
                                           length=length,
                                           center=center,
                                           name=self.name)

    def babylonjs(self, page='gm_babylonjs'):
        page += '.html'
        with open(page, 'w') as file:
            file.write(self.babylon_script())

        webbrowser.open('file://' + os.path.realpath(page))

#    def BabylonJS(self):
#        center,length=self.SceneCaracteristics()
#        s="""
#        <!doctype html>
# <html>
# <head>
#   <meta charset="utf-8">
#   <title>Babylon - Basic scene</title>
#   <style>
#      html, body {
#         overflow: hidden;
#         width: 100%;
#         height: 100%;
#         margin: 0;
#         padding: 0;
#      }
#      #renderCanvas {
#         width: 100%;
#         height: 100%;
#         touch-action: none;
#      }
#   </style>
#   <script src="https://cdnjs.cloudflare.com/ajax/libs/babylonjs/2.5.0/babylon.js"></script>
#   <script src="https://cdnjs.cloudflare.com/ajax/libs/handjs/1.3.11/hand.js"></script>
#   <script src="https://cdnjs.cloudflare.com/ajax/libs/cannon.js/0.6.2/cannon.min.js"></script> <!-- optional physics engine -->
# </head>
# <body>
#   <canvas id="renderCanvas"></canvas>
#   <script type="text/javascript">
#      // Get the canvas element from our HTML below
#      var canvas = document.querySelector("#renderCanvas");
#      // Load the BABYLON 3D engine
#      var engine = new BABYLON.Engine(canvas, true);
#      // -------------------------------------------------------------
#      // Here begins a function that we will 'call' just after it's built
#      var createScene = function () {
#         // Now create a basic Babylon Scene object
#         var scene = new BABYLON.Scene(engine);
#         // This creates and positions a free camera
#         //var plane = BABYLON.Mesh.CreatePlane("plane", 10.0, scene);
#         """
#
#
##        s+='var camera = new BABYLON.FreeCamera("camera1", new BABYLON.Vector3({}, {}, {}), scene);'.format(*center)
#        s+='var camera = new BABYLON.ArcRotateCamera("ArcRotateCamera", 0., 0., {}, new BABYLON.Vector3({}, {}, {}), scene);\n'.format(4*length,*center)
##        s+='var camera = new BABYLON.ArcRotateCamera("ArcRotateCamera", 0., 0., 15, new BABYLON.Vector3(0,0,0), scene);'
#        s+='camera.panningSensibility={};\n'.format(30)
#        s+='camera.pinchPrecision={};\n'.format(20)
#        s+='camera.wheelPrecision={};\n'.format(20)
#        s+="""
#        camera.attachControl(canvas, false);
#         // This creates a light, aiming 0,1,0 - to the sky.
#         var light = new BABYLON.HemisphericLight("light1", new BABYLON.Vector3(0, 1, 0), scene);
#         // Dim the light a small amount
#         light.intensity = .5;
#         """
#
# // Let's try our built-in 'sphere' shape. Params: name, subdivisions, size, scene
# var sphere = BABYLON.Mesh.CreateSphere("sphere1", 16, 2, scene);
# // Move the sphere upward 1/2 its height
##         sphere.position.y = 1;
# // Let's try our built-in 'ground' shape. Params: name, width, depth, subdivisions, scene
# var ground = BABYLON.Mesh.CreateGround("ground1", 6, 6, 2, scene);
#
#        for linkage in self.linkages:
#            s+='var sphere = BABYLON.Mesh.CreateSphere("sphere1", 15., {}, scene);\n'.format(length/30)
#            s+="sphere.position=new BABYLON.Vector3({},{},{});\n".format(*linkage.position)
# s+="var lines = BABYLON.Mesh.CreateLines("lines", [    new BABYLON.Vector3(-10, 0, 0)"
#        s+="""
#         // Leave this function
#         return scene;
#      }; // End of createScene function
#      // -------------------------------------------------------------
#      // Now, call the createScene function that you just finished creating
#      var scene = createScene();
#      // Register a render loop to repeatedly render the scene
#      engine.runRenderLoop(function () {
#         scene.render();
#      });
#      // Watch for browser/canvas resize events
#      window.addEventListener("resize", function () {
#         engine.resize();
#      });
#   </script>
# </body>
# </html>
#
#        """
#        with open('gm_babylonjs.html','w') as file:
#            file.write(s)
#
#        webbrowser.open('file://' + os.path.realpath('gm_babylonjs.html'))
