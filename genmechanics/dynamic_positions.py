#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""
import os
import webbrowser

import math
import random
import numpy as npy
import cma

from matplotlib.colors import hsv_to_rgb
import matplotlib.pyplot as plt

from matplotlib.patches import Arrow

import networkx as nx

from jinja2 import Environment, PackageLoader, select_autoescape

from dessia_common.core import DessiaObject
#from numpy import zeros
from scipy.optimize import minimize
import volmdlr as vm
from genmechanics.core import Part, Mechanism


class Parameter(DessiaObject):
    def __init__(self, lower_bound, upper_bound, periodicity=None):
        DessiaObject.__init__(self,
                              lower_bound=lower_bound,
                              upper_bound=upper_bound,
                              periodicity=periodicity)
    def random_value(self):
        return random.uniform(self.lower_bound, self.upper_bound)

    def are_values_equal(self, value1, value2, tol=1e-3):
        if self.periodicity is not None:
            value1 = value1 % self.periodicity
            value2 = value2 % self.periodicity

        return math.isclose(value1, value2, abs_tol=tol)


    def optimizer_bounds(self):
        if self.periodicity is not None:
            return (self.lower_bound-0.5*self.periodicity,
                    self.upper_bound+0.5*self.periodicity)
        else:
            return (self.lower_bound, self.upper_bound)


class Linkage(DessiaObject):
    _non_serializable_attributes = ['part1_position_function',
                                    'part2_position_function',
                                    'part1_basis_function',
                                    'part2_basis_function']
    
    def __init__(self,
                 part1, part1_position_function, part1_basis_function,
                 part2, part2_position_function, part2_basis_function,
                 positions_require_kinematic_parameters,
                 basis_require_kinematic_parameters,
                 kinematic_parameters,
                 name=''):
        """

        """
        DessiaObject.__init__(self,
                              part1=part1,
                              part1_position_function=part1_position_function,
                              part1_basis_function=part1_basis_function,
                              part2=part2,
                              part2_position_function=part2_position_function,
                              part2_basis_function=part2_basis_function,
                              positions_require_kinematic_parameters=positions_require_kinematic_parameters,
                              basis_require_kinematic_parameters=basis_require_kinematic_parameters,
                              kinematic_parameters=kinematic_parameters,
                              number_kinematic_parameters=len(kinematic_parameters),
                              name=name)
                
        
    def is_equivalent(self, other_linkage):
        if self.__class__ != other_linkage.__class__:
            return False

        if hasattr(self, 'part1_position'):
            if self.part1_position != other_linkage.part1_position:
                return False
            
        if hasattr(self, 'part2_position'):
            if self.part2_position != other_linkage.part2_position:
                return False
            
        if hasattr(self, 'part1_basis'):
            if self.part1_basis != other_linkage.part1_basis:
                return False
            
        if hasattr(self, 'part2_basis'):
            if self.part2_basis != other_linkage.part2_basis:
                return False
            
        return True
        


    def frame(self, linkage_parameters_values, side):
        if side:

            part1_frame = self.part1_basis_function(linkage_parameters_values)\
                .to_frame(self.part1_position_function(linkage_parameters_values))

            part2_frame = -self.part2_basis_function(linkage_parameters_values)\
                .to_frame(-self.part2_position_function(linkage_parameters_values))

            return part1_frame + part2_frame

        else:
            part1_frame = -self.part1_basis_function(linkage_parameters_values)\
                .to_frame(-self.part1_position_function(linkage_parameters_values))

            part2_frame = self.part2_basis_function(linkage_parameters_values)\
                .to_frame(self.part2_position_function(linkage_parameters_values))

            return part2_frame + part1_frame

    def babylonjs(self, initial_linkage_parameters,
                  part1_parent=None, part2_parent=None):
        part1_position = self.part1_position_function(initial_linkage_parameters)
        part2_position = self.part2_position_function(initial_linkage_parameters)

        s = ''
        if part1_parent is not None:
            s += 'var linkage_part1_mesh = BABYLON.MeshBuilder.CreateSphere("default_linkage part1", {diameter: 0.02}, scene);\n'
            s += 'linkage_part1_mesh.position = new BABYLON.Vector3({}, {}, {});\n'.format(*part1_position.vector)
            s += 'linkage_part1_mesh.parent = {};\n'.format(part1_parent)

        if part2_parent:
            s += 'var linkage_part2_mesh = BABYLON.MeshBuilder.CreateSphere("default_linkage part2", {diameter: 0.02}, scene);\n'
            s += 'linkage_part2_mesh.position = new BABYLON.Vector3({}, {}, {});\n'.format(*part2_position.vector)
            s += 'linkage_part2_mesh.parent = {};\n'.format(part2_parent)

        return s


class RevoluteLinkage(Linkage):
    holonomic = True

    def __init__(self,
                 part1, part1_position, part1_basis,
                 part2, part2_position, part2_basis, name='RevoluteLinkage'):
        """
        :param part2_basis: a basis defining orientation of linkage on part2

        """

        def part1_basis_f(q):
            return part1_basis.Rotation(part1_basis.u, q[0], copy=True)

        def part2_basis_f(q):
            return part2_basis

        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

        Linkage.__init__(self,
                         part1, lambda q: part1_position, part1_basis_f,
                         part2, lambda q: part2_position, part2_basis_f,
                         False, True,
                         [Parameter(0., 2*math.pi, 2*math.pi)], name=name)



    def babylonjs(self, initial_linkage_parameters,
                  part1_parent=None, part2_parent=None):
        s = ''
        if part1_parent is not None:

            s += 'var path1 = [new BABYLON.Vector3({}, {}, {}), new BABYLON.Vector3({}, {}, {})];\n'.format(*(self.part1_position-0.03*self.part1_basis.u),
                                                                                                           *(self.part1_position+0.03*self.part1_basis.u))
            s += 'var linkage_part1_mesh = BABYLON.MeshBuilder.CreateTube("revolute part1", {path: path1, radius: 0.01, sideOrientation:BABYLON.Mesh.DOUBLESIDE}, scene);\n'
            s += 'linkage_part1_mesh.enableEdgesRendering();\n'
            s += 'linkage_part1_mesh.edgesWidth = 0.4;\n'
            s += 'linkage_part1_mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 1);\n'
            s += 'linkage_part1_mesh.parent = {};\n'.format(part1_parent)


        if part2_parent is not None:
            s += 'var path2 = [new BABYLON.Vector3({}, {}, {}), new BABYLON.Vector3({}, {}, {})];\n'.format(*(self.part2_position-0.03*self.part2_basis.u),
                                                                                                           *(self.part2_position+0.03*self.part2_basis.u))
            s += 'var linkage_part2_mesh = BABYLON.MeshBuilder.CreateTube("revolute part2", {path: path2, radius: 0.015, sideOrientation:BABYLON.Mesh.DOUBLESIDE}, scene);\n'
            s += 'linkage_part2_mesh.enableEdgesRendering();\n'
            s += 'linkage_part2_mesh.edgesWidth = 0.4;\n'
            s += 'linkage_part2_mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 1);\n'
            s += 'linkage_part2_mesh.parent = {};\n'.format(part2_parent)

        return s


class SlidingRevoluteLinkage(Linkage):
    holonomic = True

    def __init__(self,
                 part1, part1_position, part1_basis,
                 part2, part2_position, part2_basis, name='SlidingRevoluteLinkage'):
        """
        :param part2_basis: a basis defining orientation of linkage on part2
        The first kineamtic parameter is the translation, the second the rotation
        """

        def part1_position_f(q):
            return part1_position + q[0]*part1_basis.u

        def part2_position_f(q):
            return part2_position

        def part1_basis_f(q):
            return part1_basis.Rotation(part1_basis.u, q[1], copy=True)

        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

        Linkage.__init__(self,
                         part1, part1_position_f, part1_basis_f,
                         part2, part2_position_f, lambda q: part2_basis,
                         True, True,
                         [Parameter(0., 2*math.pi, 2*math.pi),
                          Parameter(-1., 1., None)], name)



    def babylonjs(self, initial_linkage_parameters,
                  part1_parent=None, part2_parent=None):
#        part1_position = self.part1_position_function(initial_linkage_parameters)
#        part1_basis = self.part1_position(initial_linkage_parameters)

#        part2_position = self.part2_position_function(initial_linkage_parameters)
#        part2_basis = self.part2_position(initial_linkage_parameters)
        s = ''
        if part1_parent is not None:

            s += 'var path1 = [new BABYLON.Vector3({}, {}, {}), new BABYLON.Vector3({}, {}, {})];\n'.format(*(self.part1_position-0.1*self.part1_basis.u),
                                                                                                            *(self.part1_position+0.1*self.part1_basis.u))
            s += 'var linkage_part1_mesh = BABYLON.MeshBuilder.CreateTube("revolute part1", {path: path1, radius: 0.01, sideOrientation:BABYLON.Mesh.DOUBLESIDE}, scene);\n'
            s += 'linkage_part1_mesh.enableEdgesRendering();\n'
            s += 'linkage_part1_mesh.edgesWidth = 0.4;\n'
            s += 'linkage_part1_mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 1);\n'
            s += 'linkage_part1_mesh.parent = {};\n'.format(part1_parent)


        if part2_parent is not None:
            s += 'var path2 = [new BABYLON.Vector3({}, {}, {}), new BABYLON.Vector3({}, {}, {})];\n'.format(*(self.part2_position-0.03*self.part2_basis.u),
                                                                                                           *(self.part2_position+0.03*self.part2_basis.u))
            s += 'var linkage_part2_mesh = BABYLON.MeshBuilder.CreateTube("revolute part2", {path: path2, radius: 0.015, sideOrientation:BABYLON.Mesh.DOUBLESIDE}, scene);\n'
            s += 'linkage_part2_mesh.enableEdgesRendering();\n'
            s += 'linkage_part2_mesh.edgesWidth = 0.4;\n'
            s += 'linkage_part2_mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 1);\n'
            s += 'linkage_part2_mesh.parent = {};\n'.format(part2_parent)

        return s

class PrismaticLinkage(Linkage):
    holonomic = True


    def __init__(self,
                 part1, part1_position, part1_basis,
                 part2, part2_position, part2_basis, name='PrismaticLinkage'):
        """
        :param part2_basis: a basis defining orientation of linkage on part2

        """

        def part1_position_f(q):
            return part1_position + q[0]*part1_basis.u

        def part2_position_f(q):
            return part2_position

        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

        Linkage.__init__(self,
                         part1, part1_position_f, lambda q: part1_basis,
                         part2, part2_position_f, lambda q: part2_basis,
                         True, False,
                         [Parameter(-1, 1, None)], name)




#    def __init__(self,
#                 part1, part1_initial_position, part1_axis,
#                 part2, part2_position, part2_basis=vm.xyz.copy(),
#                 name=''):
#        """
#        sliding on part1, fixed on part2
#        """
#
#        def basis(q):
#            return part2_basis
#
#        Linkage.__init__(self,
#                         part1, lambda q: part1_initial_position+q[0]*part1_axis,
#                         part2, lambda q: part2_position,
#                         True, False,
#                         basis, 1, name)
#
#        DessiaObject.__init__(self,
#                              part1_axis=part1_axis,
#                              part2_basis=part2_basis)

    def babylonjs(self, initial_linkage_parameters, part1_parent=None, part2_parent=None):

        bp1 = self.part1_basis_function(initial_linkage_parameters)
        bp2 = self.part2_basis_function(initial_linkage_parameters)
        s = ''
        if part1_parent is not None:
            s += 'var linkage_part1_mesh = BABYLON.MeshBuilder.CreateBox("prismatic part1", {depth:0.015, height:0.015, width:0.25}, scene);\n'
            s += 'linkage_part1_mesh.parent = {};\n'.format(part1_parent)
            s += 'linkage_part1_mesh.position = new BABYLON.Vector3({}, {}, {});\n'.format(*self.part1_position_function(initial_linkage_parameters))
            s += 'linkage_part1_mesh.rotation = BABYLON.Vector3.RotationFromAxis(new BABYLON.Vector3({}, {}, {}),new BABYLON.Vector3({}, {}, {}), new BABYLON.Vector3({}, {}, {}));\n'.format(*bp1.u, *bp1.v, *bp1.w)
            s += 'linkage_part1_mesh.enableEdgesRendering();\n'
            s += 'linkage_part1_mesh.edgesWidth = 0.3;\n'
            s += 'linkage_part1_mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 1);\n'

        if part2_parent is not None:
            s += 'var linkage_part2_mesh = BABYLON.MeshBuilder.CreateBox("prismatic part2", {depth:0.03, height:0.03, width:0.06}, scene);\n'
            s += 'linkage_part2_mesh.parent = {};\n'.format(part2_parent)
            s += 'linkage_part2_mesh.position = new BABYLON.Vector3({}, {}, {});\n'.format(*self.part2_position_function(initial_linkage_parameters))
            s += 'linkage_part2_mesh.rotation = BABYLON.Vector3.RotationFromAxis(new BABYLON.Vector3({}, {}, {}),new BABYLON.Vector3({}, {}, {}), new BABYLON.Vector3({}, {}, {}));\n'.format(*bp2.u, *bp2.v, *bp2.w)
            s += 'linkage_part2_mesh.enableEdgesRendering();\n'
            s += 'linkage_part2_mesh.edgesWidth = 0.3;\n'
            s += 'linkage_part2_mesh.edgesColor = new BABYLON.Color4(0, 0, 0, 1);\n'

        return s

class LimitedBallLinkage(Linkage):
    holonomic = True

    def __init__(self,
                 part1, part1_position, part1_basis,
                 part2, part2_position, part2_basis,
                 name='LimitedBallLinkage'):
        """
        Allowed movements are:
            - a rotation around part1 basis u
            - a rotation around part1 basis v
            
        """

        def part1_basis_f(q):
            return part1_basis.Rotation(part1_basis.u, q[0], copy=True)\
                              .Rotation(part1_basis.v, q[1], copy=True)

        def part2_basis_f(q):
            return part2_basis


        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

        Linkage.__init__(self,
                         part1, lambda q: part1_position, part1_basis_f,
                         part2, lambda q: part2_position, part2_basis_f,
                         False, True,
                         [Parameter(0., 2*math.pi, 2*math.pi),
                          Parameter(0., 2*math.pi, 2*math.pi)],
                         name)

class BallLinkage(Linkage):
    holonomic = True

    def __init__(self,
                 part1, part1_position, part1_basis,
                 part2, part2_position, part2_basis,
                 name='BallLinkage'):
        """

        """

        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

        part1_basis_f, part2_basis_f = self.basis_functions()

        Linkage.__init__(self,
                         part1, lambda q: part1_position, part1_basis_f,
                         part2, lambda q: part2_position, part2_basis_f,
                         False, True,
                         [Parameter(0., 2*math.pi, 2*math.pi),
                          Parameter(0., 2*math.pi, 2*math.pi),
                          Parameter(0., 2*math.pi, 2*math.pi)], name)


    def update_part1_point(self, new_position):
        self.part1_position = new_position
        self.part1_position_function = lambda q: new_position


    def update_part2_point(self, new_position):
        self.part2_position = new_position
        self.part2_position_function = lambda q: new_position

    def basis_functions(self):
        
        def part1_basis_f(q):
            return self.part1_basis.Rotation(self.part1_basis.u, q[0], copy=True)\
                              .Rotation(self.part1_basis.v, q[1], copy=True)\
                              .Rotation(self.part1_basis.w, q[2], copy=True)

        def part2_basis_f(q):
            return self.part2_basis
        
        return part1_basis_f, part2_basis_f

class NoConfigurationFoundError(Exception):
    pass

class MovingMechanism(Mechanism):
    def __init__(self, linkages, ground, name=''):

        Mechanism.__init__(self,
                           linkages,
                           ground,
                           {},
                           None,
                           None,
                           name=name)

        # self.parts_setting_path = {}
        self._settings_path = {}

        # Settings
        self.settings_graph()
        n_kp = 0
        self.kinematic_parameters_mapping = {}
        for linkage in self.linkages_kinematic_setting:
            for i in range(linkage.number_kinematic_parameters):
                self.kinematic_parameters_mapping[linkage, i] = n_kp + i
            n_kp += linkage.number_kinematic_parameters

    def settings_graph(self):
        graph = self.holonomic_graph.copy()
        self.opened_linkages = []
        # print('cycles deleting')
        # graph_has_cycles = True
        graph_cycles = nx.cycle_basis(graph)
        
        while len(graph_cycles) != 0:
            # Deleting first cycle of graph
            
            ground_distance = [(l, len(nx.shortest_path(graph, l, self.ground)))\
                               for l in graph_cycles[0]\
                               if l in self.linkages\
                                   and not l in self.opened_linkages\
                                   and not l.positions_require_kinematic_parameters
                               ]

            linkage_to_delete = max(ground_distance, key=lambda x:x[1])[0]
            # print(linkage_to_delete.name)
            self.opened_linkages.append(linkage_to_delete)
            graph.remove_node(linkage_to_delete)
            graph_cycles = nx.cycle_basis(graph)

        self.linkages_kinematic_setting = [l for l in self.linkages if l not in self.opened_linkages]
        self.settings_graph = graph


    def plot_settings_graph(self):
        s="""<html>
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
        index={}
        for ipart,part in enumerate(self.parts+[self.ground]):
            index[part]=ipart
            s+="{{id: {}, label: '{}'}},\n".format(ipart,part.name)
#        s+=']);\n'
        n=len(self.parts)+1
#        index[self.ground]=n
#        n+=1
        for il,linkage in enumerate(self.linkages_kinematic_setting):
            index[linkage] = n+il
            s+="{{id: {}, label: '{}'}},\n".format(n+il,linkage.name)
        s+=']);\n'

        s+="var edges = new vis.DataSet(["
        for linkage in self.linkages_kinematic_setting:
            s+='{{from: {}, to: {}}},\n'.format(index[linkage],index[linkage.part1])
            s+='{{from: {}, to: {}}},\n'.format(index[linkage],index[linkage.part2])
        s+=']);'

        s+="""
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

        with open('gm_graph_viz.html','w') as file:
            file.write(s)

        webbrowser.open('file://' + os.path.realpath('gm_graph_viz.html'))


    def settings_path(self, part1, part2):
        if (part1, part2) in self._settings_path:
            return self._settings_path[part1, part2]
        elif (part2, part1) in self._settings_path:
            path = [(p2, linkage, not linkage_side, p1) for (p1, linkage, linkage_side, p2) in self._settings_path[part2, part1][::-1]]
            self._settings_path[part1, part2] = path
            return self._settings_path[part1, part2]
        else:
            path = []
            try:
                raw_path = list(nx.shortest_path(self.settings_graph, part1, part2))
            except nx.NetworkXNoPath:
                self.plot_settings_graph()
                print(part1.name, part2.name)
                raise nx.NetworkXError
            for path_part1, linkage, path_part2 in zip(raw_path[:-2:2], raw_path[1::2], raw_path[2::2]+[part2]):
                path.append((path_part1, linkage, linkage.part1==path_part1, path_part2))

            self._settings_path[part1, part2] = path
            return path

    def part_global_frame(self, part, kinematic_parameters_values):
        frame = vm.OXYZ
        for part1, linkage, linkage_side, part2 in self.settings_path(self.ground, part):
            linkage_parameters_values = self.extract_linkage_parameters_values(linkage, kinematic_parameters_values)
            linkage_frame = linkage.frame(linkage_parameters_values, side=linkage_side)
            frame = frame + linkage_frame

        return frame


    def part_relative_frame(self, part, reference_part, kinematic_parameters_values):
        frame = vm.OXYZ
        for part1, linkage, linkage_side, part2 in self.settings_path(reference_part, part):
            linkage_parameters_values = self.extract_linkage_parameters_values(linkage, kinematic_parameters_values)
            linkage_frame = linkage.frame(linkage_parameters_values, side=linkage_side)
            frame = frame + linkage_frame

        return frame    

    def linkage_global_position(self, linkage, global_parameter_values):
        if linkage.positions_require_kinematic_parameters:
            ql = self.extract_linkage_parameters_values(linkage,
                                                        global_parameter_values)
        else:
            ql = []


        part1_frame = self.part_global_frame(linkage.part1, global_parameter_values)
        return part1_frame.OldCoordinates(linkage.part1_position_function(ql))


    def extract_linkage_parameters_values(self, linkage, global_parameter_values):
        linkage_parameters = [global_parameter_values[self.kinematic_parameters_mapping[linkage, i]]\
               for i in range(linkage.number_kinematic_parameters)]
        return linkage_parameters

    def opened_linkage_gap(self, linkage, global_parameter_values):
        if linkage.positions_require_kinematic_parameters:
            ql = self.extract_linkage_parameters_values(linkage, global_parameter_values)
        else:
            ql = []
        position1 = self.part_global_frame(linkage.part1, global_parameter_values).OldCoordinates(linkage.part1_position_function(ql))
        position2 = self.part_global_frame(linkage.part2, global_parameter_values).OldCoordinates(linkage.part2_position_function(ql))
        return position2 - position1

    def opened_linkage_misalignment(self, linkage, global_parameter_values):
        ql = self.extract_linkage_parameters_values(linkage, global_parameter_values)

        basis1 = self.part_global_frame(linkage.part1, global_parameter_values).Basis()
        basis2 = self.part_global_frame(linkage.part2, global_parameter_values).Basis()
        basis = basis2 - basis1 - linkage.basis(ql)
        return basis

    def opened_linkages_residue(self, q):
        residue = 0.
        for linkage in self.opened_linkages:
            residue += self.opened_linkage_gap(linkage, q).Norm()
        return residue

    def reduced_x_to_full_x(self, xr, basis_vector, free_parameters_dofs):
        x = basis_vector[:]

        for qrv, i in zip(xr, free_parameters_dofs):
            x[i] = qrv
        return x

    def full_x_to_reduced_x(self, x, free_parameters_dofs):
        return [x[i] for i in free_parameters_dofs]

    def geometric_closing_residue_function(self, basis_vector,
                                           free_parameters_dofs):


        def residue_function(xr):
            x = self.reduced_x_to_full_x(xr, basis_vector, free_parameters_dofs)


            return self.opened_linkages_residue(x)

        return residue_function


    def _optimization_settings(self, imposed_parameters):
        # Free parameter identification
        free_parameters_dofs = []
        free_parameters = []
        n_free_parameters = 0
        n_parameters = len(self.kinematic_parameters_mapping.items())
        basis_vector = [0] * n_parameters
        for i in range(n_parameters):
            if i in imposed_parameters:
                basis_vector[i] = imposed_parameters[i]
            else:
                free_parameters_dofs.append(i)
                n_free_parameters += 1

        bounds = []
        for (linkage, iparameter), idof in self.kinematic_parameters_mapping.items():
            if idof in free_parameters_dofs:
                parameter = linkage.kinematic_parameters[iparameter]
                free_parameters.append(parameter)
                bounds.append(parameter.optimizer_bounds())

            bounds_cma = [[], []]
            for bmin, bmax in bounds:
               bounds_cma[0].append(bmin)
               bounds_cma[1].append(bmax)

        return basis_vector, free_parameters_dofs, free_parameters, n_free_parameters, bounds, bounds_cma


    def find_configurations(self, imposed_parameters,
                           number_max_configurations,
                           number_starts=10, tol=1e-5,
                           starting_point=None):

        # initial_imposed_parameters = {k: v[0] for k,v in steps_imposed_parameters.items()}


        basis_vector, free_parameters_dofs, free_parameters, n_free_parameters, bounds, bounds_cma\
            = self._optimization_settings(imposed_parameters)

        geometric_closing_residue = self.geometric_closing_residue_function(basis_vector,
                                                                            free_parameters_dofs)


        # Starting n times
        starting_points = []

        for istart in range(number_starts):
            if starting_point is None:
                xr0 = [0]*n_free_parameters
                for i, parameter in enumerate(free_parameters):
                    xr0[i] = parameter.random_value()
            else:
                xr0 = [starting_point[i] for i in free_parameters_dofs]

            # result = minimize(geometric_closing_residue, xr0, bounds=bounds,
            #                   tol=0.1*tol)
            # fopt = result.fun
            # if fopt < tol:
            #     xr_opt = result.x
            # else:

            xr_opt, fopt = cma.fmin(geometric_closing_residue, xr0, 0.2,
                                  options={'bounds':bounds_cma,
                                           'ftarget': tol,
                                           'verbose': -9,
                                           'maxiter': 2000})[0:2]


            if fopt <= tol:
#                print('converged')
                found_x = False
                for x in starting_points:
                    equal = True
                    for parameter, xi1, xi2 in zip(free_parameters, x, xr_opt):
                        if not parameter.are_values_equal(xi1, xi2):
                            equal = False
                            break
                    if equal:
#                        print(x, result.x, 'are equal')
                        found_x = True

                if not found_x:
                    xopt = self.reduced_x_to_full_x(xr_opt, basis_vector, free_parameters_dofs)
                    starting_points.append(xopt[:])
                    yield xopt
                    if len(starting_points) >= number_max_configurations:
                        break

        print('Found {} configurations'.format(len(starting_points)))
        raise NoConfigurationFoundError


    def solve_from_initial_configuration(self, initial_parameter_values,
                                         steps_imposed_parameters,
                                         number_step_retries=5,
                                         max_failed_steps=3,
                                         tol=1e-4):
        """
        returns a MechanismConfigurations object. The initial point deduced from initial_parameter_values
        is the first step of the MechanismConfigurations object.
        """


        x0 = initial_parameter_values
        step_imposed_parameters = {k: v[0] for k, v in steps_imposed_parameters.items()}
        basis_vector, free_parameters_dofs, free_parameters, n_free_parameters, bounds, bounds_cma\
                = self._optimization_settings(step_imposed_parameters)
        xr0 = self.full_x_to_reduced_x(x0, free_parameters_dofs)

        n_steps = len(list(steps_imposed_parameters.values())[0])

        qs = [x0]

        number_failed_steps = 0
        failed_step = False
        for istep in range(n_steps):
            step_imposed_parameters = {k: v[istep] for k, v in steps_imposed_parameters.items()}

            # basis vector needs update at each time step!
            basis_vector, free_parameters_dofs, free_parameters, n_free_parameters, bounds, bounds_cma\
                = self._optimization_settings(step_imposed_parameters)

            geometric_closing_residue = self.geometric_closing_residue_function(basis_vector,
                                                                                free_parameters_dofs)

            if n_free_parameters > 0:
                step_converged = False
                n_tries_step = 1
                while (not step_converged) and (n_tries_step<= number_step_retries):
                    result = minimize(geometric_closing_residue,
                                      npy.array(xr0)+0.01*(npy.random.random(n_free_parameters)-0.5),
                                      tol=0.1*tol, bounds=bounds)
                    xr_opt = result.x
                    fopt = result.fun
                    if fopt > tol:
                        xr_opt, fopt = cma.fmin(geometric_closing_residue, xr0, 0.1,
                                              options={'bounds':bounds_cma,
    #                                              'tolfun':0.5*tol,
                                                       'verbose':-9,
                                                       'ftarget': tol,
                                                       'maxiter': 500})[0:2]
                    n_tries_step += 1
                    step_converged = fopt < tol
                if step_converged:
                    xr0 = xr_opt[:]
                    x = self.reduced_x_to_full_x(xr_opt, basis_vector, free_parameters_dofs)

                    qs.append(x[:])
                else:
                    print('@istep {}: residue: {}'.format(istep, fopt))
                    number_failed_steps += 1
                    if number_failed_steps >= max_failed_steps:
                        print('Failed {} steps, stopping configuration computation'.format(max_failed_steps))
                        failed_step = True
                        break


            else:
                qs.append(basis_vector)
        if not failed_step:
            return MechanismConfigurations(self,steps_imposed_parameters, qs)

#    def solve_configurations(steps_imposed_parameters,
#                             number_max_configurations)
#
#        for configuration in self.find_initial_configurations(steps_imposed_parameters,
#                                    number_max_configurations,
#                                    number_starts=10, tol=1e-5):
#
#            yield self.solve_from_initial_configuration(self, initial_parameter_values,
#                                         steps_imposed_parameters,
#                                         number_step_retries=5,
#                                         max_failed_steps=3,
#                                         tol=1e-4)

def istep_from_value_on_list(list_, value):
    for ipoint, (point1, point2) in enumerate(zip(list_[:-1],
                                                  list_[1:])):
        interval = sorted((point1, point2))
        if (interval[0] <= value) and (value <= interval[1]):
            alpha = (value-point1)/(point2-point1)
            if alpha < 0 or alpha > 1:
                raise ValueError
            return ipoint + alpha
    values = [p for p in list_]
    min_values = min(values)
    max_values = max(values)
    raise ValueError('Specified value not found in list_: {} not in [{}, {}]'.format(value, min_values, max_values))


def istep_from_value_on_trajectory(trajectory, value, axis):
    for ipoint, (point1, point2) in enumerate(zip(trajectory[:-1],
                                                  trajectory[1:])):
        interval = sorted((point1[axis], point2[axis]))
        if (interval[0] <= value) and (value <= interval[1]):
            alpha = (value-point1[axis])/(point2[axis]-point1[axis])
            if alpha < 0 or alpha > 1:
                raise ValueError
            return ipoint + alpha
    values = [p[axis] for p in trajectory]
    min_values = min(values)
    max_values = max(values)
    raise ValueError('Specified value not found in trajectory: {} not in [{}, {}]'.format(value, min_values, max_values))

def point_from_istep_on_trajectory(trajectory, istep):
    istep1 = int(istep)
    if istep1 == istep:
        # No interpolation needed
        return trajectory[int(istep)]
    else:
        alpha = istep - istep1
        point1 = trajectory[istep1]
        point2 = trajectory[istep1+1]
        return (1-alpha)*point1+(alpha)*point2

def trajectory_point_from_value(trajectory, value, axis):
    for ipoint, (point1, point2) in enumerate(zip(trajectory[:-1],
                                                  trajectory[1:])):
        interval = sorted((point1[axis], point2[axis]))
        if (interval[0] <= value) and (value < interval[1]):
            alpha = (value - point1[axis])/(point2[axis] - point1[axis])
            return (1-alpha)*point1 + alpha*point2
    return None

def trajectory_derivative(trajectory, istep, delta_istep):
    istep1 = istep-0.5*delta_istep
    istep2 = istep+0.5*delta_istep
    if istep1 < 0:
        istep1 = 0
        istep2 = delta_istep
    if istep2 > len(trajectory)-1:
        istep2 = len(trajectory)-1
        istep1 = istep2 - delta_istep
        if istep1 < 0:
            raise ValueError('Delta istep is too large!')
    point1 = point_from_istep_on_trajectory(trajectory, istep1)
    point2 = point_from_istep_on_trajectory(trajectory, istep2)
    return (point2-point1)

class MechanismConfigurations(DessiaObject):

    def __init__(self,
                 mechanism,
                 steps_imposed_parameters,
                 steps):

        number_steps = len(steps)


        DessiaObject.__init__(self,
                              mechanism=mechanism,
                              steps_imposed_parameters=steps_imposed_parameters,
                              steps=steps,
                              number_steps=number_steps)

        if not self.is_valid():
            raise ValueError

        self.trajectories = {}

    def is_valid(self):
        nparam_mechanism = len(self.mechanism.kinematic_parameters_mapping)
        for istep, step in enumerate(self.steps):
            if len(step) != nparam_mechanism:
                print('Step n°{} has incorrect length'.format(istep))
                return False
        return True

    def opened_linkages_residue(self):
        residues = []
        for step in self.steps:
            residues.append(self.mechanism.opened_linkages_residue(step))
        return residues

    def interpolate_step(self, istep):
        """
        :param istep: can be a float
        """

        istep1 = int(istep)
        
        alpha = istep - istep1
        if alpha == 0.:
            return self.steps[istep1]

        return [(1-alpha)*s1+alpha*s2 for s1, s2 in zip(self.steps[istep1],
                                                        self.steps[istep1+1])]
        

    def plot_kinematic_parameters(self,
                                  linkage1, kinematic_parameter1,
                                  linkage2, kinematic_parameter2
                                  ):
        x = []
        y = []
        dof1 = self.mechanism.kinematic_parameters_mapping[linkage1, kinematic_parameter1]
        dof2 = self.mechanism.kinematic_parameters_mapping[linkage2, kinematic_parameter2]
        for step in self.steps:
            x.append(step[dof1])
            y.append(step[dof2])
        fig, ax = plt.subplots()
        ax.plot(x, y, marker='o')
        ax.set_xlabel('Parameter {} of linkage {}'.format(kinematic_parameter1+1, linkage1.name))
        ax.set_ylabel('Parameter {} of linkage {}'.format(kinematic_parameter2+1, linkage2.name))
        ax.grid()
        return fig, ax

    def trajectory(self, point, part, reference_part):
        if (point, part, reference_part) in self.trajectories:
            return self.trajectories[point, part, reference_part]

        trajectory = []
        for step in self.steps:
            frame1 = self.mechanism.part_global_frame(part, step)
            frame2 = self.mechanism.part_global_frame(reference_part, step)
            frame = frame1 - frame2
            trajectory.append(frame.OldCoordinates(point))

        self.trajectories[point, part, reference_part] = trajectory

        return trajectory

    def plot2D_trajectory(self, point, part, reference_part, x=vm.X3D, y=vm.Y3D, equal_aspect=True):
        xt = []
        yt = []
        for traj_point in self.trajectory(point, part, reference_part):
            xp, yp = traj_point.PlaneProjection2D(x, y)
            xt.append(xp)
            yt.append(yp)

        fig, ax = plt.subplots()
        ax.plot(xt, yt, marker='o')
        ax.grid()
        ax.set_xlabel(str(x))
        ax.set_ylabel(str(y))
        ax.set_title('Trajectory of point {} on part {} relatively to part {}'.format(str(point), part.name, reference_part.name))

        if equal_aspect:
            ax.set_aspect('equal')

        return fig, ax

    def plot_trajectory(self, point, part, reference_part, equal_aspect=True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        xt = []
        yt = []
        zt = []
        for point in self.trajectory(point, part, reference_part):
            xp, yp, zp = point
            xt.append(xp)
            yt.append(yp)
            zt.append(zp)

#        fig, ax = plt.subplots()
        ax.plot(xt, yt, zt, marker='o')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Trajectory of point {} on part {} relatively to part {}'.format(str(point), part.name, reference_part.name))


        if equal_aspect:
            ax.set_aspect('equal')
#        fig.canvas.set_window_title('Trajectory')
        return fig, ax



    def part_local_point_global_speed(self, part, point, istep):
        """

        """
        if (istep < 0) or (istep > self.number_steps-1):
            raise ValueError('istep outside of bounds: {}'.format(istep))
        elif istep < 0.5:
            # Backward extrapolation from speeds 1 and 2
            frame1 = self.mechanism.part_global_frame(part, self.steps[0])
            frame2 = self.mechanism.part_global_frame(part, self.steps[1])
            frame3 = self.mechanism.part_global_frame(part, self.steps[2])
            p1 = frame1.OldCoordinates(point)
            p2 = frame2.OldCoordinates(point)
            p3 = frame3.OldCoordinates(point)
            v1 = p2 - p1
            v2 = p3 - p2
            alpha = istep - 0.5
            return (1-alpha)*v1 + alpha*v2

        elif istep > self.number_steps-1.5:
            # forward extrapolation from speeds n-1 and n
            i1 = int(istep-0.5)
            frame1 = self.mechanism.part_global_frame(part, self.steps[-3])
            frame2 = self.mechanism.part_global_frame(part, self.steps[-2])
            frame3 = self.mechanism.part_global_frame(part, self.steps[-1])
            p1 = frame1.OldCoordinates(point)
            p2 = frame2.OldCoordinates(point)
            p3 = frame3.OldCoordinates(point)
            v1 = p2 - p1
            v2 = p3 - p2
            alpha = istep - (self.number_steps - 2.5)
            return (1-alpha)*v1 + alpha*v2
        else:

            int_istep = int(istep)
            if int_istep+0.5 == istep:
                # Using exact derivative
                frame1 = self.mechanism.part_global_frame(part, self.steps[int_istep])
                frame2 = self.mechanism.part_global_frame(part, self.steps[int_istep+1])
                p1 = frame1.OldCoordinates(point)
                p2 = frame2.OldCoordinates(point)
                return p2 - p1
            else:
                # interpolation in between
                i1 = int(istep-0.5)
                frame1 = self.mechanism.part_global_frame(part, self.steps[i1])
                frame2 = self.mechanism.part_global_frame(part, self.steps[i1+1])
                frame3 = self.mechanism.part_global_frame(part, self.steps[i1+2])
                p1 = frame1.OldCoordinates(point)
                p2 = frame2.OldCoordinates(point)
                p3 = frame3.OldCoordinates(point)
                v1 = p2 - p1
                v2 = p3 - p2
                alpha = istep - i1 - 0.5
                return (1-alpha)*v1 + alpha*v2



    def part_global_rotation_vector(self, part, istep):

        step = self.interpolate_step(istep)
        frame = self.mechanism.part_global_frame(part, step)

        point1 = vm.O3D
        point1_speed = self.part_local_point_global_speed(part, point1, istep)
        for point2 in [vm.X3D, vm.Y3D, vm.Z3D]:

            point2_speed = self.part_local_point_global_speed(part, point2, istep)
            delta_speeds = point2_speed - point1_speed
            if not math.isclose(delta_speeds.Norm(), 0, abs_tol=1e-8):
                break

        p21 = frame.OldCoordinates(point2) - frame.OldCoordinates(point1)
        R = delta_speeds.Cross(p21)#/d_p21_2
        return R


    def part_instant_rotation_global_axis_point(self, part, istep):
        w = self.part_global_rotation_vector(part, istep)
        w2 = w.Dot(w)
        if math.isclose(w2, 0, abs_tol=1e-8):
            return None
        step = self.interpolate_step(istep)
        frame = self.mechanism.part_global_frame(part, step)

        for point in [vm.O3D, 0.1*vm.X3D, 0.1*vm.Y3D, 0.1*vm.Z3D]:
            vp = self.part_local_point_global_speed(part, point, istep)

            if not math.isclose(vp.Norm(), 0, abs_tol=1e-6):
                return frame.OldCoordinates(point) - w.Cross(vp)/w2
        raise ValueError


    def plot2D(self, x=vm.X3D, y=vm.Y3D, isteps=None, plot_frames=False,
               plot_rotation_axis=False):
        fig, ax = plt.subplots()

        # Linkage colors
        np = len(self.mechanism.parts)
        colors = {p: hsv_to_rgb((ip / np, 0.78, 0.87)) for ip, p in enumerate(self.mechanism.parts)}
        colors[self.mechanism.ground] = (0,0,0)

#            i: to_hex(
#                ) for i in range(nlines)}

        if isteps == None:
            steps = self.steps[:]
        else:
            steps = [self.steps[i] for i in isteps]

#        # Fetching wireframes lines
#        wireframes = {}
#        for part in self.mechanism.parts:
#            # Fetching local points
#            part_points = []
#            for linkage in self.mechanism.part_linkages:


        for istep, step in enumerate(steps):
            linkage_positions = {}
            part_frames = {}
            for linkage in self.mechanism.linkages:

    #            flp1.origin.PlaneProjection2D(x, y).MPLPlot(ax=ax)
                if linkage.positions_require_kinematic_parameters:
                    ql = self.mechanism.extract_linkage_parameters_values(linkage,
                                                                          step)
                else:
                    ql = []

                part1_frame = self.mechanism.part_global_frame(linkage.part1,
                                                        step)
                part_frames[linkage.part1] = part1_frame
    #
                linkage_position1 = part1_frame.OldCoordinates(linkage.part1_position_function(ql))
                linkage_position1_2D = linkage_position1.PlaneProjection2D(x, y)

                part2_frame = self.mechanism.part_global_frame(linkage.part2,
                                                        step)
                part_frames[linkage.part2] = part2_frame
    #
                linkage_position2 = part2_frame.OldCoordinates(linkage.part2_position_function(ql))
                linkage_position2_2D = linkage_position1.PlaneProjection2D(x, y)

                if linkage_position1 != linkage_position2:
                    ax.text(*linkage_position1_2D, linkage.name+' position1')
                    ax.text(*linkage_position2_2D, linkage.name+' position2')
                    error = linkage_position2_2D - linkage_position1_2D
                    ax.add_patch(Arrow(*linkage_position1_2D,
                                       *error, 0.05))
                else:
                    if istep == 0:
                        ax.text(*linkage_position1_2D, linkage.name)

                linkage_positions[linkage, linkage.part1] = linkage_position1
                linkage_positions[linkage, linkage.part2] = linkage_position2


            part_linkages = self.mechanism.part_linkages()
            del part_linkages[self.mechanism.ground]

            for ipart, (part, linkages) in enumerate(part_linkages.items()):
#                middle_point = vm.o2D
#                for linkage in linkages:
#                    middle_point += linkage_positions[linkage, part]
#                for point in part.interest_points:
#                    middle_point += point
#                middle_point /= (len(linkages) + len(part.interest_points))
#                xm, ym = middle_point.vector
                points = []
                for linkage in linkages:
                    points.append(linkage_positions[linkage, part])
                points.extend([part_frames[part].OldCoordinates(p) for p in part.interest_points])
                xm, ym = vm.Point3D.mean_point(points).PlaneProjection2D(x, y).vector

                if istep == 0:
                    ax.text(xm, ym, part.name + ' step 0',
                            ha="center", va="center",
                            bbox=dict(boxstyle="square",
                                      ec=colors[part],
                                      fc=(1., 1, 1),
                                      ))
                else:
                    if ipart == 0:
                        ax.text(xm, ym, 'step {}'.format(istep),
                                ha="center", va="center",
                                bbox=dict(boxstyle="square",
                                      ec=colors[part],
                                      fc=(1., 1, 1),
                                      ))

#                for linkage in linkages:
#                    x1, y1 = linkage_positions[linkage, part]
#                    ax.plot([x1, xm], [y1, ym], color=colors[part])

                for line in Part.wireframe_lines(points):
                    line.MPLPlot2D(x, y, ax, color=colors[part], width=5)

                part_frame = self.mechanism.part_global_frame(part, step)
                for point in part.interest_points:
                    x1, y1 = part_frame.OldCoordinates(point).PlaneProjection2D(x, y)
                    ax.plot([x1, xm], [y1, ym], color=colors[part])


                if plot_frames:
                    part_frame = self.mechanism.part_global_frame(part, step)
                    part_frame.plot2d(x=x, y=y, ax=ax)

                if plot_rotation_axis:
                    axis = self.part_global_rotation_vector(part, istep)
                    point = self.part_instant_rotation_global_axis_point(part, istep)
                    if point is not None:
                        axis.Normalize()
                        line = vm.Line3D(point-axis, point+axis)
                        line.PlaneProjection2D(x, y).MPLPlot(ax=ax, color=colors[part], dashed=True)


        ax.set_aspect('equal')
        ax.set_xlabel(str(x))
        ax.set_ylabel(str(y))
        ax.margins(.1)

    def babylonjs(self, page='gm_babylonjs', plot_frames=False,
                  plot_trajectories=True, plot_instant_rotation_axis=False,
                  use_cdn=False):


        page+='.html'

        env = Environment(loader=PackageLoader('genmechanics', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon.html')

        np = len(self.mechanism.parts)
        colors = {p: hsv_to_rgb((ip / np, 0.78, 0.87)) for ip, p in enumerate(self.mechanism.parts)}

        part_points = {p: [] for p in self.mechanism.parts}
        part_points[self.mechanism.ground] = []
#        part_frames = {}

        for part, linkages in self.mechanism.part_linkages().items():

            for linkage in linkages:
                if linkage.positions_require_kinematic_parameters:
                    ql = self.mechanism.extract_linkage_parameters_values(linkage,
                                                                          self.steps[0])
                else:
                    ql = []

                if part == linkage.part1:
                    linkage_position = linkage.part1_position_function(ql)
                else:
                    linkage_position = linkage.part2_position_function(ql)

                part_points[part].append(linkage_position)

            for point in part.interest_points:
                part_points[part].append(point)


        meshes_string = 'var parts_parent = [];\n'



        for part in self.mechanism.parts:
            meshes_string += 'var part_children = [];\n'
            lines = part.wireframe_lines(part_points[part])
            meshes_string += lines[0].Babylon(name='part_parent', color=colors[part])
            meshes_string += 'parts_parent.push(part_parent);\n'
            for l in lines[1:]:
                meshes_string += l.Babylon(color=colors[part], parent='part_parent')
#                meshes_string += 'part_meshes.push(line);\n'

#            # Adding interest points
#            for point in part.interest_points:
#                meshes_string += 'var point = BABYLON.MeshBuilder.CreateSphere("interest_point", {diameter: 0.01}, scene);\n'
#                meshes_string += 'point.position = new BABYLON.Vector3({}, {}, {});'.format(*point.vector)
#                meshes_string += 'part_meshes.push(point);'


            if plot_frames:
                meshes_string += vm.OXYZ.babylonjs(parent='part_parent', size=0.1)

            # if plot_instant_rotation_axis:
            #     rotation_axis = self.par

        if plot_instant_rotation_axis:
            for part in self.mechanism.parts:
                line = vm.LineSegment3D(-0.5*vm.X3D, 0.5*vm.X3D)
                meshes_string += line.Babylon(name='rotation_axis',  color=colors[part], type_='dashed')
                meshes_string += 'parts_parent.push(rotation_axis);\n'


        linkages_string = ''
        for linkage in self.mechanism.linkages:
            if linkage not in self.mechanism.opened_linkages:

                ql = self.mechanism.extract_linkage_parameters_values(linkage,
                                                                      self.steps[0])
            else:
                ql = []

            if linkage.part1 in self.mechanism.parts:
                part1_parent = 'parts_parent[{}]'.format(self.mechanism.parts.index(linkage.part1))
            else:
                part1_parent = None

            if linkage.part2 in self.mechanism.parts:
                part2_parent = 'parts_parent[{}]'.format(self.mechanism.parts.index(linkage.part2))
            else:
                part2_parent = None

            linkages_string += linkage.babylonjs(ql,
                                                 part1_parent=part1_parent,
                                                 part2_parent=part2_parent)


        # Computing positions and orientations
        positions = []
        orientations = []
        linkage_positions = []

        # n_steps = len(self.steps)

        for istep, step in enumerate(self.steps):
            step_positions = []
            step_orientations = []
            step_linkage_positions = []
#            step_linkage_positions = []
            for part in self.mechanism.parts:

                frame = round(self.mechanism.part_global_frame(part,
                                                  step))

                step_positions.append(list(frame.origin))
                step_orientations.append([list(frame.u),
                                          list(frame.v),
                                          list(frame.w)])


            if plot_instant_rotation_axis:
                # if istep == n_steps-1:
                #     istep_speed = n_steps-2
                #     # step_speed = self.steps[istep_speed]
                # else:
                #     istep_speed = istep
                #     # step_speed = step
                for part in self.mechanism.parts:

                    # print('unorm', u, u.Norm())
                    axis_point = self.part_instant_rotation_global_axis_point(part, istep)
                    if axis_point is None:
                        u = vm.X3D.copy()
                        v = vm.Y3D.copy()
                        w = vm.Z3D.copy()
                        axis_point = vm.Point3D((100, 100, 100))
                    else:
                        u = self.part_global_rotation_vector(part, istep)
                        u.Normalize()
                        v = u.RandomUnitNormalVector()
                        w = u.Cross(v)

                    step_positions.append(list(axis_point))
                    step_orientations.append([list(u),
                                              list(v),
                                              list(w)])
                    # print(step_positions)
                    # print(step_orientations)

            for linkage in self.mechanism.linkages:
                step_linkage_positions.append(list(self.mechanism.linkage_global_position(linkage, step)))

            positions.append(step_positions)
            orientations.append(step_orientations)
            linkage_positions.append(step_linkage_positions)




        trajectories = []
        if plot_trajectories:
            for trajectory in self.trajectories.values():
                trajectories.append([list(p) for p in trajectory])
#            for point, part in self.mechanism.parts:
#                for point in part.interest_points:
#                    trajectories.append([list(p.vector) for p in self.trajectory(point, part, self.mechanism.ground)])


        script = template.render(center=(0, 0, 0),
                                 length=2*0.5,
                                 meshes_string=meshes_string,
                                 linkages_string=linkages_string,
                                 positions=positions,
                                 orientations=orientations,
                                 linkage_positions=linkage_positions,
                                 trajectories=trajectories,
                                 use_cdn=use_cdn)

        with open(page,'w') as file:
            file.write(script)

        webbrowser.open('file://' + os.path.realpath(page))