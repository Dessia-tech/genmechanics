#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""
import os
import webbrowser

import numpy as npy

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



class Linkage(DessiaObject):
    def __init__(self,
                 part1, part1_position_function, part1_basis_function,
                 part2, part2_position_function, part2_basis_function,
                 positions_require_kinematic_parameters,
                 basis_require_kinematic_parameters,
                 number_kinematic_parameters,
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
                              number_kinematic_parameters=number_kinematic_parameters,
                              name=name)
        
#    def part1_frame(self, linkage_parameters_values):
#        return self.part1_basis(linkage_parameters_values)\
#            .to_frame(self.part1_position(linkage_parameters_values))
#
#    def part2_frame(self, linkage_parameters_values):
#        return self.part2_basis(linkage_parameters_values)\
#            .to_frame(self.part2_position(linkage_parameters_values))


    def frame(self, linkage_parameters_values, side):
#        print('linkage_side', side)
        if side:     
#            print('\n====')
#            print(self.part1_basis_function(linkage_parameters_values))
#            print(self.part2_basis_function(linkage_parameters_values))
#            
#            #                linkage_basis = linkage.basis(linkage_parameters_values)
#            linkage_basis = (self.part1_basis_function(linkage_parameters_values)
#                             + self.part2_basis_function(linkage_parameters_values))
#
#            print('= ', linkage_basis)
#
#            
#            origin = (self.part1_position_function(linkage_parameters_values) \
#                      - linkage_basis.OldCoordinates(self.part2_position_function(linkage_parameters_values)))
#            
#            return linkage_basis.to_frame(origin)
            
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
                 part2, part2_position, part2_basis, name=''):
        """
        :param part2_basis: a basis defining orientation of linkage on part2

        """

        def part1_basis_f(q):
#            print('q', q)
#            return part1_basis
            return part1_basis.Rotation(part1_basis.u, q[0], copy=True)

        def part2_basis_f(q):
#            return part2_basis.Rotation(part2_basis.u, q[0], copy=True)
            return part2_basis


        Linkage.__init__(self,
                         part1, lambda q: part1_position, part1_basis_f,
                         part2, lambda q: part2_position, part2_basis_f,
                         False, True,
                         1, name)
        
        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)
        


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
                 part2, part2_position, part2_basis, name=''):
        """
        :param part2_basis: a basis defining orientation of linkage on part2

        """

        def part1_position_f(q):
            return part1_position + q[0]*part1_basis.u

        def part2_position_f(q):
            return part2_position

        def part1_basis_f(q):
            return part1_basis.Rotation(part1_basis.u, q[0], copy=True)

        Linkage.__init__(self,
                         part1, part1_position_f, part1_basis_f,
                         part2, part2_position_f, lambda q: part2_basis,
                         True, True,
                         2, name)
        
        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

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
                 part2, part2_position, part2_basis, name=''):
        """
        :param part2_basis: a basis defining orientation of linkage on part2

        """

        def part1_position_f(q):
            return part1_position + q[0]*part1_basis.u

        def part2_position_f(q):
            return part2_position


        Linkage.__init__(self,
                         part1, part1_position_f, lambda q: part1_basis,
                         part2, part2_position_f, lambda q: part2_basis,
                         True, False,
                         1, name)
        
        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)
        

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


class BallLinkage(Linkage):
    holonomic = True
    
    
#        def __init__(self,
#                 part1, part1_position, part1_basis,
#                 part2, part2_position, part2_basis, name=''):
#        """
#        :param part2_basis: a basis defining orientation of linkage on part2
#
#        """
#
#        def part1_basis_f(q):
##            print('q', q)
##            return part1_basis
#            return part1_basis.Rotation(part1_basis.u, q[0], copy=True)
#
#        def part2_basis_f(q):
##            return part2_basis.Rotation(part2_basis.u, q[0], copy=True)
#            return part2_basis
#
#
#        Linkage.__init__(self,
#                         part1, lambda q: part1_position, part1_basis_f,
#                         part2, lambda q: part2_position, part2_basis_f,
#                         False, True,
#                         1, name)
#        
#        DessiaObject.__init__(self,
#                              part1_position=part1_position,
#                              part2_position=part2_position,
#                              part1_basis=part1_basis,
#                              part2_basis=part2_basis)
        

    def __init__(self,
                 part1, part1_position, part1_basis,
                 part2, part2_position, part2_basis,
                 name=''):
        """

        """

        def part1_basis_f(q):
#            print('q', q)
#            return part1_basis
            return part1_basis.Rotation(part1_basis.u, q[0], copy=True)\
                              .Rotation(part1_basis.v, q[1], copy=True)\
                              .Rotation(part1_basis.w, q[2], copy=True)

        def part2_basis_f(q):
#            return part2_basis.Rotation(part2_basis.u, q[0], copy=True)
            return part2_basis


        Linkage.__init__(self,
                         part1, lambda q: part1_position, part1_basis_f,
                         part2, lambda q: part2_position, part2_basis_f,
                         False, True,
                         3, name)
        
        DessiaObject.__init__(self,
                              part1_position=part1_position,
                              part2_position=part2_position,
                              part1_basis=part1_basis,
                              part2_basis=part2_basis)

#class PartWireFrame(Part):
#    def __init__(self, part):
#        Part.__init__(self, name=part.name, interest_points=part.interest_points)



class MovingMechanism(Mechanism):
    def __init__(self, linkages, ground, name):

        Mechanism.__init__(self,
                           linkages,
                           ground,
                           {},
                           None,
                           None,
                           name='')

        self.parts_setting_path = {}
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
        for cycle in nx.cycle_basis(graph):
            ground_distance = [(l, len(nx.shortest_path(graph, l, self.ground)))\
                               for l in cycle\
                               if l in self.linkages\
                                   and not l in self.opened_linkages\
                                   and not l.positions_require_kinematic_parameters
                               ]
#            print(ground_distance)
            linkage_to_delete = max(ground_distance, key=lambda x:x[1])[0]
#            print(linkage_to_delete)
            self.opened_linkages.append(linkage_to_delete)
            graph.remove_node(linkage_to_delete)

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
            raw_path = list(nx.shortest_path(self.settings_graph, part1, part2))
#            print('rp', raw_path)
            for part1, linkage, part2 in zip(raw_path[:-2:2], raw_path[1::2], raw_path[2::2]+[part2]):
                path.append((part1, linkage, linkage.part1==part1, part2))

            self._settings_path[part1, part2] = path
            return path

    def part_frame(self, part, kinematic_parameters_values):
        frame = vm.oxyz
        for part1, linkage, linkage_side, part2 in self.settings_path(self.ground, part):
#            print('\nj', part1.name, linkage.name, linkage_side, part2.name)
#            linkage_parameters_values = [kinematic_parameters_values[self.kinematic_parameters_mapping[linkage, i]]\
#                                         for i in range(linkage.number_kinematic_parameters)]
            linkage_parameters_values = self.extract_linkage_parameters_values(linkage, kinematic_parameters_values)
#            print('lpv', linkage_parameters_values)

#            if linkage_side:
#                linkage_basis = linkage.basis(linkage_parameters_values)
#                origin = (linkage.part1_position(linkage_parameters_values) \
#                          - linkage_basis.OldCoordinates(linkage.part2_position(linkage_parameters_values)))
#            else:
#                linkage_basis = -linkage.basis(linkage_parameters_values)
#                origin = (linkage.part2_position(linkage_parameters_values) \
#                          - linkage_basis.OldCoordinates(linkage.part1_position(linkage_parameters_values)))
#
##            print('linkage basis: ', linkage_basis)
#            linkage_frame = vm.Frame3D(origin,
#                                       linkage_basis.u,
#                                       linkage_basis.v,
#                                       linkage_basis.w,
#                                       )
            linkage_frame = linkage.frame(linkage_parameters_values, side=linkage_side)


#                origin = linkage.part2_position(linkage_parameters_values) - linkage.part1_position(linkage_parameters_values)
#                linkage_frame = vm.Frame3D(origin,
#                                           linkage_basis.u,
#                                           linkage_basis.v,
#                                           linkage_basis.w,
#                                           )

#            print('linkage_frame', linkage_frame)
#            print('frame', frame)
#            frame += linkage_frame
            frame = frame + linkage_frame
#            print('frame after', frame)


        return frame

    def linkage_global_position(self, linkage, global_parameter_values):
        if linkage.positions_require_kinematic_parameters:
            ql = self.extract_linkage_parameters_values(linkage,
                                                        global_parameter_values)
        else:
            ql = []
            
        
        part1_frame = self.part_frame(linkage.part1, global_parameter_values)
        return part1_frame.OldCoordinates(linkage.part1_position_function(ql))


    def extract_linkage_parameters_values(self, linkage, global_parameter_values):
#        if linkage in self.opened_linkages:# TODO: be sure of this!
#            return []
        linkage_parameters = [global_parameter_values[self.kinematic_parameters_mapping[linkage, i]]\
               for i in range(linkage.number_kinematic_parameters)]
        return linkage_parameters

    def opened_linkage_gap(self, linkage, global_parameter_values):
        if linkage.positions_require_kinematic_parameters:
            ql = self.extract_linkage_parameters_values(linkage, global_parameter_values)
        else:
            ql = []
        position1 = self.part_frame(linkage.part1, global_parameter_values).OldCoordinates(linkage.part1_position_function(ql))
        position2 = self.part_frame(linkage.part2, global_parameter_values).OldCoordinates(linkage.part2_position_function(ql))
        return position2 - position1

    def opened_linkage_misalignment(self, linkage, global_parameter_values):
        ql = self.extract_linkage_parameters_values(linkage, global_parameter_values)

        basis1 = self.part_frame(linkage.part1, global_parameter_values).Basis()
        basis2 = self.part_frame(linkage.part2, global_parameter_values).Basis()
        basis = basis2 - basis1 - linkage.basis(ql)# !!!!!
        return basis

    def opened_linkages_residue(self, q):
        residue = 0.
        for linkage in self.opened_linkages:
            residue += self.opened_linkage_gap(linkage, q).Norm()
#            misalignment_basis = self.opened_linkage_misalignment(linkage, q)
##            print('mb', misalignment_basis)
#            residue += (misalignment_basis.u - vm.x3D).Norm()
#            residue += (misalignment_basis.v - vm.y3D).Norm()
#            residue += (misalignment_basis.w - vm.z3D).Norm()

#        print(residue)
#            print(misalignment_basis.u.Norm()+misalignment_basis.v.Norm()+misalignment_basis.w.Norm())

        return residue

    def solve_configurations(self, imposed_parameters):

        n_parameters = len(self.kinematic_parameters_mapping.items())
        n_steps = len(list(imposed_parameters.values())[0])


        def geometric_closing_residue(qr):
            q = basis_vector[:]
            for qrv, i in zip(qr, free_parameters):
                q[i] = qrv

            return self.opened_linkages_residue(q)

        # Free parameter identification
        free_parameters = []
        for i in range(n_parameters):
            if not i in imposed_parameters:
                free_parameters.append(i)


        qs = []
#        x0 = zeros(len(free_parameters))
        x0 = npy.random.random(len(free_parameters))

        for istep in range(n_steps):

            basis_vector = [0] * n_parameters
            for i in range(n_parameters):
                if i in imposed_parameters:
                    basis_vector[i] = imposed_parameters[i][istep]

            if len(free_parameters) > 0:
                result = minimize(geometric_closing_residue, x0)
                if result.fun < 1e-5:
                    x0 = result.x
                    q = basis_vector[:]
                    for qrv, i in zip(result.x, free_parameters):
                        q[i] = qrv
                    qs.append(q)
                else:
                    print('@istep {}: residue: {}'.format(istep, result.fun))
                    print(result)
            
            else:
                qs.append(basis_vector)

        return MechanismConfigurations(self, qs)


class MechanismConfigurations(DessiaObject):

    def __init__(self, mechanism, steps):

        DessiaObject.__init__(self,
                              mechanism=mechanism,
                              steps=steps)
        
        self.trajectories = {}

    def opened_linkages_residue(self):
        return self.mechanism.opened_linkages_residue(self.kinematic_parameters_values)
    
    def interpolate_step(self, istep):
        """
        :param istep: can be a float
        """
        
        istep1 = int(istep)
        alpha = istep - istep1
        
        
        new_step = [(1-alpha)*s1+alpha*s2 for s1, s2 in zip(self.steps[istep1],
                                                            self.steps[istep1+1])]
        return new_step

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
            frame1 = self.mechanism.part_frame(part, step)
            frame2 = self.mechanism.part_frame(reference_part, step)
            frame = frame1 - frame2
            trajectory.append(frame.OldCoordinates(point))
        
        self.trajectories[point, part, reference_part] = trajectory
        
        return trajectory

    def plot2D_trajectory(self, point, part, reference_part, x=vm.x3D, y=vm.y3D, equal_aspect=True):
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

    def plot2D(self, x=vm.x3D, y=vm.y3D, isteps=None, plot_frames=False):
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

                part1_frame = self.mechanism.part_frame(linkage.part1,
                                                        step)
                part_frames[linkage.part1] = part1_frame
    #
                linkage_position1 = part1_frame.OldCoordinates(linkage.part1_position_function(ql))
                linkage_position1_2D = linkage_position1.PlaneProjection2D(x, y)

                part2_frame = self.mechanism.part_frame(linkage.part2,
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
#                print(points)
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
                    line.MPLPlot2D(x, y, ax, color=colors[part])

                part_frame = self.mechanism.part_frame(part, step)
                for point in part.interest_points:
                    x1, y1 = part_frame.OldCoordinates(point).PlaneProjection2D(x, y)
                    ax.plot([x1, xm], [y1, ym], color=colors[part])


                if plot_frames:
                    part_frame = self.mechanism.part_frame(part, step)
                    part_frame.plot2d(x=x, y=y, ax=ax)
             


        ax.set_aspect('equal')
        ax.set_xlabel(str(x))
        ax.set_ylabel(str(y))
        ax.margins(.1)

    def babylonjs(self, page='gm_babylonjs', plot_frames=False,
                  plot_trajectories=True, use_cdn=False):
        
        
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
                meshes_string += vm.oxyz.babylonjs(parent='part_parent', size=0.15)
#                meshes_string += 'part_meshes.push(line1);\n'
#                meshes_string += 'part_meshes.push(line2);\n'
#                meshes_string += 'part_meshes.push(line3);\n'

#            meshes_string += 'parts.push(part_meshes);\n'
            
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
#        frame_babylon = vm.Frame3D(vm.o3D, vm.x3D, vm.z3D, vm.y3D)
        for step in self.steps:
            step_positions = []
            step_orientations = []
            step_linkage_positions = []
#            step_linkage_positions = []
            for part in self.mechanism.parts:

                frame = round(self.mechanism.part_frame(part,
                                                  step))

                step_positions.append(list(frame.origin))
                step_orientations.append([list(frame.u),
                                          list(frame.v),
                                          list(frame.w)])

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
                    
#        print(linkage_positions)

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