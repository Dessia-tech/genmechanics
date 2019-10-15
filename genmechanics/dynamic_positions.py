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
                 part1, part1_position,
                 part2, part2_position,
                 positions_require_kinematic_parameters,
                 basis,
                 number_kinematic_parameters,
                 name=''):
        """

        """
        DessiaObject.__init__(self,
                              part1=part1,
                              part1_position=part1_position,
                              part2=part2,
                              part2_position=part2_position,
                              positions_require_kinematic_parameters=positions_require_kinematic_parameters,
                              basis=basis,
                              number_kinematic_parameters=number_kinematic_parameters,
                              name=name)


class RevoluteLinkage(Linkage):
    holonomic = True

    def __init__(self, part1, part1_position, part1_axis, part2, part2_position, name=''):
        """

        """

        def basis(q):
            return vm.xyz.Rotation(part1_axis, q[0])

        Linkage.__init__(self,
                         part1, lambda q: part1_position,
                         part2, lambda q: part2_position,
                         False,
                         basis, 1, name)


class SlidingRevoluteLinkage(Linkage):
    holonomic = True

    def __init__(self, part1, part1_initial_position, axis, part2, part2_position, name=''):
        """
        sliding on part1, fixed on part2
        """

        def basis(q):
            return vm.xyz.Rotation(axis, q[0])

        Linkage.__init__(self,
                         part1, lambda q: part1_initial_position+q[1]*axis,
                         part2, lambda q: part2_position,
                         True,
                         basis, 2, name)

class PrismaticLinkage(Linkage):
    holonomic = True

    def __init__(self, part1, part1_initial_position, axis, part2, part2_position, name=''):
        """
        sliding on part1, fixed on part2
        """

        def basis(q):
            return vm.xyz

        Linkage.__init__(self,
                         part1, lambda q: part1_initial_position+q[0]*axis,
                         part2, lambda q: part2_position,
                         True,
                         basis, 1, name)

class BallLinkage(Linkage):
    holonomic = True

    def __init__(self, part1, part1_position, part2, part2_position, name=''):
        """
        sliding on part1, fixed on part2
        """

        def basis(q):
            b = vm.xyz.Rotation(vm.x3D, q[0])
            b.Rotation(vm.y3D, q[1], False)
            b.Rotation(vm.z3D, q[2], False)
            return b

        Linkage.__init__(self,
                         part1, lambda q: part1_position,
                         part2, lambda q: part2_position,
                         False,
                         basis, 3, name)

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

            if linkage_side:
                linkage_basis = linkage.basis(linkage_parameters_values)
                origin = (linkage.part1_position(linkage_parameters_values) \
                          - linkage_basis.OldCoordinates(linkage.part2_position(linkage_parameters_values)))
            else:
                linkage_basis = -linkage.basis(linkage_parameters_values)
                origin = (linkage.part2_position(linkage_parameters_values) \
                          - linkage_basis.OldCoordinates(linkage.part1_position(linkage_parameters_values)))

#            print('linkage basis: ', linkage_basis)
            linkage_frame = vm.Frame3D(origin,
                                       linkage_basis.u,
                                       linkage_basis.v,
                                       linkage_basis.w,
                                       )


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

    def extract_linkage_parameters_values(self, linkage, global_parameter_values):

        linkage_parameters = [global_parameter_values[self.kinematic_parameters_mapping[linkage, i]]\
               for i in range(linkage.number_kinematic_parameters)]
        return linkage_parameters

    def opened_linkage_gap(self, linkage, global_parameter_values):
        if linkage.positions_require_kinematic_parameters:
            ql = self.extract_linkage_parameters_values(linkage, global_parameter_values)
        else:
            ql = []
        position1 = self.part_frame(linkage.part1, global_parameter_values).OldCoordinates(linkage.part1_position(ql))
        position2 = self.part_frame(linkage.part2, global_parameter_values).OldCoordinates(linkage.part2_position(ql))
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

        return MechanismConfigurations(self, qs)


class MechanismConfigurations(DessiaObject):

    def __init__(self, mechanism, steps):

        DessiaObject.__init__(self,
                              mechanism=mechanism,
                              steps=steps)

    def opened_linkages_residue(self):
        return self.mechanism.opened_linkages_residue(self.kinematic_parameters_values)

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
        trajectory = []
        for step in self.steps:
            frame1 = self.mechanism.part_frame(part, step)
            frame2 = self.mechanism.part_frame(reference_part, step)
            frame = frame1 - frame2
            trajectory.append(frame.OldCoordinates(point))
        return trajectory

    def plot2D_trajectory(self, point, part, reference_part, x=vm.x3D, y=vm.y3D):
        xt = []
        yt = []
        for point in self.trajectory(point, part, reference_part):
            xp, yp = point.PlaneProjection2D(x, y)
            xt.append(xp)
            yt.append(yp)

        fig, ax = plt.subplots()
        ax.plot(xt, yt, marker='o')
        ax.grid()
        ax.set_xlabel(str(x))
        ax.set_ylabel(str(y))
        ax.set_title('Trajectory of point {} on part {} relatively to part {}'.format(str(point), part.name, reference_part.name))

        return fig, ax

    def plot_trajectory(self, point, part, reference_part):
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

#        ax.set_aspect('equal')
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
            for linkage in self.mechanism.linkages:

    #            flp1.origin.PlaneProjection2D(x, y).MPLPlot(ax=ax)
                if linkage.positions_require_kinematic_parameters:
                    ql = self.mechanism.extract_linkage_parameters_values(linkage,
                                                                          step)
                else:
                    ql = []

                part1_frame = self.mechanism.part_frame(linkage.part1,
                                                        step)
    #
                linkage_position1 = part1_frame.OldCoordinates(linkage.part1_position(ql))
                linkage_position1_2D = linkage_position1.PlaneProjection2D(x, y)

                part2_frame = self.mechanism.part_frame(linkage.part2,
                                                        step)
    #
                linkage_position2 = part2_frame.OldCoordinates(linkage.part2_position(ql))
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
                points.extend(part.interest_points)
#                print(points)
                xm, ym, _ = vm.Point3D.mean_point(points).vector

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

    def babylonjs(self, page='gm_babylonjs', plot_frames=False, plot_trajectories=False):
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
#            part_frame = self.mechanism.part_frame(part,
#                                                   self.steps[0])
            for linkage in linkages:
                if linkage.positions_require_kinematic_parameters:
                    ql = self.mechanism.extract_linkage_parameters_values(linkage,
                                                                          self.steps[0])
                else:
                    ql = []


                if part == linkage.part1:
                    linkage_position = linkage.part1_position(ql)
                else:
                    linkage_position = linkage.part2_position(ql)

                part_points[part].append(linkage_position)

            for point in part.interest_points:
                part_points[part].append(point)


        meshes_string = 'var parts = [];\n'

        for part in self.mechanism.parts:
            meshes_string += 'var part_meshes = [];\n'
            for l in part.wireframe_lines(part_points[part]):
                meshes_string += l.Babylon(color=colors[part])
                meshes_string += 'part_meshes.push(line);\n'


            if plot_frames:
                meshes_string += vm.oxyz.babylonjs()
                meshes_string += 'part_meshes.push(line1);\n'
                meshes_string += 'part_meshes.push(line2);\n'
                meshes_string += 'part_meshes.push(line3);\n'

            meshes_string += 'parts.push(part_meshes);\n'
            
            
            
        # Computing positions and orientations
        positions = []
        orientations = []
#        frame_babylon = vm.Frame3D(vm.o3D, vm.x3D, vm.z3D, vm.y3D)
        for step in self.steps:
            step_positions = []
            step_orientations = []
            for part in self.mechanism.parts:

                frame = round(self.mechanism.part_frame(part,
                                                  step))

                step_positions.append(list(frame.origin))
                step_orientations.append([list(frame.u),
                                          list(frame.v),
                                          list(frame.w)])


            positions.append(step_positions)
            orientations.append(step_orientations)
            
        trajectories = []
        if plot_trajectories:
            for part in self.mechanism.parts:
                for point in part.interest_points:
                    trajectories.append([list(p.vector) for p in self.trajectory(point, part, self.mechanism.ground)])

        script = template.render(center=(0, 0, 0),
                                 length=2*0.5,
                                 meshes_string=meshes_string,
                                 positions=positions,
                                 orientations=orientations,
                                 trajectories=trajectories)

        with open(page,'w') as file:
            file.write(script)

        webbrowser.open('file://' + os.path.realpath(page))