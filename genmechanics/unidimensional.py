#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Springs models

"""
import numpy as npy
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import networkx as nx
import random

from matplotlib.patches import Arrow, Polygon

class Body:
    """
    :param y_position: The y position in plots of the body. 
    It has no effect on the computation
    """
    def __init__(self, initial_position, y_position=0., name=''):
        self.initial_position = initial_position
        self.y_position = y_position
        self.name = name
        
class Load:
    def __init__(self, body, value):
        self.body = body
        self.value = value
        
    def Plot(self, bodies_position=None, ax=None, width=0.008, intensity_factor=1e-4):
        if bodies_position is None:
            bodies_position = {self.body: self.body.initial_position}
        x = bodies_position[self.body]
        y = self.body.y_position
        ax.add_patch(Arrow(x, y,
                     self.value*intensity_factor, 0, width))
        ax.text(x+0.5*self.value*intensity_factor, y+0.2*width,
                '{} N'.format(round(self.value)),
                color='k', horizontalalignment='center', verticalalignment='bottom')
        
class ImposedDisplacement:
    def __init__(self, body, value):
        self.body = body
        self.value = value
        self.neqs = 1
        
    def Plot(self, bodies_positions=None, strains=None, ax=None, 
             width=0.01, color='k', intensity_factor=1e-6):
        if bodies_positions is None:
            bodies_positions = {self.body: self.body}
        
        x0 = self.body.initial_position
        x1 = bodies_positions[self.body]
        y = self.body.y_position

        h = 0.5*width
        ax.add_patch(Polygon([(x1, y), (x1-h, y-h), (x1+h, y-h)], True, color='grey'))
        if self.value != 0:
            ax.add_patch(Arrow(x0, y, self.value, 0, width = 0.5*h, color='g'))        
            ax.text(x0 + 0.5*self.value, y+0.2*h,
                    '{} mm'.format(round(1000*self.value)),
                    horizontalalignment='center', verticalalignment='bottom')
            
            
class Linkage:
    def __init__(self, body1, body2, n_constraints, name=''):
        self.body1 = body1
        self.body2 = body2
        self.n_constraints = n_constraints
        self.name = name 
        
class Spring(Linkage):
    def __init__(self, body1, body2, stiffness, free_length, name=''):
        Linkage.__init__(self, body1, body2, 0, name) 
        self.stiffness = stiffness
        self.free_length = free_length
        
    def StiffnessMatrices(self):
        return (self.stiffness * npy.array([[1,-1], [-1, 1]]),
                self.stiffness * self.free_length * npy.array([-1, 1]))
        
    def ConstraintsMatrices(self):
                return None, None
        
    def Plot(self, bodies_positions=None, strains=None, ax=None,
             linkages_width=0.01, color='k', intensity_factor=1e-6):
        if bodies_positions is None:
            bodies_positions = {self.body1: self.body1.initial_position,
                                self.body2: self.body2.initial_position}
            
        x1 = bodies_positions[self.body1]
        x2 = bodies_positions[self.body2]
        y1 = self.body1.y_position
        y2 = self.body2.y_position
        ax.plot([x1, x2], [y1, y2], linestyle='--', color=color)
        
        xm = 0.5*(x1 + x2)
        ym = 0.5*(y1 + y2)
        
        ax.text(xm, ym-0.5*linkages_width,
                '{} #{}'.format(self.__class__.__name__, str(self.__hash__())[-5:]),
                horizontalalignment='center')

        
    def Strains(self, positions):
        return self.stiffness*(self.free_length - (positions[1] - positions[0]))
        

class CompressionSpring(Spring):
    def __init__(self, body1, body2, stiffness, free_length, name=''):
        Linkage.__init__(self, body1, body2, 0, name) 
        self.stiffness = stiffness
        self.free_length = free_length
        
    def StiffnessMatrices(self):
        return (self.stiffness * npy.array([[1,-1], [-1, 1]]),
                self.stiffness * self.free_length * npy.array([-1, 1]))
        
    def ConstraintsMatrices(self):
        return None, None
        
    def ActivationCondition(self, positions, reactions):
        return (positions[1] - positions[0]) < self.free_length + 1e-10

    def DesactivationCondition(self, positions, reactions):
        return (positions[1] - positions[0]) > self.free_length + 1e-10

class UnilateralContact(Linkage):
    def __init__(self, body_left, body_right, contact_distance, name=''):
        Linkage.__init__(self, body_left, body_right, 1, name) 
        self.contact_distance = contact_distance
    
    def StiffnessMatrices(self):
        return None, None
    
    def ConstraintsMatrices(self):
        return npy.array([[1.,-1.]]), npy.array([0.])
        
    def ActivationCondition(self, positions, reactions):
#        print('act', self.name, p, (p[1] - p[0]) < self.contact_distance)
        return (positions[1] - positions[0]) <= self.contact_distance-1e-8

    def DesactivationCondition(self, positions, reactions):        
#        print('dis', self.name, r, r < -1e-10)
        return reactions[0] < 1e-10
        
    def Plot(self, bodies_positions=None, strains=None, ax=None,
             linkages_width=0.01, color='k', intensity_factor = 1e-6):
        x1 = bodies_positions[self.body1]
        x2 = bodies_positions[self.body2]
        y1 = self.body1.y_position
        y2 = self.body2.y_position
        ym = 0.5*(y1+y2)
#        h = max(0.6*abs(y1-y2), linkages_width)
#        l = ((x2-x1)**2+(y2-y1)**2)**0.5
        h = 0.5*linkages_width
        
        if self.contact_distance > 0:
            xplate1 = x1 + 0.5*(self.contact_distance)
            xplate2 = x2 - 0.5*(self.contact_distance)
            xt1 = [x1, xplate1, xplate1, xplate1]
            yt1 = [y1, ym, ym+h, ym-h]
            xt2 = [x2, xplate2, xplate2, xplate2]
            yt2 = [y2, ym, ym+h, ym-h]
        else:
            xplate1 = x1 + 0.5*(self.contact_distance)
            xplate2 = x2 - 0.5*(self.contact_distance)
            xturn1 = x1 + 0.55*(self.contact_distance)
            xturn2 = x2 - 0.55*(self.contact_distance)
            xt1 = [x1, xturn1, xturn1, xplate1, xplate1, xplate1]
            yt1 = [y1, y1, ym, ym, ym+h, ym-h]
            xt2 = [x2, xturn2, xturn2, xplate2, xplate2, xplate2]
            yt2 = [y2, y2, ym, ym, ym+h, ym-h]
        
            
        ax.plot(xt1, yt1, color=color)
        ax.plot(xt2, yt2, color=color)
        ax.text(0.5*(x1+x2), ym-0.5*linkages_width,
                'Unilateral contact #{}'.format(str(self.__hash__())[-5:]),
                horizontalalignment='center')
        
        if strains is not None:
            if self in strains:
                strain = strains[self][0]
           
            
                if strain < 1e-6:
                    ax.add_patch(Arrow(xplate1, ym,
                     -strain*intensity_factor, 0, linkages_width, color='r'))
                    ax.add_patch(Arrow(xplate2, ym,
                     strain*intensity_factor, 0, linkages_width, color='r'))
                elif strain > -1e-6:
                    ax.add_patch(Arrow(xplate1 - strain*intensity_factor, ym,
                     strain*intensity_factor, 0, linkages_width, color='r'))
                    ax.add_patch(Arrow(xplate2 + strain*intensity_factor, ym,
                     -strain*intensity_factor, 0, linkages_width, color='r'))
                
                if abs(strain) > 1e-6:
                    ax.text(xplate1, ym+0.5*linkages_width,
                            '{} N'.format(round(strain, 3)), color='k',
                            horizontalalignment='center')
        
class ModelSolvingError(Exception):
    def __init__(self, strain_violation):

        # Call the base class constructor with the parameters it needs
        super().__init__('Residue of least squares to high')

#        self.displacement_violation = displacement_violation
        self.strain_violation = strain_violation
        
class ModelConvergenceError(Exception):
    def __init__(self, iters):

        # Call the base class constructor with the parameters it needs
        super().__init__('Number of iteration exceeded: {}'.format(iters))

#        self.displacement_violation = displacement_violation
        self.iters = iters
    
class UnidimensionalModel:
    """
    
    """
    def __init__(self, bodies, linear_linkages, nonlinear_linkages,
                 loads, imposed_displacements):
        self.bodies = bodies
        self.linear_linkages = linear_linkages
        self.nonlinear_linkages = nonlinear_linkages
        self.loads = loads
        self.imposed_displacements = imposed_displacements
#        self.nonlinear_imposed_displacements = nonlinear_imposed_displacements
        
        
    def Graph(self, activated_nonlinar_linkages=None):
        G = nx.Graph()
        G.add_nodes_from(self.bodies)
        if activated_nonlinar_linkages is None:
            activated_nonlinar_linkages = self.nonlinear_linkages[:]
        for linkage in self.linear_linkages:
            G.add_edge(linkage.body1, linkage.body2, linear=True)
        for linkage in activated_nonlinar_linkages:
            G.add_edge(linkage.body1, linkage.body2, linear=False)
        return G
    
    def PlotGraph(self, activated_nonlinar_linkages=None):
        if activated_nonlinar_linkages is None:
            activated_nonlinar_linkages = self.nonlinear_linkages
        plt.figure()
        G = self.Graph(activated_nonlinar_linkages)
        pos={}
        labels = {}
        for body in self.bodies:
            labels[body] = body.name
            pos[body] = (body.initial_position, body.y_position)
            
        linear_edges = []
        nonlinear_edges = []
        for e in G.edges():
            if G.edges()[e]['linear']:
                linear_edges.append(e)
            else:
                nonlinear_edges.append(e)
                
        nx.draw_networkx_nodes(G, pos)
        nx.draw_networkx_edges(G, pos, linear_edges)
        nx.draw_networkx_edges(G, pos, nonlinear_edges, style='--')

        nx.draw_networkx_labels(G, pos, labels)
        
    def IsModelValid(self, activated_nonlinar_linkages):
        """
        Check if each Load is connected to at least a force exit (imposed displacement)
        """
        output_bodies = [i_dis.body for i_dis in self.imposed_displacements]
        input_bodies = [l.body for l in self.loads]
#        print('i', input_bodies)
#        print('o', output_bodies)
       
        G = self.Graph(activated_nonlinar_linkages)
        for ib in input_bodies:
            has_exit = False
            for ob in output_bodies:
                if nx.has_path(G, ib, ob):
                    has_exit = True
#                    print('exit', ib.name, ob.name)
                    break
            if not has_exit:
                return False
        return True
            
        
    def Settings(self, activated_nonlinear_linkages):
        dof = {}
        for ibody, body in enumerate(self.bodies):
            dof[body] = ibody
        ndof = len(dof)
        
        n_constraints = 0
        for linkage in activated_nonlinear_linkages:
            n_constraints += linkage.n_constraints
        for imposed_displacement in self.imposed_displacements:
            n_constraints += imposed_displacement.neqs
        
        return dof, ndof, n_constraints
        
    def LinearSolve(self, activated_nonlinear_linkages=[], strains_tol=1.):
        
        dof, ndof, n_constraints = self.Settings(activated_nonlinear_linkages)
        K = npy.zeros((ndof, ndof))
        F = npy.zeros(ndof)
        D = npy.zeros((n_constraints, ndof))
        Fd = npy.zeros(n_constraints)
#        print(D)
        
        for linkage in self.linear_linkages+activated_nonlinear_linkages:
            dof_linkage1 = dof[linkage.body1]
            dof_linkage2 = dof[linkage.body2]
            Ke, Fe,  = linkage.StiffnessMatrices()
#            print(dof_linkage1, Ke.shape)
            if Ke is not None:
                K[dof_linkage1, dof_linkage1] += Ke[0, 0]
                K[dof_linkage1, dof_linkage2] += Ke[0, 1]
                K[dof_linkage2, dof_linkage1] += Ke[1, 0]
                K[dof_linkage2, dof_linkage2] += Ke[1, 1]

            if Fe is not None:            
                F[dof_linkage1] += Fe[0]
                F[dof_linkage2] += Fe[1]
            
        for load in self.loads:
            F[dof[load.body]] += load.value
        

        i_dof_constraint = ndof
        ic = 0
        for imposed_displacement in self.imposed_displacements:
            D[ic, dof[imposed_displacement.body]] = 1
            Fd[ic] = imposed_displacement.value
            dof[imposed_displacement] = i_dof_constraint
            i_dof_constraint += 1
            ic += 1
        
        for linkage in self.linear_linkages+activated_nonlinear_linkages:
            if linkage.n_constraints > 0:
                Ce, Fce = linkage.ConstraintsMatrices()
                dof_linkage1 = dof[linkage.body1]
                dof_linkage2 = dof[linkage.body2]
                for j in range(linkage.n_constraints):
                    dof[linkage] = i_dof_constraint
                    i_dof_constraint += 1
                    D[ic, dof_linkage1] += Ce[j, 0]
                    D[ic, dof_linkage2] += Ce[j, 1]
                    Fd[ic] += Fce[j]
                    ic += 1
        
        A = npy.vstack((npy.hstack((K, D.T)), 
                        npy.hstack((D, npy.zeros((n_constraints, n_constraints))))))
#        print(K)
#        print(D)
#        print(K.shape, A.shape, linalg.det(A))
#        return A
        b = npy.hstack((F, Fd))
#        print(A, b)
        if abs(linalg.det(A)) < 1e-10:
            res = linalg.lstsq(A, b)
#            print(res)
            x = res[0]
#            print('res', )
            residue = npy.dot(A,x)-b
#            displ_violation = []
            strains_violation = []
            for ires, res in enumerate(npy.abs(residue)):
#                print(ires, res)
                if ires < ndof:
                    if res > strains_tol:
                        strains_violation.append((ires, res))
#                else:
#                    if res > strains_tol:
#                        strains_violation.append((ires, res))
#            print(strains_violation)
            if len(strains_violation)> 0:
                raise ModelSolvingError(strains_violation)
        
        else:
            x = linalg.solve(A, b)
#        print(x)
        displacements = {}
        strains = {}
        for element, dof_element in dof.items():
#            print(element, dof_element)
            if dof_element < ndof:
                displacements[element] = x[dof_element]
            else:
                if element in strains:                    
                    strains[element].append(x[dof_element])
                else:
                    strains[element] = [x[dof_element]]

        return UnidimensionalModelResults(self, activated_nonlinear_linkages, displacements, strains)
    
    def Solve(self, max_iters=100):
        iters = 0
        activated_nonlinear_linkages = []
        new_activated_nonlinear_linkages = self.nonlinear_linkages[:]

        while (set(new_activated_nonlinear_linkages) != set(activated_nonlinear_linkages)) and iters < max_iters:

            activated_nonlinear_linkages = new_activated_nonlinear_linkages[:]
            additional_linkages = []
            while not self.IsModelValid(activated_nonlinear_linkages+additional_linkages):
                
                disactivated_linkages = [l for l in self.nonlinear_linkages if not l in activated_nonlinear_linkages]
                nla = random.choice(range(1, len(disactivated_linkages)))
                additional_linkages = list(random.sample(disactivated_linkages, nla))
            activated_nonlinear_linkages.extend(additional_linkages)
            
            try:
                result = self.LinearSolve(activated_nonlinear_linkages)
            except ModelSolvingError as e:
                raise RuntimeError
            
            new_activated_nonlinear_linkages = []
            disactivated_nonlinear_linkages = []
            for linkage in self.nonlinear_linkages:
                positions = (result.positions[linkage.body1], result.positions[linkage.body2])
                if linkage in result.strains:
                    reactions = result.strains[linkage]
                else:
                    reactions = None
                if linkage in activated_nonlinear_linkages:
                    # Checking if linkage should be desactivated
                    if not linkage.DesactivationCondition(positions, reactions):
                        new_activated_nonlinear_linkages.append(linkage)
                    else:
                        disactivated_nonlinear_linkages.append(linkage)
                else:
                    
                    if linkage.ActivationCondition(positions, reactions):
                        new_activated_nonlinear_linkages.append(linkage)
            iters += 1

        if max_iters == iters:
#            self.Plot()
#            self.PlotGraph()
#            print(self.linear_linkages, self.nonlinear_linkages)
#            print(self.imposed_displacements)
            raise ModelConvergenceError(iters)
            
        result = self.LinearSolve(activated_nonlinear_linkages)
        return result
        
            
    
    def Plot(self, activated_nonlinear_linkages = None, bodies_positions=None,
             strains={}, ax=None, color='k', linkages_width=0.006, intensity_factor=1e-6):
        
        if activated_nonlinear_linkages is None:
            activated_nonlinear_linkages = self.nonlinear_linkages
        
        if bodies_positions is None:
            bodies_positions = {b:b.initial_position for b in self.bodies}
            
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = None
            
        for body in self.bodies:
            ax.plot([bodies_positions[body]], [body.y_position], 'o', color=color)
            
        for linear_linkages in self.linear_linkages+activated_nonlinear_linkages:
            linear_linkages.Plot(bodies_positions, strains=strains, ax=ax, color=color,
                                 linkages_width=linkages_width, intensity_factor=intensity_factor)
            
        for load in self.loads:
            load.Plot(bodies_positions, ax=ax, width=linkages_width,
                      intensity_factor=intensity_factor)
            
        for imposed_displacement in self.imposed_displacements:
            imposed_displacement.Plot(bodies_positions, ax=ax, width=linkages_width,
                      intensity_factor=intensity_factor)
            
            
        ax.autoscale()
        return fig, ax
        
class UnidimensionalModelResults:
    def __init__(self, model, activated_nonlinear_linkages, displacements, strains):
        self.model = model
        self.activated_nonlinear_linkages = activated_nonlinear_linkages
        self.displacements = displacements
        self.strains = strains
        
        self.positions = {b:b.initial_position+displacements[b] for b in self.model.bodies}
        
    def Plot(self, intensity_factor=1e-6):
        fig, ax = plt.subplots()
        self.model.Plot(self.activated_nonlinear_linkages,
                        bodies_positions=self.positions,
                        strains=self.strains,
                        ax=ax, intensity_factor=intensity_factor)
        
        
    