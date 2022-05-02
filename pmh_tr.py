import os
import argparse
import time
import networkx as nx
import gurobipy as gp
from gurobipy import GRB

def process_args():
    parser = argparse.ArgumentParser(description='PMH-TR')

    parser.add_argument('clone_tree', type=str, help='Input clone tree')
    parser.add_argument('leaf_labeling', type=str, help='Input leaf labeling')

    parser.add_argument('-p', '--primary', type=str, help='Primary anatomical site')
    parser.add_argument('-c', '--colormap', metavar='COLORMAP', type=str, help='Color map file', action='store')
    parser.add_argument('--log', action='store_true', default=False, help='Outputs Gurobi logging')
    parser.add_argument('-o', '--output', action='store', default=None, help='Output folder')
    parser.add_argument('-N', '--nsolutions', type=int, help='Maximum number of solutions retained', default=10)
    parser.add_argument('-C', '--count_solutions', action='store_true', default=False, help='Only prints the number of solutions (default=False)')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads')
    parser.add_argument('-s', '--suboptimal', action='store_true', default=False, help='Returns suboptimal solutions without duplicates, \
        may be slow (default=False)')
    parser.add_argument('--print_model', action='store_true', default=False, help='Prints all the ILP constraints in human readable \
        form (experimental, default=False)')

    return parser.parse_args()


class CloneTree:
    
    def __init__(self, clone_tree_filename = None, leaf_labeling_filename = None, edges = None, labels = None, primary_site = None):
        if clone_tree_filename is not None:
            with open(clone_tree_filename, 'r') as clone_tree_file:
                self.tree = nx.DiGraph()
                for edge in clone_tree_file:
                    try:
                        a, b = edge.split()
                    except:
                        raise ValueError('Ill-formatted input clone tree file')
                    self.tree.add_edge(a,b)
                if not nx.is_tree(self.tree):
                    raise ValueError('Input clone tree is not a valid tree')

            self.leaves = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
            self.root = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]
            self.sites = []
            with open(leaf_labeling_filename, 'r') as leaf_labeling_file:
                attrs = dict()
                for label in leaf_labeling_file:
                    try:
                        a, l = label.split()
                    except:
                        raise ValueError('Ill-formatted input leaf labeling file')
                    attrs[a] = {'label': l}
                    if l not in self.sites:
                        self.sites.append(l)
            nx.set_node_attributes(self.tree, attrs)
            for leaf in self.leaves:
                if 'label' not in self.tree.nodes[leaf]:
                    raise ValueError(f'Site label of clone {leaf} is missing')
            if primary_site is not None and primary_site not in self.sites:
                raise ValueError('Primary site not found')
            else:
                self.primary_site = primary_site
            self.n_sites = len(self.sites)
        else:
            self.tree = nx.DiGraph()
            for edge in edges:
                self.tree.add_edge(edge[0], edge[1])
            if not nx.is_tree(self.tree):
                raise RuntimeError('Output is not a tree. Possible bug')
            nx.set_node_attributes(self.tree, labels)
            self.primary_site = primary_site
            self.leaves = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
            self.root = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]

        self.nodes = [v for v in self.tree.nodes]
        self.edges = [e for e in self.tree.edges]

        self.n_nodes = len(self.nodes)
        self.n_edges = len(self.edges)
        self.n_leaves = len(self.leaves)

    def get_label(self, node):
        return self.tree.nodes[node]['label']

    def get_parent_arc(self, node):
        return [k for k in self.tree.in_edges(nbunch=node)][0]

    def get_children_arcs(self, node):
        return [k for k in self.tree.out_edges(nbunch=node)]

    def infer_paths(self):
        self.paths = {}
        self.max_height = 0
        for u in self.leaves:
            self.paths[u] = []
            path_u = nx.shortest_path(self.tree, self.root, u)
            for ii in range(len(path_u) - 1):
                self.paths[u].append( (path_u[ii], path_u[ii+1]) )
            if len(self.paths[u]) > self.max_height:
                self.max_height = len(self.paths[u])

    def write_dot(self, filename, colormap=None):
        with open(filename, 'w+') as f:
            f.write('digraph T {\n\t{\n\t\trank=same\n')
            node_edge_index = {}
            ind = 0
            for l in self.leaves:
                f.write(f'\t\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(l)]},label="{l}\\n{self.get_label(l)}"]\n')
                node_edge_index[l] = ind
                ind += 1
            f.write('\t}\n')
            for v in self.tree.nodes:
                if v not in self.leaves:
                    f.write(f'\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(v)]},label="{v}"]\n')
                    node_edge_index[v] = ind
                    ind += 1
            for i,j in self.tree.edges:
                f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' + 
                    f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\"]\n')
            f.write('}\n')

    def write_tree(self, filename):
        with open(filename, 'w+') as f:
            for edge in self.tree.edges:
                f.write(f'{edge[0]} {edge[1]}\n')

    def write_labeling(self, filename):
        with open(filename, 'w+') as f:
            for node in self.tree.nodes:
                f.write(f'{node} {self.get_label(node)}\n')
        

class ILPSolver:

    def __init__(self, clone_tree, nSolutions, logfile=None, n_threads=None, suboptimal_mode = False):
        self.clone_tree = clone_tree
        self.nSolutions = nSolutions
        self.clone_tree.infer_paths()

        self.node_edge_index = { p:i for i,p in enumerate(clone_tree.nodes + clone_tree.edges) }
        self.site_index = { p:i for i,p in enumerate(clone_tree.sites) }
        if suboptimal_mode:
            self.leaf_index = { p:i for i,p in enumerate(clone_tree.leaves) }

        self.m = gp.Model('ILP')
        self.m._logfile = logfile
        self.m.setParam(GRB.param.LogToConsole, 0)
        self.m.setParam(GRB.param.LogFile, logfile)
        if n_threads is not None:
            self.m.setParam(GRB.Param.Threads, n_threads)

        self.add_ILP_vars()
        if suboptimal_mode:
            self.add_extra_vars_suboptimal()
        self.add_leaf_constraints()
        self.add_polytomy_resolution_compatibility_constraint()
        self.add_original_edges_compatibility_constraint()
        self.add_constraints_for_p()
        if not suboptimal_mode:
            self.add_constraints_for_z()
            self.add_constraints_for_q()
        else:
            self.add_constraints_for_z_b_R_suboptimal()
            self.add_constraints_for_q_suboptimal()
        self.add_constraints_binary()
        if clone_tree.primary_site != None:
            self.add_primary_site_constraints()

        self.set_optimization_function()

        self.m.setParam(GRB.Param.PoolSolutions, nSolutions)
        self.m.setParam(GRB.Param.PoolSearchMode, 2)
        self.m.setParam(GRB.Param.MIPGap, 0)

        start_time = time.time()
        self.m.optimize()
        self.total_time = time.time() - start_time

        self.nActualSolutions = self.m.SolCount
        self.n_migrations, self.n_comigrations, self.n_seeding_sites = self.compute_summary()

    def add_ILP_vars(self):
        self.l = self.m.addMVar((self.clone_tree.n_nodes, self.clone_tree.n_sites), vtype=GRB.BINARY, name="l")
        self.g = self.m.addMVar((self.clone_tree.n_sites, self.clone_tree.n_sites, self.clone_tree.n_edges + self.clone_tree.n_nodes), \
            vtype=GRB.BINARY, name="g")
        self.z = self.m.addMVar((self.clone_tree.n_sites, self.clone_tree.n_sites), vtype=GRB.INTEGER , name="z")
        self.q = self.m.addMVar(self.clone_tree.n_sites, vtype=GRB.BINARY, name="q")
        self.p = self.m.addMVar((self.clone_tree.n_sites, self.clone_tree.n_sites, self.clone_tree.n_edges), vtype=GRB.BINARY, name="p")

    def add_extra_vars_suboptimal(self):
        self.R = self.m.addMVar((self.clone_tree.n_sites, self.clone_tree.n_sites, len(self.clone_tree.leaves) ), vtype=GRB.BINARY, name="R")
        self.b = self.m.addMVar((self.clone_tree.n_sites, self.clone_tree.n_sites, len(self.clone_tree.leaves) ), vtype=GRB.BINARY, name="b")

    def add_leaf_constraints(self):
        for i in range(self.clone_tree.n_nodes):
            if self.clone_tree.nodes[i] in self.clone_tree.leaves:
                hat_ell_i = self.site_index[self.clone_tree.get_label(self.clone_tree.nodes[i])]
                self.m.addConstr( self.l[ i, hat_ell_i ] == 1 )
                sum1 = 0
                for s in range(self.clone_tree.n_sites):
                    if s != hat_ell_i:
                        sum1 += self.l[ i, s ]
                self.m.addConstr(sum1 == 0)

    def add_constraints_for_p(self):
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                for ij in range(self.clone_tree.n_edges):
                    i = self.node_edge_index[self.clone_tree.edges[ij][0]]
                    sum1 = 0
                    sum2 = 0
                    for t_prime in range(self.clone_tree.n_sites):
                        sum1 += self.g[t,t_prime,self.clone_tree.n_nodes + ij]
                        sum2 += self.p[t, t_prime, ij]
                    self.m.addConstr( self.p[s,t,ij] >= self.g[s,t,i] + sum1 - 1)
                    self.m.addConstr( self.p[s,t,ij] >= self.g[s,t,i] + sum2 - 1)
                    self.m.addConstr( self.p[s,t,ij] <= self.g[s,t,i])
                    self.m.addConstr( self.p[s,t,ij] <= sum1 + sum2)

    def add_polytomy_resolution_compatibility_constraint(self):
        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    if s != t:
                        self.m.addConstr( self.g[s,t,i] <= self.l[i,s] )
                        self.m.addConstr( self.g[s,t,i] <= self.l[i,t] )
                    else:
                        self.m.addConstr( self.g[s,t,i] == 0)

        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                if self.clone_tree.nodes[i] != self.clone_tree.root:
                    pii = self.node_edge_index[self.clone_tree.get_parent_arc(self.clone_tree.nodes[i])]
                sum1 = 0
                for t in range(self.clone_tree.n_sites):
                    sum1 += self.g[t,s,i]
                    if self.clone_tree.nodes[i] != self.clone_tree.root:
                        sum1 += self.g[t,s,pii]
                self.m.addConstr(sum1 <= self.l[i,s])

        for i in range(self.clone_tree.n_nodes):
            sum1 = 0
            sum2 = 0
            for s in range(self.clone_tree.n_sites):
                sum2 += self.l[i,s]
                for t in range(self.clone_tree.n_sites):
                    sum1 += self.g[s,t,i]
            self.m.addConstr( sum1 == sum2 - 1 )

        for s in range(self.clone_tree.n_sites):
            for ij in range(self.clone_tree.n_edges):
                sum1 = 0
                sum2 = 0
                for t in range(self.clone_tree.n_sites):
                    sum1 += self.g[s,t,self.clone_tree.n_nodes + ij]
                    sum2 += self.p[s, t, ij]
                self.m.addConstr( sum1 + sum2 <= 1 )

        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                if self.clone_tree.nodes[i] != self.clone_tree.root:
                    pii = self.node_edge_index[self.clone_tree.get_parent_arc(self.clone_tree.nodes[i])]
                delta_i = self.clone_tree.get_children_arcs(self.clone_tree.nodes[i])
                sum1 = 0
                for t in range(self.clone_tree.n_sites):
                    if self.clone_tree.nodes[i] != self.clone_tree.root:
                        sum1 += self.g[ t,s,pii ]
                    for ij in delta_i:
                        sum1 += self.g[s, t, self.node_edge_index[ij]]
                self.m.addConstr( sum1 >= self.l[i, s] )

        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                sum1 = 0
                sum2 = 0
                for t in range(self.clone_tree.n_sites):
                    if s != t:
                        sum1 += self.g[s,t,i] + self.g[t,s,i]
                        sum2 += self.l[i,t]
                self.m.addConstr( self.clone_tree.n_sites * sum1 >=  sum2 + self.clone_tree.n_sites * self.l[i,s] - self.clone_tree.n_sites)

    def add_original_edges_compatibility_constraint(self):
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                for ij in range(self.clone_tree.n_edges):
                    i = self.node_edge_index[self.clone_tree.edges[ij][0]]
                    j = self.node_edge_index[self.clone_tree.edges[ij][1]]
                    self.m.addConstr( self.g[s,t,self.clone_tree.n_nodes + ij] <= self.l[i,s] )
                    self.m.addConstr( self.g[s,t,self.clone_tree.n_nodes + ij] <= self.l[j,t] )

        for ij in range(self.clone_tree.n_edges):
            sum1 = 0
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    sum1 += self.g[s,t,self.clone_tree.n_nodes + ij]
            self.m.addConstr( sum1 == 1 )

    def add_refinement_constraint(self):
        pass

    def add_constraints_for_z(self):
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for k in range(self.clone_tree.n_leaves):
                        sum1 = 0
                        path_k = self.clone_tree.paths[self.clone_tree.leaves[k]]
                        for uv in path_k:
                            sum1 += self.g[ s,t,self.node_edge_index[ uv ] ]
                            sum1 += self.p[ s,t,self.node_edge_index[ uv ] - self.clone_tree.n_nodes ]
                        self.m.addConstr( self.z[s,t] >= sum1 )
                else:
                    self.m.addConstr(self.z[s,t]==0)

    def add_constraints_for_z_b_R_suboptimal(self):
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    sum2 = 0
                    sum3 = 0
                    for k in range(self.clone_tree.n_leaves):
                        sum1 = 0
                        path_k = self.clone_tree.paths[self.clone_tree.leaves[k]]
                        for uv in path_k:
                            sum1 += self.g[ s,t,self.node_edge_index[ uv ] ]
                            sum1 += self.p[ s,t,self.node_edge_index[ uv ] - self.clone_tree.n_nodes ]
                        self.m.addConstr( self.z[s,t] >= sum1 )
                        self.m.addConstr( self.z[s,t] <= sum1 +  self.clone_tree.max_height * (1 - self.b[ s,t,k ]))
                        self.m.addConstr( self.z[ s,t ] - sum1 >= 1 - self.R[ s,t,k ] )
                        self.m.addConstr( self.z[ s,t ] - sum1 <= self.clone_tree.max_height * (1 - self.R[ s,t,k ]) )

                        sum2 += self.b[ s,t,k ]
                        self.m.addConstr( self.b[ s,t,k ] >= self.R[ s,t,k ] - sum3 )
                        sum3 += self.R[ s,t,k ]
                    self.m.addConstr( sum2 == 1 )
                else:
                    self.m.addConstr(self.z[s,t]==0)
                    for k in range(self.clone_tree.n_leaves):
                        self.m.addConstr(self.R[ s,t,k ] == 0)
                        self.m.addConstr(self.b[ s,t,k ] == 0)

    def add_constraints_for_q(self):
        for s in range(self.clone_tree.n_sites):
            self.m.addConstr( self.q[s] <= 1 )
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for i in range(self.clone_tree.n_nodes + self.clone_tree.n_edges):
                        self.m.addConstr( self.q[s] >= self.g[ s,t,i ] )

    def add_constraints_for_q_suboptimal(self):
        for s in range(self.clone_tree.n_sites):
            self.m.addConstr( self.q[s] <= 1 )
            sum1 = 0
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for i in range(self.clone_tree.n_nodes + self.clone_tree.n_edges):
                        self.m.addConstr( self.q[s] >= self.g[ s,t,i ] )
                        sum1 += self.g[ s,t,i ]
            self.m.addConstr( self.q[s] <= sum1)

    def add_constraints_binary(self):
        for i in range(self.clone_tree.n_nodes):
            delta_i = self.clone_tree.get_children_arcs(self.clone_tree.nodes[i])
            if len(delta_i) <= 2:
                sum1 = 0
                for s in range(self.clone_tree.n_sites):
                    sum1 += self.l[i,s]
                self.m.addConstr(sum1  == 1)
            else:
                for s in range(self.clone_tree.n_sites):
                    sum1 = 0
                    for t in range(self.clone_tree.n_sites):
                        sum1 += self.g[s,t,i]
                        for ij in delta_i:
                            sum1 += self.g[s,t,self.node_edge_index[ij]]
                    self.m.addConstr(sum1 >= 2 * self.l[i,s])

    def add_primary_site_constraints(self):
        self.m.addConstr(self.l[self.node_edge_index[clone_tree.root], self.site_index[clone_tree.primary_site]] == 1)
        sum1 = 0
        for s in range(self.clone_tree.n_sites):
            sum1 += self.g[s, self.site_index[clone_tree.primary_site], self.node_edge_index[clone_tree.root]]
        self.m.addConstr(sum1 == 0)

    def set_optimization_function(self):
        sum1 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for k in range(self.clone_tree.n_nodes+self.clone_tree.n_edges):
                        sum1 += self.g[ s,t,k ]
        sum2 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    sum2 += self.z[ s,t ]
        sum3 = 0
        for s in range(self.clone_tree.n_sites):
            sum3 += self.q[s]

        self.m.setObjective( ( (self.clone_tree.n_sites ** 2) * (self.clone_tree.n_sites + 1) ) * sum1 + (self.clone_tree.n_sites + 1) \
            * sum2 + sum3, GRB.MINIMIZE)

    def compute_summary(self, soln_ind=None):
        if soln_ind == None:
            self.m.setParam(GRB.Param.SolutionNumber, 0)
        else:
            self.m.setParam(GRB.Param.SolutionNumber, soln_ind)
        sum1 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for k in range(self.clone_tree.n_nodes + self.clone_tree.n_edges):
                        sum1 += self.g[ s,t,k ].Xn
        sum2 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    sum2 += self.z[ s,t ].Xn
        sum3 = 0
        for s in range(self.clone_tree.n_sites):
            sum3 += self.q[s].Xn
        return round(sum1), round(sum2), round(sum3)

    def count_optimal_solution(self):
        self.m.setParam(GRB.Param.SolutionNumber, 0)
        best_obj_val = self.m.PoolObjVal
        for e in range(self.nActualSolutions):
            self.m.setParam(GRB.Param.SolutionNumber, e)
            if self.m.PoolObjVal - best_obj_val > 0.5:
                return e
        return self.nSolutions

    def print_model(self):
        for i in range(self.clone_tree.n_nodes):
            if self.clone_tree.nodes[i] in self.clone_tree.leaves:
                hat_ell_i = self.site_index[self.clone_tree.get_label(self.clone_tree.nodes[i])]
                print( f'l[ {self.clone_tree.nodes[i]}, {self.clone_tree.sites[hat_ell_i]} ] == 1' )
                conststr = ''
                for s in range(self.clone_tree.n_sites):
                    if s != hat_ell_i:
                        conststr += f'+ l[ {self.clone_tree.nodes[i]}, {self.clone_tree.sites[s]} ]'
                print(conststr + ' == 0')

        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    if s != t:
                        print( f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}] <= '
                            f'l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[s]}]' )
                        print( f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}] <= '
                            f'l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[t]}]' )
                    else:
                        print( f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}] == 0')

        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                conststr = ''
                conststr2 = ''
                for t in range(self.clone_tree.n_sites):
                    if s != t:
                        conststr += f'+ g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}] + ' +\
                            f'g[{self.clone_tree.sites[t]},{self.clone_tree.sites[s]},{self.clone_tree.nodes[i]}]'
                        conststr2 += f'+ l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[t]}]'
                print( f'{self.clone_tree.n_sites} * ({conststr}) >=  ({conststr2}) + {self.clone_tree.n_sites} * '
                    f'l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[s]}] - {self.clone_tree.n_sites}')

        for i in range(self.clone_tree.n_nodes):
            conststr = ''
            conststr2 = ''
            for s in range(self.clone_tree.n_sites):
                conststr2 += f'+ l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[s]}]'
                for t in range(self.clone_tree.n_sites):
                    conststr += f'+ g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}]'
            print( f'{conststr} == {conststr2} - 1' )

        for i in range(self.clone_tree.n_nodes):
            if self.clone_tree.nodes[i] != self.clone_tree.root:
                for s in range(self.clone_tree.n_sites):
                    conststr = ''
                    pii = self.node_edge_index[self.clone_tree.get_parent_arc(self.clone_tree.nodes[i])]
                    for t in range(self.clone_tree.n_sites):
                        conststr += f'+ g[{self.clone_tree.sites[t]},{self.clone_tree.sites[s]},{self.clone_tree.edges[pii-self.clone_tree.n_nodes]}]'
                        conststr += f'+ g[{self.clone_tree.sites[t]},{self.clone_tree.sites[s]},{self.clone_tree.nodes[i]}]'
                    print(f'{conststr} == l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[s]}]')

        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                for ij in range(self.clone_tree.n_edges):
                    i = self.node_edge_index[self.clone_tree.edges[ij][0]]
                    j = self.node_edge_index[self.clone_tree.edges[ij][1]]
                    print( f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}] <= '
                        f'l[{self.clone_tree.nodes[i]},{self.clone_tree.sites[s]}]' )
                    print( f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}] <= '
                        f'l[{self.clone_tree.nodes[j]},{self.clone_tree.sites[t]}]' )
        
        for ij in range(self.clone_tree.n_edges):
            conststr = ''
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    conststr += f'+ g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}]'
            print( f'{conststr} == 1' )

        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                for ij in range(self.clone_tree.n_edges):
                    i = self.node_edge_index[self.clone_tree.edges[ij][0]]
                    conststr = ''
                    conststr2 = ''
                    for t_prime in range(self.clone_tree.n_sites):
                        conststr += f'+ g[{self.clone_tree.sites[t]},{self.clone_tree.sites[t_prime]},{self.clone_tree.edges[ij]}]'
                        conststr2 += f'+ p[{self.clone_tree.sites[t]}, {self.clone_tree.sites[t_prime]}, {self.clone_tree.edges[ij]}]'
                    print( f'p[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}] >= '
                        f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}] + {conststr} - 1')
                    print( f'p[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}] >= '
                        f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}] + {conststr2} - 1')
                    print( f'p[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}] <= '
                        f'g[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]}]')
                    print( f'p[{self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[ij]}] '
                        f'<= {conststr} + {conststr2}')

        for s in range(self.clone_tree.n_sites):
            conststr = ''
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for i in range(self.clone_tree.n_nodes + self.clone_tree.n_edges):
                        if i < self.clone_tree.n_nodes:
                            print( f'q[{self.clone_tree.sites[s]}] >= g[ {self.clone_tree.sites[s]},{self.clone_tree.sites[t]},'
                                f'{self.clone_tree.nodes[i]} ]' )
                            conststr += f'+ g[ {self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.nodes[i]} ]'
                        else:
                            print( f'q[{self.clone_tree.sites[s]}] >= g[ {self.clone_tree.sites[s]},{self.clone_tree.sites[t]},'
                                f'{self.clone_tree.edges[i - self.clone_tree.n_nodes]} ]' )
                            conststr += f'+ g[ {self.clone_tree.sites[s]},{self.clone_tree.sites[t]},{self.clone_tree.edges[i - self.clone_tree.n_nodes]} ]'
            print( f'q[{self.clone_tree.sites[s]}] <= {conststr}')
        for s in range(self.clone_tree.n_sites):
            print( f'q[{self.clone_tree.sites[s]}] <= 1' )
        
        print(f'l[{clone_tree.root}, {clone_tree.primary_site}] == 1')
        conststr = ''
        for s in range(self.clone_tree.n_sites):
            conststr += f'+ g[{self.clone_tree.sites[s]}, {clone_tree.primary_site}, {clone_tree.root}]'
        print(f'{conststr} == 0')
        
    def print_solution(self, e):
        print('soln num - '+ str(e))
        self.m.setParam(GRB.Param.SolutionNumber, e)
        print(self.m.PoolObjVal)
        print('----l----')
        for i in self.clone_tree.nodes:
            for s in self.clone_tree.sites:
                if self.l[self.node_edge_index[i],self.site_index[s]].Xn > 0:
                    print(i, s, self.l[self.node_edge_index[i],self.site_index[s]].Xn)
        print('----g----')
        for s in self.clone_tree.sites:
            for t in self.clone_tree.sites:
                for i in self.clone_tree.nodes:
                    if self.g[self.site_index[s],self.site_index[t],self.node_edge_index[i]].Xn > 0:
                        print(s,t,i,self.g[self.site_index[s],self.site_index[t],self.node_edge_index[i]].Xn)
                for i in clone_tree.edges:
                    if self.g[self.site_index[s],self.site_index[t],self.node_edge_index[i]].Xn > 0:
                        print(s,t,i,self.g[self.site_index[s],self.site_index[t],self.node_edge_index[i]].Xn)
        print('----z----')
        for s in self.clone_tree.sites:
            for t in self.clone_tree.sites:
                if self.z[self.site_index[s],self.site_index[t]].Xn > 0:
                    print(s, t, self.z[self.site_index[s],self.site_index[t]].Xn)
        print('----q----')
        for s in self.clone_tree.sites:
            if self.q[self.site_index[s]].Xn > 0:
                print(s, self.q[self.site_index[s]].Xn)
        print('---p----')
        for s in self.clone_tree.sites:
            for t in self.clone_tree.sites:
                for ij in range(clone_tree.n_edges):
                    if self.p[self.site_index[s],self.site_index[t],ij].Xn > 0:
                        print(s,t,self.clone_tree.edges[ij],self.p[self.site_index[s],self.site_index[t],ij].Xn)

    def refine_tree(self, soln):
        self.m.setParam(GRB.Param.SolutionNumber, soln)
        edges = []
        labels = {}
        node_to_str = {}
        is_refined = [0] * self.clone_tree.n_nodes
        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                if self.l[i,s].Xn > 0.5:
                    is_refined[i] += 1
        for i in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                if self.l[i,s].Xn > 0.5:
                    if is_refined[i] > 1:
                        node_to_str[(i,s)] = f'{self.clone_tree.nodes[i]}^{self.clone_tree.sites[s]}'
                        labels[f'{self.clone_tree.nodes[i]}^{self.clone_tree.sites[s]}'] = {'label': self.clone_tree.sites[s]}
                    else:
                        node_to_str[(i,s)] = f'{self.clone_tree.nodes[i]}'
                        labels[f'{self.clone_tree.nodes[i]}'] = {'label': self.clone_tree.sites[s]}
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                for i in range(self.clone_tree.n_nodes):
                    if self.g[s,t,i].Xn > 0.5:
                        edges.append((node_to_str[(i,s)], node_to_str[(i,t)]))
                for ij in range(self.clone_tree.n_edges):
                    if self.g[s,t,self.clone_tree.n_nodes+ij].Xn > 0.5:
                        i = self.node_edge_index[self.clone_tree.edges[ij][0]]
                        j = self.node_edge_index[self.clone_tree.edges[ij][1]]
                        edges.append((node_to_str[(i,s)], node_to_str[(j,t)]))

        return CloneTree(edges=edges, labels=labels)

    def write_G_dot(self, soln, filename, colormap):
        self.m.setParam(GRB.Param.SolutionNumber, soln)
        with open(filename, 'w+') as f:
            f.write('digraph G {\n')
            for s in range(self.clone_tree.n_sites):
                f.write(f'\t{s} [shape=box,penwidth=3,colorscheme=set19,color={colormap[self.clone_tree.sites[s]]},' + 
                    f'label="{self.clone_tree.sites[s]}"]\n')
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    if s != t:
                        for i in range(self.clone_tree.n_nodes + self.clone_tree.n_edges):
                            if self.g[s,t,i].Xn > 0.5:
                                f.write(f'\t{s} -> {t} [penwidth=3,colorscheme=set19,' +
                                    f'color="{colormap[self.clone_tree.sites[s]]};0.5:{colormap[self.clone_tree.sites[t]]}"]\n')
            f.write('}\n')


def process_colormap(colormap_filename = None, clone_tree = None):
    colormap = {}
    if colormap_filename is not None:
        with open(colormap_filename, 'r') as f:
            for line in f:
                site, color = line.split()
                colormap[site] = color
    else:
        colormap = {s:i for i,s in enumerate(clone_tree.sites)}

    return colormap


if __name__ == '__main__':
    args = process_args()
    clone_tree = CloneTree(clone_tree_filename=args.clone_tree, leaf_labeling_filename=args.leaf_labeling, primary_site=args.primary)

    if args.output is None:
        output_str = '.'
    else:
        output_str = args.output
        if not os.path.exists(output_str):
            os.makedirs(output_str)

    if args.primary is None:
        primary_str = 'ALL'
    else:
        primary_str = args.primary

    if args.log == True:
        logfile = f'{output_str}/{primary_str}-log.txt'
    else:
        logfile = ''

    soln = ILPSolver(clone_tree, args.nsolutions, logfile=logfile, n_threads=args.threads, suboptimal_mode = args.suboptimal)

    if args.print_model:
        soln.print_model()

    if args.count_solutions == True:
        print(f'{primary_str}-\t{int(soln.n_migrations)}\t{int(soln.n_comigrations)}\t{int(soln.n_seeding_sites)}\t\
                {soln.count_optimal_solution()}')
    else:
        colormap = process_colormap(colormap_filename=args.colormap, clone_tree=clone_tree)
        optimal_count = soln.count_optimal_solution()

        if args.suboptimal:
            padding = len(str(soln.nActualSolutions))
        else:
            padding = len(str(optimal_count))

        for e in range(soln.nActualSolutions):
            if args.suboptimal or e < optimal_count:
                T_prime = soln.refine_tree(e)
                primary_str = T_prime.get_label(T_prime.root)
                T_prime.write_dot(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.dot', colormap=colormap)
                T_prime.write_tree(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.tree')
                T_prime.write_labeling(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.labeling')
                soln.write_G_dot(e, f'{output_str}/{primary_str}-G-{str(e).zfill(padding)}.dot', colormap)
                summary = soln.compute_summary(soln_ind=e)
                if e < optimal_count:
                    print(f'{primary_str}-\t{e}\t{summary[0]}\t{summary[1]}\t{summary[2]}\tOptimal\t\t{soln.total_time}')
                else:
                    print(f'{primary_str}-\t{e}\t{summary[0]}\t{summary[1]}\t{summary[2]}\tSuboptimal\t{soln.total_time}')
            else:
                break