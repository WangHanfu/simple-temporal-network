import copy
import math
import random
import matplotlib.pyplot as plt
import networkx as nx
# =============================
#  FILE:    stn.py
#  AUTHOR:  Sudais Moorad / Muhammad Furrukh Asif
#  DATE:    June 2020
# =============================


# STN is a distance graph with timepoints as nodes and temporal constraints as edges
class STN(nx.DiGraph):
    def __init__(self):
        super().__init__()
        self.add_node(0) # ztp is the zero time point
        self.distance_matrix = {}
        self.dispatchable_apsp_graph = None
        self.dist_up_to_date = False

    def insert_stn_vertex(self, tp):
        self.add_node(tp)
        self.dist_up_to_date = False
    
    # 新增一个边，相当于新增一条temporal constraint
    # tp2-tp1 <= weight ----------  edge: tp1->tp2
    def insert_stn_edge(self, tp1, tp2, w):
        if tp1 in self.nodes and tp2 in self.nodes:
            self.update_distance_matrix() # 内部自动判断是否需要更新距离矩阵
            d_ji = self.distance_matrix[tp2, tp1]
            
            if w < -d_ji: # 新增的约束为inconsistent
                raise ValueError("Inconsistent edge added")
            # elif w >= d_ij: # 新增的约束为redundant
            #     raise ValueError("Redundant edge added")
            elif w == -d_ji: #新增的约束为rigid-connected timepoints
                self.add_edge(tp1, tp2, weight=w)
                self.dist_up_to_date = False
            else: # -d_ji<w<d_ij,新增的约束为non-rigid, consistent
                self.add_edge(tp1, tp2, weight=w)
                self.dist_up_to_date = False
        else:
            self.add_edge(tp1, tp2, weight=w)
            self.dist_up_to_date = False

    def insert_temporal_constraint(self, tp1, tp2, w):
        self.insert_stn_edge(tp1, tp2, weight=w)

    # 若lb=-ub,则tp1,tp2为一个rigid-connected timepoints
    # dij = dji,则tp1,tp2为一个rigid-connected timepoints
    def insert_interval_constraint(self, tp1, tp2, lb,ub):
        self.insert_stn_edge(tp1, tp2, weight=ub)
        self.insert_stn_edge(tp2, tp1, weight=-lb)
        self.dist_up_to_date = False

    def delete_stn_edge(self, tp1, tp2):
        self.remove_edge(tp1, tp2)
        self.dist_up_to_date = False

    def delete_stn_vertex(self, tp):
        self.remove_node(tp)
        self.dist_up_to_date = False

    # 检查一个solution是否满足STN的所有时间约束.遍历所有边，检查每个边上的时间约束是否满足
    def check_solution(self, distances:dict):
        for e in self.edges():
            if distances[e[1]] - distances[e[0]] > self.edges[e]['weight']:
                return False
        return True

    # 更新距离矩阵，使用Floyd-Warshall算法或者Johnson算法。后者更快。
    # equiv., compute APSP graph
    def compute_distance_matrix(self):
        paths = nx.johnson(self,weight='weight') # 使用johnson算法，或者使用Floyd-Warshall算法
        for i in self.nodes:
            for j in self.nodes:
                if i == j:
                    self.distance_matrix[i, j] = 0
                elif j in paths[i].keys():
                    self.distance_matrix[i, j] = nx.path_weight(self,paths[i][j],weight='weight')
                else:
                    self.distance_matrix[i, j] = float('inf')
        
    # TODO 有新的约束添加时，可以用本算法更快的计算出新的距离矩阵
    def incremental_compute_distance_matrix(self, tp1, tp2, weight):
    # 增量更新距离矩阵。2个方法：1, naive O(n^2). 2, constraint propagation.O(n^2).
        pass

    def update_distance_matrix(self):
        if not self.dist_up_to_date:
            self.compute_distance_matrix()
        # self.incremental_compute_distance_matrix()
        self.dist_up_to_date = True

    # Consistency checking: Does an STN have a solution?
    # An STN is consistent if it has no negative cycles; D has non-negative values down its main diagonal
    # method 1: check distance matrix
    # method 2: Directed Path Consistency (DPC) Algorithm. # Dechter,Meiri, andPearl (1991) proposed DPC, for checking whether an STP instance is consistent (i.e. the graph contains no negative cycles).
    def is_consistent(self):
        # 通过检查APSP矩阵（或complete d-graph）对角线上的值是否为0来判断STN是否consistent
        self.update_distance_matrix()
        for i in self.nodes:
            if not self.distance_matrix[i,i] == 0: #对角线上的值为0，则STN是consistent
                return False
        return True
        # TODO DPC算法


    # 输入：shortest_path_array
    def is_consistent_by_cycle(self, shortest_path_array):
        """The STN is consistent if：
         1. it has no negative-length cycles or
         2. it has non-negative values down its main diagonal """

        consistent = True
        for node, val in shortest_path_array.items():
            # Check if the tolerance is too large. Maybe it is better to use
            # only integers and change the resolution to seconds
            if not math.isclose(val[node], 0.0, abs_tol=1e-01):
                consistent = False
        return consistent

    # 利用d-graph，在区间中找到一个可行解，无需担心时间点的求解顺序
    # Decomposability guarantees that such a value can always be found, regardless of the order of assignment
    def generate_a_solution_from_d_graph(self):
        if not self.dist_up_to_date:
            self.compute_distance_matrix()
        sol =  {}
        for i in self.nodes:
            lb = -self.distance_matrix[i,0]
            ub = self.distance_matrix[0,i]
            sol[i] = [lb,ub] # time windows
        return sol

    def get_potential_function(self):
        return dict(nx.single_source_bellman_ford_path_length(self, 0))

    # Tarjan's algorithm for finding strongly connected components
    # 用在dispatchability算法中，用于构建dispatcbility graph的第二个方法。
    def tarjan(self):
        a = nx.strongly_connected_components(self)
        return a
    
    # full path consistency
    def get_dispatchable_apsp_graph(self):
        apsp = copy.deepcopy(self)

        shortest_path_array = nx.floyd_warshall(self)
        if self.is_consistent_by_cycle(shortest_path_array):
            # Get minimal stn by updating the edges of the stn to reflect the shortest path distances
            for i, v in shortest_path_array.items():
                for j,w in v.items():
                    apsp.insert_stn_edge(i, j, w)
            self.dispatchable_apsp_graph = apsp
        else:
            raise("The minimal network is inconsistent. STP could not be solved")

    def visualize(self):
        pos = nx.shell_layout(self)
        nx.draw_networkx_nodes(self, pos, node_size=700)
        nx.draw_networkx_edges(self, pos, edgelist=self.edges, arrowstyle="->", connectionstyle='arc3, rad = 0.1', arrowsize=20, width=3, edge_color='r', alpha=1)
        labels = nx.get_edge_attributes(self, 'weight')
        nx.draw_networkx_edge_labels(self, pos, edge_labels=labels, label_pos=0.3, font_color="b")
        nx.draw_networkx_labels(self, pos, font_size=20, font_family='sans-serif')
        plt.axis('off')
        plt.show()


if __name__ == "__main__":
    network = STN()
    for node_idx in [0,1,2,3]:
        network.insert_stn_vertex(node_idx)

    network.insert_stn_edge(0, 1, 10)
    network.insert_stn_edge(1, 0, -5)
    network.insert_stn_edge(1, 2, 20)
    network.insert_stn_edge(2, 1, -20)
    network.insert_stn_edge(2, 3, 10)
    network.insert_stn_edge(3, 2, -5)
    network.insert_stn_edge(0, 3, 30)
    network.insert_stn_edge(3, 0, -30)
    
    # network.visualize()
    print(network.edges)
    
    network.visualize()
    network.compute_distance_matrix()
    print(network.distance_matrix)
    print(network.is_consistent())
    network.get_dispatchable_apsp_graph()
    print(network.dispatchable_apsp_graph)

