from re import A
from xml.dom.minicompat import NodeList
import networkx as nx
import matplotlib.pyplot as plt
from numpy import true_divide
import my_networkx as my_nx
from tkinter.filedialog import askopenfilename
import csv
# import json

damping = 0.85
startingNodeWeight = 1

class Node:
    def __init__(self, name, importance):
        self.preimportance = 0
        self.importance = 0
        self.edges = {}
        self.name = name
        self.importance = importance
    
    def add_edge(self, node, weight):
        self.edges.update({node:weight})

    def add_preimportance(self, num):
        self.preimportance = self.preimportance + num

    def shift_importance(self, numnodes):
        if(self.importance != 0):
            print("Importance leak")
        self.importance = self.preimportance + ((1-damping)*(startingNodeWeight*numnodes)/numnodes)
        self.preimportance = 0

    def distribute_importance_strategy_one(self):
        for key in self.edges:
            key.add_preimportance((self.edges.get(key) * self.importance * damping))
        self.importance = 0

    def is_disconnected(self):
        if(len(self.edges) == 0):
            return True
        return False

    def take_importance(self):
        tmp = self.importance
        self.importance = 0
        return tmp
    
    def name_with_weight(self):
        return self.name + ":" + str(round(self.importance,2))

class Graph:

    def __init__(self):
        self.pagerank_cycles=0
        self.nodelist = []
        print("Graph Initialized")
        self.visual = []

    #for now, all importance initialized to one
    def add_node(self, name):
        self.nodelist.append(Node(name, 1))

    def add_explicit_node(self, node):
        self.nodelist.append(node)

    #from n1 to n2
    def add_edge(self, node1, node2, weight):
        if not (node1 in self.nodelist and node2 in self.nodelist):
            print("Node not found")
            return
        node1.add_edge(node2, weight)
        temp = [node1.name, node2.name]
        self.visual.append(temp)

    def print_graph(self):
        print("PageRank iterations:" + str(self.pagerank_cycles))
        for node in self.nodelist:
            print(node.name + " has importance " + str(node.importance))

    def pagerank_cycle(self):
        self.pagerank_cycles = self.pagerank_cycles + 1
        for node in self.nodelist:
            if(node.is_disconnected()):
                imptospread = node.take_importance() * damping
                imptospread = imptospread / (len(self.nodelist)-1)
                for node2 in self.nodelist:
                    if(node2.name != node.name):
                        node2.preimportance = node2.preimportance + imptospread
            else:
                node.distribute_importance_strategy_one()
        for node in self.nodelist:
            node.shift_importance(len(self.nodelist))
    
    def getWeightList(self):
        wl = []
        for node in self.nodelist:
            wl.append(node.importance)
        return wl

    def pagerankToPrecision(self, precision):
        delta = 99999
        while(delta > precision):
            prewl = self.getWeightList()
            self.pagerank_cycle()
            postwl = self.getWeightList()

            difwl = []
            for i in range(len(prewl)):
                difwl.append(abs(prewl[i] - postwl[i]))

            delta = max(difwl)

    def visualize(self):
        namelist = [node.name_with_weight() for node in self.nodelist]
        longestname = max(namelist, key=len)
        nsize = len(longestname)**2 * 50
        G = nx.DiGraph()
        # for node in self.nodelist:
        for node in self.nodelist:
            for edge in node.edges:
                G.add_edge(node.name_with_weight(), edge.name_with_weight(), weight=node.edges.get(edge))
        pos=nx.circular_layout(G) # pos = nx.nx_agraph.graphviz_layout(G)
        nx.draw_networkx_nodes(G,pos,node_color='b',node_size=nsize,alpha=1)
        nx.draw_networkx_labels(G, pos,font_color='w', font_size=8)
        curved_edges = [edge for edge in G.edges() if reversed(edge) in G.edges()]
        straight_edges = list(set(G.edges()) - set(curved_edges))
        nx.draw_networkx_edges(G, pos,node_size=nsize, edgelist=straight_edges)
        arc_rad = 0.15
        nx.draw_networkx_edges(G, pos,node_size=nsize, edgelist=curved_edges, connectionstyle=f'arc3, rad = {arc_rad}')
        edge_weights = nx.get_edge_attributes(G,'weight')
        curved_edge_labels = {edge: round(edge_weights[edge],2) for edge in curved_edges}
        straight_edge_labels = {edge: round(edge_weights[edge],2) for edge in straight_edges}
        my_nx.my_draw_networkx_edge_labels(G, pos, edge_labels=curved_edge_labels,rotate=False,rad = arc_rad)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=straight_edge_labels,rotate=False)
        plt.show()

#test graph
# def importFile():
#     votedata = []
#     filename = askopenfilename()
#     with open(filename) as file:
#         Lines = file.readlines()
#         for line in Lines:
#             votedata.append(json.loads(line.replace_all("'",'"')))
#     return votedata

# votedata = importFile()
def getSwindex(votematrix, sw):
    rows = len(votematrix)
    for i in range(1, rows):
        if(votematrix[i][0]==sw):
            return i
    return -1

def getTotalAgrees(votematrix, swindex):
    numagrees = 0
    for col in range(1, len(votematrix[0])):
        for row in range(1, len(votematrix)):
            if(row==swindex):
                continue
            if(votematrix[row][col]==votematrix[swindex][col]):
                if(votematrix[row][col]=="NC"):
                    continue
                numagrees = numagrees + 1
    return numagrees

def calcAgreeValue(votematrix, sw1, sw2):
    sw1index = getSwindex(votematrix, sw1)
    sw2index = getSwindex(votematrix, sw2)

    sw1agrees = getTotalAgrees(votematrix, sw1index)
    # sw2agrees = getTotalAgrees(sw2index)

    numagrees = 0
    for i in range(1, len(votematrix)):
        if(votematrix[sw1index][i]==votematrix[sw2index][i]):
            if(votematrix[sw1index][i]=="NC"):
                continue
            numagrees = numagrees + 1
    if(sw1agrees == 0):
        return 0
    return numagrees / sw1agrees
    

def constructGraphFromVotematrix(votematrix):
    #Find the softwares in our vote matrix
    softwares = []
    for i in range(1, len(votematrix)):
        softwares.append(votematrix[i][0])

    #Add the nodes to the graph
    g = Graph()
    swnodes = []
    for sw in softwares:
        node = Node(sw, 1)
        g.add_explicit_node(node)
        swnodes.append(node)
    
    #Add the edges to the graph
    for node1 in swnodes:
        for node2 in swnodes:
            if(node1.name == node2.name):
                continue
            edgeweight = calcAgreeValue(votematrix, node1.name, node2.name)
            #don't add edges with 0 weight
            if(edgeweight > 0.001):
                g.add_edge(node1, node2, edgeweight)
    
    return g

votematrix = [["TopLeft", "M1", "M2", "M3", "M4", "M5"],
["SW1", "C1", "C2", "C3", "NC", "NC"],
["SW2", "C1", "NC", "C3", "NC", "NC"],
["SW3", "NC", "C2", "NC", "NC", "NC"],
["SW4", "C2", "NC", "C1", "C4", "NC"],
["SW5", "C2", "NC", "C1", "C4", "NC"]]


g = constructGraphFromVotematrix(votematrix)

g.pagerankToPrecision(0.01)
print(g.pagerank_cycles)
g.visualize()

