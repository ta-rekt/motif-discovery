# this is for testing implementations of graph algorithms to speed up motif seach in a network.
# after refactoring, testing and verification, make this available as open-source.
# Tarek Tohme, June 2019

import numpy as np
import pandas as pd
import networkx as nx


######################################
#### useful functions and classes ####
######################################

# function f: V -> V that maps nodes in H to nodes in G
# f: ordered list of tuples representing the partial map between D and R, respectively the domain and range of f
class Map():

    def __init__(self, init):   # init is a tuple containing the initial nodes of the partial map
        self.map = np.array([init], dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

    def extend(self, extension):   # extends the partial map by a list of tuples called extension
        temp = np.array(extension, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])
        self.map = np.concatenate((self.map, temp))

    def getDomain(self):
        return self.map['domainNode']

    def getRange(self):
        return self.map['rangeNode']

    def getMap(self):
        return self.map

    def applyMap(self, node):
        for i in self.map['domainNode']:
            if (self.map['domainNode'] == node):
                return self.map['rangeNode'][i]
        return -1


# returns true if g can support h (according to node degree and neighbor degree sequence) and false otherwise.
# g, h: indices of corresponding nodes
# H: query graph
# G: network to be queried
def canSupport(h, g, H, G):

    if(G.degree[g] >= H.degree(h)):
        neighborsH = sortDegrees(H)

        for i in range(len(neighborsH)):
            if(neighborsG[i] < neighborsH[i]):
                # print("rejected because of sequence")
                return False
        return True

    # print("rejected because of degree")
    return False


# returns an ordered list of tuples (node, degree) sorted by degree
# H: network
def sortDegrees(H):
    sortedDegreeH = H.degree()
    sortedDegreeH = [(i, sortedDegreeH(i)) for i in range(len(sortedDegreeH))]
    sortedDegreeH = np.array(sortedDegreeH, dtype=[('node', 'i4'), ('degree', 'i4')])
    sortedDegreeH = np.sort(sortedDegreeH, order='degree')
    return sortedDegreeH


# finds the nodes in H\D with the most neighbors in D and among those, selects the node
# with highest degree and degree sequence
# D: domain of the partial map (list of nodes)
# H: query graph (of type Graph)
def mostConstrainedNode(D, H):
    max = 0
    count = 0
    candidates = np.array([], dtype=[('node', 'i4'), ('neighborsInD', 'i4')])
    neighborsD = {}

    # list of neighbors of D
    for d in D:
        neighborsD = neighborsD.union(set(H[d]))

    # neighbors of D exclusively in H with the number of neighbors in D
    neighborsD = neighborsD.difference(set(D))

    for n in neighborsD:
        for i in H[n]:
            if (i in D):
                count += 1

        count = np.array((n, count), dtype=[('node', 'i4'), ('neighborsInD', 'i4')])
        candidates = np.concatenate(candidates, count)

    # select the candidates with the most neighbors in D
    max = max(candidates['neighborsInD'])
    indices = [i for i, val in enumerate(candidates['neighborsInD']) if val == max]
    indices = [len(H[i]) for i in candidates['node'][indices]] # indices contains degrees of nodes with most neighbors in D

    # select those with highest degree
    max = np.argmax(indices)
    return candidates['node'][max]

    # select those with highest degree sequence (later)


###############################################
#### grochow-kellis motif search algorithm ####
###############################################

# finds all instances of query graph H in network G
# H: query graph
# G: network to be queried
def findSubgraphInstances(H, G):
    instances = np.array([], dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

    # sort nodes of G by degree
    sortedDegreeG = sortDegrees(G)

    # sort by degree sequence. larger sequence is pointwise larger for all nodes (later)

    for g in sortedDegreeG['nodes']:
        for h in H.nodes():
            if(canSupport(h, g, H, G)):
                f = Map((h, g))  # initialize partial map with f(h) = g
                f.extend(isomorphicExtensions(f, H, G))

                if(len(f) != 0):
                    instances = np.concatenate((instances, f.getMap()))
        G.remove_node(g)

    return instances


# finds all isomorphic extensions of a partial map [satisfying the symmetry breaking condition C at h]
# returns them in a list of tuples [(a1, b1), ... ,(ak, bk)]
# f: partial map to be extended
def isomorphicExtensions(f, H, G): #, C, h)
    isomorphisms = np.array([])
    neighborsR = {}
    neighborsD = {}
    D = f.getDomain()
    R = f.getRange()

    if(len(D) == len(H.nodes())):  # can do this because nodes in f are guaranteed to be distinct
        return f

    m = mostConstrainedNode(D, H)

    # list of neighbors of f(D)
    for r in R:
        neighborsR = neighborsR.union(set(G[r]))

    # list of neighbors of D
    for d in D:
        neighborsD = neighborsD.union(set(H[d]))

    # exclusive neighborhood
    neighborsR = neighborsR.difference(set(R))
    neighborsD = neighborsD.difference(set(D))

    while i < len(neighborsR):
        sig = True
        n = neighborsR[i]

        for d in neighborsD:
            if ((d in H[m] & f.applyMap(d) not in G[n]) | (d not in H[m] & f.applyMap(d) in G[n])):
                sig = False
                break

        if(!sig):
            i += 1
        else:
            fp = f
            newNode = np.array((m, n), dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])
            fp.extend(newNode)

            np.append(isomorphisms, isomorphicExtensions(fp, H, G))

    return isomorphisms

# finds symmetry-breaking conditions for H given HE and Aut(H)
# HE: ?
# Aut(H): what is life
def symmetryConditions(HE, Aut)





#######################################
#### onion decomposition algorithm ####
#######################################

# labels each node with its onion layer and coreness
# G: network to be decomposed
def onionDecompose(G):
    D = sortDegrees(G)
    labels = np.array([], dtype=[('node', 'i4'), ('layer', 'i4'), ('coreness', 'i4')])

    core = 1
    layer = 1
    thisLayer = []

    while (len(G.nodes()) > 0):
        for v in G.nodes():
            if (G[v] < core):
                thisLayer.append(v)

        for v in thisLayer:
            newLabel = np.array([(v, core, layer)], dtype=[('node', 'i4'), ('layer', 'i4'), ('coreness', 'i4')])
            labels = np.concatenate(newLabel, labels)

            for w in G[v]:
                D['degree'][w] = D['degree'][w] - 1

            G.remove_node(v)
            # delete from D

        layer = layer + 1

        if (D['degree'][0] > core + 1):
            core = D['degree'][0]

    return labels
